"""Distances between GO terms (Section 2.4 of the paper).

* :func:`d_sp`   shortest-path distance in the GO DAG / GO diameter (metric)
* :func:`d_lca`  (b+c)/(2a+b+c) with the lowest common ancestor (semimetric)
* :func:`d_lin`  reciprocal of Lin's similarity (semimetric)
* :func:`d_kb`   average network distance between training proteins (semimetric)
* :func:`combine` (1-alpha) d + alpha d_KB   (Eq. 16)

All functions return dense symmetric ``(k, k)`` matrices with a zero diagonal.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import shortest_path

from .dataset import GOInfo, graph_diameter


def _finalize(D: np.ndarray) -> np.ndarray:
    D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)
    return D


def normalize_max(D: np.ndarray) -> np.ndarray:
    m = float(np.max(D))
    return D / m if m > 0 else D.copy()


# --------------------------------------------------------------------------- #
def d_sp(go: GOInfo) -> np.ndarray:
    """Shortest-path distance in the (undirected) GO DAG divided by its diameter."""
    D = np.array(go.sp_len, dtype=float)
    if not np.all(np.isfinite(D)):
        D[~np.isfinite(D)] = go.go_diameter
    return _finalize(D / go.go_diameter)


def d_lca(go: GOInfo) -> np.ndarray:
    """LCA distance ``(b + c) / (2a + b + c)``.

    ``a`` is the depth (shortest root path) of the lowest common ancestor
    ``u``, ``b``/``c`` the shortest upward distances from ``x``/``y`` to ``u``.
    In a DAG several common ancestors exist; we take the one giving the
    smallest distance (i.e. the largest Wu-Palmer similarity ``2a/(2a+b+c)``).
    """
    up = go.up  # (k, |A|)
    depth = go.anc_depth[None, :]
    k = up.shape[0]
    D = np.ones((k, k))
    for i in range(k):
        b = up[i][None, :]  # (1, |A|)
        bc = b + up  # (k, |A|), inf where not common
        den = 2 * depth + bc
        with np.errstate(invalid="ignore", divide="ignore"):
            val = np.where(np.isfinite(bc) & (den > 0), bc / den, np.inf)
        D[i] = val.min(axis=1)
    D[~np.isfinite(D)] = 1.0
    return _finalize(D)


def term_probabilities(go: GOInfo, train: np.ndarray, pseudocount: float = 1.0) -> np.ndarray:
    """Empirical probability that a training protein is annotated to each ancestor.

    ``train`` indexes the proteins whose annotations may be used.  Only
    proteins with at least one BP annotation enter the denominator.  A
    pseudocount avoids ``log(0)`` for terms unseen in a training fold.
    """
    train = np.asarray(train)
    train = train[go.has_bp[train]]
    counts = np.asarray(go.P[train].sum(axis=0)).ravel().astype(float)
    N = float(len(train))
    return (counts + pseudocount) / (N + pseudocount)


def d_lin(go: GOInfo, prob: np.ndarray, mode: str = "reciprocal") -> np.ndarray:
    """Distance from Lin's similarity ``2 log Pr(lca) / (log Pr(x) + log Pr(y))``.

    The paper takes the reciprocal ``(log Pr(x) + log Pr(y)) / (2 log Pr(lca))``
    and we set ``d(x,x) = 0``.  ``lca`` is the most informative common
    ancestor (smallest ``Pr``).  When that ancestor is the root
    (``Pr = 1``) the reciprocal is infinite; such entries are set to the
    largest finite value.  ``mode="one_minus"`` gives ``1 - sim`` instead.
    The result is scaled to a maximum of 1.
    """
    is_anc = go.is_ancestor  # (k, |A|)
    logp = np.log(np.clip(prob, 1e-300, 1.0))  # (|A|,)
    term_idx = np.array([go.anc_ids.index(t) for t in go.term_ids])
    lx = logp[term_idx]  # (k,)
    k = is_anc.shape[0]
    # most informative common ancestor: min log Pr over common ancestors
    S = np.zeros((k, k))
    for i in range(k):
        common = is_anc[i][None, :] & is_anc  # (k, |A|)
        mica = np.where(common, logp[None, :], 0.0).min(axis=1)  # <= 0
        den = lx[i] + lx
        with np.errstate(invalid="ignore", divide="ignore"):
            S[i] = np.where(den < 0, 2.0 * mica / den, 0.0)
    S = np.clip(S, 0.0, 1.0)
    np.fill_diagonal(S, 1.0)
    if mode == "one_minus":
        D = 1.0 - S
    else:
        with np.errstate(divide="ignore"):
            D = 1.0 / S
        finite = np.isfinite(D)
        cap = D[finite].max() if finite.any() else 1.0
        D[~finite] = cap
    return normalize_max(_finalize(D))


# --------------------------------------------------------------------------- #
def network_distances(adj: sp.csr_matrix, sources: np.ndarray, chunk: int = 1024) -> np.ndarray:
    """Unweighted BFS distances from ``sources`` to all nodes (float32)."""
    out = np.empty((len(sources), adj.shape[0]), dtype=np.float32)
    for s in range(0, len(sources), chunk):
        idx = sources[s : s + chunk]
        out[s : s + len(idx)] = shortest_path(adj, unweighted=True, directed=False, indices=idx)
    return out


def d_kb(
    adj: sp.csr_matrix,
    Y_train: np.ndarray,
    diameter: Optional[float] = None,
    sp_cache: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> np.ndarray:
    """Knowledge-based distance from training annotations.

    ``d_KB(x, y) = sum_{p in P_x, q in P_y} sp(p, q) / (diam * |P_x| |P_y|)``
    where ``P_x`` are training proteins annotated with ``x``.  Terms without
    training proteins get distance 1 to everything.

    ``sp_cache`` may hold ``(nodes, distances)`` from :func:`network_distances`
    computed for a superset of the training proteins.
    """
    Y = np.asarray(Y_train, dtype=float)
    labeled = np.flatnonzero(Y.any(axis=1))
    if diameter is None:
        diameter = graph_diameter(adj)
    if sp_cache is not None:
        nodes, Dall = sp_cache
        pos = {int(p): i for i, p in enumerate(nodes)}
        rows = np.array([pos[int(p)] for p in labeled])
        Dsub = Dall[rows][:, labeled].astype(float)
    else:
        Dsub = network_distances(adj, labeled)[:, labeled].astype(float)
    Dsub[~np.isfinite(Dsub)] = diameter
    Yl = Y[labeled]  # (L, k)
    sizes = Yl.sum(axis=0)
    tot = Yl.T @ Dsub @ Yl  # (k, k)
    with np.errstate(invalid="ignore", divide="ignore"):
        D = tot / (diameter * np.outer(sizes, sizes))
    D[~np.isfinite(D)] = 1.0
    return _finalize(D)


def combine(d: np.ndarray, d_kb_: np.ndarray, alpha: float) -> np.ndarray:
    """Eq. (16): ``(1 - alpha) d + alpha d_KB``."""
    return _finalize((1.0 - alpha) * d + alpha * d_kb_)


def uniform_metric(k: int) -> np.ndarray:
    D = np.ones((k, k))
    np.fill_diagonal(D, 0.0)
    return D


# --------------------------------------------------------------------------- #
def triangle_violations(D: np.ndarray, tol: float = 1e-9) -> Tuple[int, float]:
    """Number of violated triangle inequalities and the largest violation.

    Counts ordered triples ``(i, j, l)`` (distinct) with
    ``D[i, j] > D[i, l] + D[l, j] + tol``.
    """
    D = np.asarray(D, dtype=float)
    k = D.shape[0]
    count = 0
    worst = 0.0
    for l in range(k):
        V = D - (D[:, l][:, None] + D[l, :][None, :])
        V[l, :] = -np.inf
        V[:, l] = -np.inf
        np.fill_diagonal(V, -np.inf)
        count += int((V > tol).sum())
        worst = max(worst, float(V.max()))
    return count // 2, max(worst, 0.0)


def is_metric(D: np.ndarray, tol: float = 1e-9) -> bool:
    D = np.asarray(D, dtype=float)
    if np.any(D < -tol) or not np.allclose(D, D.T, atol=tol) or np.any(np.abs(np.diag(D)) > tol):
        return False
    return triangle_violations(D, tol)[0] == 0


def distance_from_spec(
    name: str,
    go: GOInfo,
    adj: Optional[sp.csr_matrix] = None,
    Y_train: Optional[np.ndarray] = None,
    train_idx: Optional[np.ndarray] = None,
    alpha: float = 0.0,
    net_diameter: Optional[float] = None,
    sp_cache=None,
    lin_mode: str = "reciprocal",
) -> np.ndarray:
    """Build one of the paper's distances by name (``sp``/``lca``/``lin``/``kb``/``uniform``),
    optionally combined with ``d_KB`` (Eq. 16) using weight ``alpha``."""
    name = name.lower()
    k = len(go.term_ids)
    if name in ("sp", "shortest", "shortest_path"):
        base = d_sp(go)
    elif name == "lca":
        base = d_lca(go)
    elif name == "lin":
        if train_idx is None:
            raise ValueError("d_Lin needs the training proteins")
        base = d_lin(go, term_probabilities(go, train_idx), mode=lin_mode)
    elif name == "kb":
        base = None
        alpha = 1.0
    elif name in ("uniform", "potts"):
        return uniform_metric(k)
    else:
        raise ValueError(f"unknown distance {name!r}")
    if alpha <= 0:
        return base
    if adj is None or Y_train is None:
        raise ValueError("d_KB needs the network and the training labels")
    dkb = d_kb(adj, Y_train, diameter=net_diameter, sp_cache=sp_cache)
    if base is None:
        return dkb
    return combine(base, dkb, alpha)
