"""Graph-cut move-making algorithms.

* :func:`binary_cut` - exact minimisation of a submodular binary energy via
  s-t min cut (PyMaxflow / Boykov-Kolmogorov; SciPy fallback).
* :func:`alpha_expansion` - Boykov, Veksler & Zabih (2001) for metric
  energies; used for the GenMultiCut (Potts) baseline.
* :func:`hierarchical_moves` - the hierarchical graph cuts of Kumar & Koller
  (UAI 2009) for r-HST metrics, the core of the semimetric MAP baseline.

All functions work on a :class:`~metriclabeling.lp.Reduced` instance (free
nodes only, fixed neighbours folded into the linear costs).
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .hst import HST
from .lp import Reduced

try:  # pragma: no cover - import guard
    import maxflow as _maxflow

    HAVE_MAXFLOW = True
except Exception:  # pragma: no cover
    HAVE_MAXFLOW = False


def _big(lin: np.ndarray) -> float:
    fin = lin[np.isfinite(lin)]
    return 1e3 * (1.0 + (np.abs(fin).max() if fin.size else 1.0))


def binary_cut(U0: np.ndarray, U1: np.ndarray, eu: np.ndarray, ev: np.ndarray,
               A: np.ndarray, B: np.ndarray, C: np.ndarray, Dd: np.ndarray) -> np.ndarray:
    """Minimise ``sum_u U_u(b_u) + sum_e V_e(b_u, b_v)`` over ``b in {0,1}^n``.

    ``V_e(0,0)=A, V_e(0,1)=B, V_e(1,0)=C, V_e(1,1)=Dd``; terms must be
    submodular (``B + C >= A + Dd``); tiny negative excesses are truncated.
    """
    n = len(U0)
    u0 = U0.astype(float).copy()
    u1 = U1.astype(float).copy()
    # E = A + (C-A) b_u + (Dd-C) b_v + (B+C-A-Dd)(1-b_u) b_v
    np.add.at(u1, eu, C - A)
    np.add.at(u1, ev, Dd - C)
    lam = np.maximum(B + C - A - Dd, 0.0)
    shift = np.minimum(u0, u1)
    u0 -= shift
    u1 -= shift
    if HAVE_MAXFLOW:
        g = _maxflow.Graph[float](n, len(eu))
        nodes = g.add_nodes(n)
        g.add_grid_tedges(np.asarray(nodes), u1, u0)  # source cap = cost of b=1, sink cap = cost of b=0
        m = lam > 0
        if m.any():
            nodes = np.asarray(nodes)
            g.add_edges(nodes[eu[m]], nodes[ev[m]], lam[m], np.zeros(int(m.sum())))
        g.maxflow()
        return np.asarray(g.get_grid_segments(np.asarray(nodes)), dtype=bool)
    return _binary_cut_scipy(u0, u1, eu, ev, lam)


def _binary_cut_scipy(u0, u1, eu, ev, lam):  # pragma: no cover - fallback path
    import scipy.sparse as sp
    from scipy.sparse.csgraph import maximum_flow

    n = len(u0)
    s, t = n, n + 1
    scale = 1e6 / max(1e-12, max(u0.max(initial=0), u1.max(initial=0), lam.max(initial=0)))
    rows = np.r_[np.full(n, s), np.arange(n), eu]
    cols = np.r_[np.arange(n), np.full(n, t), ev]
    caps = np.r_[u1, u0, lam] * scale
    caps = np.round(caps).astype(np.int64)
    G = sp.csr_matrix((caps, (rows, cols)), shape=(n + 2, n + 2))
    G.sum_duplicates()
    res = maximum_flow(G.astype(np.int32), s, t)
    R = (G - res.flow).tocsr()
    R.data = np.maximum(R.data, 0)
    R.eliminate_zeros()
    # nodes reachable from s in the residual graph are on the source side (b = 0)
    from scipy.sparse.csgraph import breadth_first_order

    reach = breadth_first_order(R, s, directed=True, return_predecessors=False)
    b = np.ones(n, dtype=bool)
    b[reach[reach < n]] = False
    return b


def reduced_energy(red: Reduced, labels: np.ndarray, D: np.ndarray) -> float:
    lin = red.lin[np.arange(red.n_free), labels]
    return float(lin.sum() + np.sum(red.fw * D[labels[red.fu], labels[red.fv]]) + red.const)


def initial_labels(red: Reduced) -> np.ndarray:
    lin = np.where(np.isfinite(red.lin), red.lin, np.inf)
    return np.argmin(lin, axis=1)


def alpha_expansion(red: Reduced, D: np.ndarray, labels: Optional[np.ndarray] = None,
                    max_sweeps: int = 10, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """alpha-expansion (Boykov et al. 2001) on the free nodes of ``red``."""
    rng = np.random.default_rng(rng)
    n, k = red.n_free, red.k
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    big = _big(red.lin)
    lin = np.where(np.isfinite(red.lin), red.lin, big)
    f = initial_labels(red) if labels is None else labels.copy()
    best = reduced_energy(red, f, D)
    idx = np.arange(n)
    for _ in range(max_sweeps):
        improved = False
        for a in rng.permutation(k):
            U0 = lin[idx, f]
            U1 = lin[:, a]
            fu, fv = f[red.fu], f[red.fv]
            w = red.fw
            A = w * D[fu, fv]
            B = w * D[fu, a]
            C = w * D[a, fv]
            Dd = np.zeros_like(A)
            b = binary_cut(U0, U1, red.fu, red.fv, A, B, C, Dd)
            g = np.where(b, a, f)
            e = reduced_energy(red, g, D)
            if e < best - 1e-12:
                best, f, improved = e, g, True
        if not improved:
            break
    return f


def hierarchical_moves(red: Reduced, tree: HST) -> np.ndarray:
    """Kumar & Koller (2009) hierarchical graph cuts for an r-HST metric.

    Bottom-up over the tree: a leaf's labeling assigns its label to every
    node; an internal node fuses its children's labelings one child at a time
    with a binary (keep / switch) move solved exactly by a min cut.  On an
    HST every such move is submodular.
    """
    n = red.n_free
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    Dt = tree.label_distance()
    big = _big(red.lin)
    lin = np.where(np.isfinite(red.lin), red.lin, big)
    idx = np.arange(n)
    lab = {}
    for t in range(tree.n_nodes - 1, -1, -1):
        ch = tree.children[t]
        if not ch:
            lab[t] = np.full(n, tree.label_of_leaf[t], dtype=np.int64)
            continue
        f = lab.pop(ch[0])
        for c in ch[1:]:
            g = lab.pop(c)
            U0 = lin[idx, f]
            U1 = lin[idx, g]
            w = red.fw
            A = w * Dt[f[red.fu], f[red.fv]]
            B = w * Dt[f[red.fu], g[red.fv]]
            C = w * Dt[g[red.fu], f[red.fv]]
            Dd = w * Dt[g[red.fu], g[red.fv]]
            b = binary_cut(U0, U1, red.fu, red.fv, A, B, C, Dd)
            f = np.where(b, g, f)
        lab[t] = f
    return lab[0]
