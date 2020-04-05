"""Shared helpers for the prediction methods."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp

from ..folds import Fold
from ..lp import LabelingProblem, Reduced, reduce_problem


def row_normalize(Y: np.ndarray) -> np.ndarray:
    Y = np.asarray(Y, dtype=float)
    s = Y.sum(axis=1, keepdims=True)
    return np.divide(Y, s, out=np.zeros_like(Y), where=s > 0)


def assignment_costs(fold: Fold, nodes: np.ndarray, mode: str = "nonuniform",
                     adj: Optional[sp.csr_matrix] = None) -> np.ndarray:
    """Assignment costs ``c(p, i)`` of the unknown proteins (Section 2.4).

    ``uniform``: 1 for every label.  ``nonuniform``: ``1 / n_pi`` where
    ``n_pi`` is the number of neighbours of ``p`` in the training set with
    function ``i`` (edge-weighted for weighted networks, capped at 2), and 2
    when no neighbour has the function.
    """
    k = fold.k
    if mode == "uniform":
        return np.ones((len(nodes), k))
    if mode != "nonuniform":
        raise ValueError(f"unknown cost mode {mode!r}")
    A = fold.adj if adj is None else adj
    cnt = np.asarray(A[nodes] @ fold.Y_train.astype(float))
    with np.errstate(divide="ignore"):
        c = np.where(cnt > 0, 1.0 / np.maximum(cnt, 1e-12), 2.0)
    return np.minimum(c, 2.0)


def split_reduced(fold: Fold, unary_u: np.ndarray, D: np.ndarray,
                  adj: Optional[sp.csr_matrix] = None) -> Tuple[Reduced, np.ndarray]:
    """Reduced Metric Labeling instance on the unknown proteins ("split" semantics).

    A training protein ``p`` with label set ``B_p`` is replaced by ``|B_p|``
    single-label copies sharing its edges with weight ``w / |B_p|``.  Every
    known node is then fixed, so its edges become the linear terms
    ``w(u,p) * mean_{b in B_p} d(i, b)`` of the unknown neighbour ``u``;
    edges between two known proteins only add a constant and are dropped.
    """
    A = (fold.adj if adj is None else adj).tocsr()
    U = fold.unknown
    known = fold.known
    Yn = row_normalize(fold.Y_train)  # zeros on unknown rows
    lin = unary_u + np.asarray(A[U][:, np.flatnonzero(known)] @ (Yn[known] @ D))
    pos = -np.ones(fold.n, dtype=np.int64)
    pos[U] = np.arange(len(U))
    T = sp.triu(A, k=1).tocoo()
    m = (pos[T.row] >= 0) & (pos[T.col] >= 0)
    red = Reduced(
        n=len(U), k=fold.k, free=np.arange(len(U)), fixed=np.zeros(0, dtype=np.int64),
        fixed_label=-np.ones(len(U), dtype=np.int64), lin=lin,
        fu=pos[T.row[m]], fv=pos[T.col[m]], fw=T.data[m].astype(float), const=0.0,
    )
    return red, U


def choice_reduced(fold: Fold, unary_u: np.ndarray, D: np.ndarray,
                   adj: Optional[sp.csr_matrix] = None) -> Tuple[Reduced, np.ndarray, LabelingProblem]:
    """Reduced instance with the literal semantics of Section 2.1: a training
    protein with label set ``B_p`` has cost 0 on ``B_p`` and infinity elsewhere
    (it may pick any of its known labels).  Returns the reduced instance, the
    dataset index of each free node, and the full problem."""
    A = fold.adj if adj is None else adj
    n, k = fold.n, fold.k
    unary = np.full((n, k), np.inf)
    U = fold.unknown
    unary[U] = unary_u
    Kn = fold.train
    unary[Kn] = np.where(fold.Y_train[Kn], 0.0, np.inf)
    # training proteins whose (possibly corrupted) label set is empty: allow everything
    empty = Kn[~fold.Y_train[Kn].any(axis=1)]
    unary[empty] = 1.0
    T = sp.triu(A, k=1).tocoo()
    P = LabelingProblem(unary, T.row.astype(np.int64), T.col.astype(np.int64), T.data.astype(float), D)
    red = reduce_problem(P)
    return red, red.free, P


def local_scores(fold: Fold, nodes: np.ndarray, unary_u: np.ndarray, marg: np.ndarray, D: np.ndarray,
                 temperature: Optional[float] = None, adj: Optional[sp.csr_matrix] = None) -> np.ndarray:
    """Local min-marginal scores derived from the Metric Labeling objective.

    For node ``u`` and label ``i``::

        E_u(i) = c(u, i) + sum_v w(u, v) sum_j m_v(j) d(i, j)

    is the value of ``u``'s terms of ``Q`` (Eq. 5) when ``u`` alone takes label
    ``i`` while its neighbours stay at their (fractional) solution ``m_v``.
    The score is the regret relative to the best label, scaled by the largest
    possible regret::

        s(u, i) = 1 - (E_u(i) - min_j E_u(j)) / (max c + W_u max d)

    (``W_u`` = weighted degree), so the best label scores 1 and scores are
    comparable across proteins.  If ``temperature`` is given, the exponential
    form ``exp(-(E_u(i) - min_j E_u(j)) / (temperature (1 + W_u)))`` is used.
    """
    A = (fold.adj if adj is None else adj).tocsr()
    E = unary_u + np.asarray(A[nodes] @ (marg @ D))
    W = np.asarray(A[nodes].sum(axis=1)).ravel()
    E = E - E.min(axis=1, keepdims=True)
    if temperature:
        return np.exp(-E / (temperature * (1.0 + W))[:, None])
    cmax = float(np.max(unary_u)) if unary_u.size else 0.0
    scale = cmax + W * float(np.max(D))
    return np.clip(1.0 - E / np.maximum(scale, 1e-12)[:, None], 0.0, 1.0)


def perturb_weights(w: np.ndarray, rng: np.random.Generator, eps: float = 1e-5) -> np.ndarray:
    """Uniform perturbation in (-w_max eps, w_max eps), as in Nabieva et al."""
    if len(w) == 0:
        return w
    wm = float(np.max(w))
    return np.maximum(w + rng.uniform(-wm * eps, wm * eps, size=len(w)), 0.0)


def perturbed_adj(adj: sp.csr_matrix, rng: np.random.Generator, eps: float = 1e-5) -> sp.csr_matrix:
    T = sp.triu(adj, k=1).tocoo()
    w = perturb_weights(T.data.astype(float), rng, eps)
    n = adj.shape[0]
    return sp.csr_matrix((np.r_[w, w], (np.r_[T.row, T.col], np.r_[T.col, T.row])), shape=(n, n))
