"""Hierarchically well-separated trees (HSTs).

* :func:`frt_tree` - random 2-HST of Fakcharoenphol, Rao & Talwar (2003):
  dominates the input metric and stretches distances by ``O(log k)`` in
  expectation.  Used by the Chekuri et al. / Kleinberg-Tardos approximation
  for Metric Labeling and by the semimetric MAP baseline of Kumar & Koller.
* :func:`star_tree` - the uniform (Potts) metric as a tree.
* :func:`best_tree` - derandomisation by sampling: keep the sampled tree with
  the smallest cost under a user-supplied criterion.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional

import numpy as np


@dataclass
class HST:
    """Rooted tree whose leaves are the ``k`` labels.

    ``parent[t]`` (``-1`` for the root) and ``length[t]`` (edge to parent)
    are indexed by tree node; nodes are stored in topological order
    (parents before children).  ``leaf[i]`` is the node of label ``i``.
    """

    parent: np.ndarray
    length: np.ndarray
    leaf: np.ndarray

    def __post_init__(self) -> None:
        self.parent = np.asarray(self.parent, dtype=np.int64)
        self.length = np.asarray(self.length, dtype=np.float64)
        self.leaf = np.asarray(self.leaf, dtype=np.int64)
        n = len(self.parent)
        self.children: List[List[int]] = [[] for _ in range(n)]
        for t in range(n):
            p = self.parent[t]
            if p >= 0:
                if p >= t:
                    raise ValueError("HST nodes must be in topological order")
                self.children[p].append(t)
        self.label_of_leaf = -np.ones(n, dtype=np.int64)
        self.label_of_leaf[self.leaf] = np.arange(len(self.leaf))
        # height (distance from root)
        h = np.zeros(n)
        for t in range(1, n):
            h[t] = h[self.parent[t]] + self.length[t]
        self.height = h

    @property
    def n_nodes(self) -> int:
        return len(self.parent)

    @property
    def k(self) -> int:
        return len(self.leaf)

    @property
    def root(self) -> int:
        return 0

    def incidence(self) -> np.ndarray:
        """``(k, n_nodes)`` 0/1 matrix: label ``i`` lies in the subtree of node ``t``."""
        M = np.zeros((self.k, self.n_nodes))
        for i, t in enumerate(self.leaf):
            u = int(t)
            while u >= 0:
                M[i, u] = 1.0
                u = int(self.parent[u])
        return M

    def label_distance(self) -> np.ndarray:
        """Tree distance between all pairs of labels."""
        M = self.incidence()  # (k, n)
        # d(i,j) = sum over non-root nodes t of length[t] * [exactly one of i,j under t]
        L = self.length.copy()
        L[0] = 0.0
        A = M * L[None, :]
        common = A @ M.T  # sum of lengths of shared ancestor edges
        hi = A.sum(axis=1)
        D = hi[:, None] + hi[None, :] - 2 * common
        np.fill_diagonal(D, 0.0)
        return np.maximum(D, 0.0)

    def arcs(self):
        """Directed arcs of the tree (both directions) with lengths: (src, dst, len)."""
        t = np.arange(1, self.n_nodes)
        p = self.parent[1:]
        L = self.length[1:]
        return np.r_[t, p], np.r_[p, t], np.r_[L, L]


def _compress(parent: List[int], length: List[float], leaf: np.ndarray) -> HST:
    """Remove internal nodes with a single child and re-index topologically."""
    n = len(parent)
    children: List[List[int]] = [[] for _ in range(n)]
    for t, p in enumerate(parent):
        if p >= 0:
            children[p].append(t)
    is_leaf = np.zeros(n, dtype=bool)
    is_leaf[leaf] = True
    parent = list(parent)
    length = list(length)
    # bottom-up: splice out unary internal nodes
    root = 0
    order = _topo(parent)
    for t in reversed(order):
        if is_leaf[t]:
            continue
        if len(children[t]) == 1:
            c = children[t][0]
            p = parent[t]
            if p >= 0:
                parent[c] = p
                length[c] += length[t]
                children[p] = [c if x == t else x for x in children[p]]
            else:  # root with one child: child becomes root
                parent[c] = -1
                length[c] = 0.0
                root = c
            children[t] = []
            parent[t] = -2  # removed
    # re-index by BFS from root
    new_id = {}
    order = [root]
    q = 0
    while q < len(order):
        u = order[q]
        q += 1
        for c in children[u]:
            if parent[c] != -2:
                order.append(c)
    for i, t in enumerate(order):
        new_id[t] = i
    par = np.array([new_id[parent[t]] if parent[t] >= 0 else -1 for t in order], dtype=np.int64)
    ln = np.array([length[t] if parent[t] >= 0 else 0.0 for t in order])
    lf = np.array([new_id[int(t)] for t in leaf], dtype=np.int64)
    return HST(par, ln, lf)


def _topo(parent: List[int]) -> List[int]:
    n = len(parent)
    children: List[List[int]] = [[] for _ in range(n)]
    roots = []
    for t, p in enumerate(parent):
        if p >= 0:
            children[p].append(t)
        elif p == -1:
            roots.append(t)
    order = list(roots)
    q = 0
    while q < len(order):
        for c in children[order[q]]:
            order.append(c)
        q += 1
    return order


def frt_tree(D: np.ndarray, rng: Optional[np.random.Generator] = None, eps: float = 1e-6) -> HST:
    """Sample an FRT 2-HST for the (pseudo-)metric ``D``.

    Points at distance 0 are separated at the bottom by edges of length
    ``eps * min positive distance``.  Edge lengths are in the units of ``D``.
    """
    rng = np.random.default_rng(rng)
    D = np.asarray(D, dtype=float)
    k = D.shape[0]
    if k == 1:
        return HST(np.array([-1]), np.array([0.0]), np.array([0]))
    pos = D[D > 0]
    scale = float(pos.min()) if pos.size else 1.0
    Dn = D / scale
    delta = int(np.ceil(np.log2(max(Dn.max(), 1.0)))) + 1
    beta = rng.uniform(1.0, 2.0)
    perm = rng.permutation(k)
    Dp = Dn[:, perm]  # columns in permutation order

    parent: List[int] = [-1]
    length: List[float] = [0.0]
    current = [(0, np.arange(k))]
    for i in range(delta - 1, -1, -1):
        r = beta * 2.0 ** (i - 1)
        nxt = []
        for node, members in current:
            within = Dp[members] <= r
            first = within.argmax(axis=1)
            for c in np.unique(first):
                sub = members[first == c]
                parent.append(node)
                length.append(2.0 ** (i + 1))
                nxt.append((len(parent) - 1, sub))
        current = nxt
    leaf = np.empty(k, dtype=np.int64)
    for node, members in current:
        if len(members) == 1:
            leaf[members[0]] = node
        else:
            for x in members:
                parent.append(node)
                length.append(eps)
                leaf[x] = len(parent) - 1
    hst = _compress(parent, length, leaf)
    hst.length = hst.length * scale
    hst.__post_init__()
    return hst


def star_tree(k: int, edge: float = 0.5) -> HST:
    """Star with ``k`` leaves; leaf-to-leaf distance ``2 * edge`` (uniform metric)."""
    parent = np.r_[-1, np.zeros(k, dtype=np.int64)]
    length = np.r_[0.0, np.full(k, edge)]
    return HST(parent, length, np.arange(1, k + 1))


def tree_from_metric_if_uniform(D: np.ndarray, tol: float = 1e-12) -> Optional[HST]:
    off = D[~np.eye(D.shape[0], dtype=bool)]
    if off.size and np.ptp(off) <= tol:
        return star_tree(D.shape[0], off[0] / 2.0)
    return None


def stretch_cost(D: np.ndarray, weights: Optional[np.ndarray] = None) -> Callable[[HST], float]:
    """Criterion: weighted average stretch sum_{i<j} W_ij d_T(i,j) / D(i,j)."""
    k = D.shape[0]
    iu = np.triu_indices(k, 1)
    d = D[iu]
    w = np.ones_like(d) if weights is None else np.asarray(weights)[iu]
    ok = d > 0

    def cost(T: HST) -> float:
        dt = T.label_distance()[iu]
        return float(np.sum(w[ok] * dt[ok] / d[ok]))

    return cost


def lsq_cost(D: np.ndarray) -> Callable[[HST], float]:
    """Criterion: residual of the best least-squares fit ``min_s sum (s d_T - D)^2``."""
    iu = np.triu_indices(D.shape[0], 1)
    d = D[iu]

    def cost(T: HST) -> float:
        dt = T.label_distance()[iu]
        den = float(dt @ dt)
        if den <= 0:
            return float(d @ d)
        return float(d @ d - (d @ dt) ** 2 / den)

    return cost


def rescale_to_fit(T: HST, D: np.ndarray) -> HST:
    """Scale all edge lengths by the least-squares factor ``s = <D, d_T> / <d_T, d_T>``.

    FRT trees dominate ``D`` and therefore overstate distances by a constant
    factor on average; a global rescaling keeps the tree structure (and the
    ``O(log k)`` distortion) but restores the balance between assignment
    and separation costs in the Metric Labeling objective.
    """
    iu = np.triu_indices(D.shape[0], 1)
    dt = T.label_distance()[iu]
    den = float(dt @ dt)
    if den <= 0:
        return T
    s = float(D[iu] @ dt) / den
    if s <= 0:
        return T
    return HST(T.parent.copy(), T.length * s, T.leaf.copy())


def sample_trees(D: np.ndarray, n: int, seed=None) -> List[HST]:
    rng = np.random.default_rng(seed)
    uni = tree_from_metric_if_uniform(D)
    if uni is not None:
        return [uni]
    return [frt_tree(D, rng) for _ in range(n)]


def best_tree(D: np.ndarray, n_samples: int = 32, seed=None, criterion: Optional[Callable[[HST], float]] = None,
              rescale: bool = True) -> HST:
    """Sample ``n_samples`` FRT trees and return the one minimising ``criterion``
    (default: least-squares fit to ``D``), optionally rescaled to fit ``D``."""
    trees = sample_trees(D, n_samples, seed)
    if len(trees) == 1:
        return trees[0]
    crit = criterion or lsq_cost(D)
    costs = [crit(T) for T in trees]
    T = trees[int(np.argmin(costs))]
    return rescale_to_fit(T, D) if rescale else T
