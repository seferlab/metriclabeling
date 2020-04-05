"""Cross-validation folds and per-fold caches (distances, LSD metrics, kernels)."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp

from .dataset import Dataset, graph_diameter
from .distances import distance_from_spec, network_distances
from .lsd import LSDResult, lsd

log = logging.getLogger(__name__)


@dataclass
class Fold:
    """One cross-validation split.

    ``train`` are labeled proteins whose labels are visible (possibly
    corrupted, see :func:`misannotate`); ``test`` are labeled proteins whose
    labels are hidden and evaluated.  Proteins outside both sets have no
    selected label and are neither known nor evaluated.
    """

    ds: Dataset
    adj: sp.csr_matrix
    train: np.ndarray
    test: np.ndarray
    Y_train: np.ndarray
    fold_id: int = 0
    seed: int = 0
    shared: Dict = field(default_factory=dict)  # cache shared by folds on the same network
    cache: Dict = field(default_factory=dict)  # fold-local cache

    # ------------------------------------------------------------------ basics
    @property
    def n(self) -> int:
        return self.ds.n

    @property
    def k(self) -> int:
        return self.ds.k

    @property
    def known(self) -> np.ndarray:
        m = np.zeros(self.n, dtype=bool)
        m[self.train] = True
        return m

    @property
    def unknown(self) -> np.ndarray:
        return np.flatnonzero(~self.known)

    def edges(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        A = sp.triu(self.adj, k=1).tocoo()
        return A.row.astype(np.int64), A.col.astype(np.int64), A.data.astype(np.float64)

    # ------------------------------------------------------------- shared data
    def net_diameter(self) -> float:
        if "diameter" not in self.shared:
            self.shared["diameter"] = graph_diameter(self.adj)
        return self.shared["diameter"]

    def sp_cache(self):
        """BFS distances from all labeled proteins (superset of every fold's training set)."""
        if "sp" not in self.shared:
            nodes = self.ds.labeled
            self.shared["sp"] = (nodes, network_distances(self.adj, nodes))
        return self.shared["sp"]

    # ---------------------------------------------------------------- distances
    def semimetric(self, name: str, alpha: float = 0.0, lin_mode: str = "reciprocal") -> np.ndarray:
        key = ("semi", name, float(alpha), lin_mode)
        if key not in self.cache:
            need_kb = alpha > 0 or name.lower() == "kb"
            self.cache[key] = distance_from_spec(
                name,
                self.ds.go,
                adj=self.adj if need_kb else None,
                Y_train=self.Y_train if need_kb else None,
                train_idx=self.train,
                alpha=alpha,
                net_diameter=self.net_diameter() if need_kb else None,
                sp_cache=self.sp_cache() if need_kb else None,
                lin_mode=lin_mode,
            )
        return self.cache[key]

    def metric(self, name: str, alpha: float = 0.0, lin_mode: str = "reciprocal",
               lsd_method: str = "triangle_fixing") -> Tuple[np.ndarray, LSDResult]:
        """The LSD metric approximation of the requested (semi)metric (cached)."""
        key = ("lsd", name, float(alpha), lin_mode, lsd_method)
        if key not in self.cache:
            S = self.semimetric(name, alpha, lin_mode)
            res = lsd(S, method=lsd_method)
            log.debug("LSD %s alpha=%.2f: %s", name, alpha, res.as_dict())
            self.cache[key] = (res.M, res)
        return self.cache[key]


def make_folds(
    ds: Dataset,
    n_folds: int = 5,
    seed: int = 0,
    adj: Optional[sp.csr_matrix] = None,
    shared: Optional[Dict] = None,
) -> List[Fold]:
    """Stratification-free k-fold split of the labeled proteins."""
    rng = np.random.default_rng(seed)
    lab = ds.labeled.copy()
    rng.shuffle(lab)
    parts = np.array_split(lab, n_folds)
    adj = ds.adj if adj is None else adj
    shared = {} if shared is None else shared
    folds = []
    for f in range(n_folds):
        test = np.sort(parts[f])
        train = np.sort(np.concatenate([parts[g] for g in range(n_folds) if g != f]))
        Y_train = np.zeros_like(ds.Y)
        Y_train[train] = ds.Y[train]
        fold_adj = adj
        if "sources" in ds.meta:  # integrated network: per-fold noisy-or reliabilities
            from .integrated import noisy_or_adjacency

            known = np.zeros(ds.n, dtype=bool)
            known[train] = True
            sources = ds.meta["sources"]
            if adj is not ds.adj:  # edge-removal experiments: restrict sources to the kept edges
                mask = (adj > 0).astype(float)
                sources = {k: M.multiply(mask).tocsr() for k, M in sources.items()}
            fold_adj, rel = noisy_or_adjacency(sources, Y_train, known)
            log.info("fold %d source reliabilities: %s", f, {k: round(v, 3) for k, v in rel.items()})
        folds.append(Fold(ds, fold_adj, train, test, Y_train, fold_id=f, seed=seed, shared=shared))
    return folds


# --------------------------------------------------------------------------- #
# perturbations for the robustness experiments (Section 3.3)
# --------------------------------------------------------------------------- #
def remove_edges(adj: sp.csr_matrix, frac: float, seed: int = 0) -> sp.csr_matrix:
    """Remove a random fraction of the (undirected) edges."""
    if frac <= 0:
        return adj
    rng = np.random.default_rng(seed)
    A = sp.triu(adj, k=1).tocoo()
    keep = rng.random(A.nnz) >= frac
    r, c, w = A.row[keep], A.col[keep], A.data[keep]
    n = adj.shape[0]
    return sp.csr_matrix((np.r_[w, w], (np.r_[r, c], np.r_[c, r])), shape=(n, n))


def misannotate(fold: Fold, frac: float, seed: int = 0) -> Fold:
    """Replace the label sets of a random fraction of training proteins by random labels
    (same number of labels); test labels stay correct."""
    if frac <= 0:
        return fold
    rng = np.random.default_rng(seed + 7919 * (fold.fold_id + 1))
    Y = fold.Y_train.copy()
    chosen = rng.choice(fold.train, size=int(round(frac * len(fold.train))), replace=False)
    k = Y.shape[1]
    for p in chosen:
        m = int(Y[p].sum())
        Y[p] = False
        Y[p, rng.choice(k, size=max(m, 1), replace=False)] = True
    return Fold(fold.ds, fold.adj, fold.train, fold.test, Y, fold.fold_id, fold.seed, fold.shared, {})
