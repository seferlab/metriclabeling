"""Metric Labeling for protein function prediction (the paper's method)."""

from __future__ import annotations

import logging
import time
from typing import Dict, Optional, Tuple

import numpy as np

from ..folds import Fold
from ..graphcut import reduced_energy
from ..hst import best_tree, frt_tree, rescale_to_fit, sample_trees, tree_from_metric_if_uniform
from ..lp import emd_tree_cost, round_on_tree, solve_relaxation, transport_lp_size
from .common import assignment_costs, choice_reduced, local_scores, row_normalize, split_reduced

log = logging.getLogger(__name__)


class MetricLabeling:
    """Metric Labeling with (optional) LSD metric approximation.

    Parameters
    ----------
    distance, alpha:
        GO distance (``sp``/``lca``/``lin``/``kb``/``uniform``) and the weight of
        ``d_KB`` in Eq. (16).
    costs: ``nonuniform`` or ``uniform`` assignment costs (Section 2.4).
    metric: ``lsd`` (LSD approximation of semimetrics, the paper's method),
        ``raw`` (use the semimetric as is; no guarantees).
    lp: ``auto`` | ``transport`` | ``tree``.  ``transport`` is the exact LP of
        Chekuri et al. (Eqs. 6-11); ``tree`` solves the LP on an FRT tree
        (Kleinberg-Tardos).  ``auto`` uses ``transport`` when it has at most
        ``transport_max_vars`` variables.
    backend: LP solver: ``highs`` | ``gurobi`` | ``cplex`` | ``pdhg`` | ``auto``
        (``auto``: HiGHS, or PDHG for large tree LPs).
    score: how to turn the solution into scores in [0, 1]:
        ``lp`` (LP marginals), ``local`` (local min-marginals, see
        :func:`local_scores`), ``lp+local`` (max of both; default) or
        ``frequency`` (fraction of ``rounds`` randomized HST roundings).
    known: ``split`` (multi-labeled training proteins act as one fixed copy per
        label) or ``choice`` (they may pick one of their labels in the LP).
    """

    name = "metric_labeling"

    def __init__(
        self,
        distance: str = "lca",
        alpha: float = 0.0,
        costs: str = "nonuniform",
        metric: str = "lsd",
        lsd_method: str = "triangle_fixing",
        lp: str = "auto",
        backend: str = "auto",
        transport_max_vars: int = 1_500_000,
        tree_samples: int = 32,
        score: str = "lp+local",
        temperature: Optional[float] = None,
        rounds: int = 50,
        known: str = "split",
        lin_mode: str = "reciprocal",
        pdhg: Optional[Dict] = None,
        seed: int = 0,
        **unused,
    ) -> None:
        self.distance = distance
        self.alpha = float(alpha)
        self.costs = costs
        self.metric_mode = metric
        self.lsd_method = lsd_method
        self.lp = lp
        self.backend = backend
        self.transport_max_vars = int(transport_max_vars)
        self.tree_samples = int(tree_samples)
        self.score = score
        self.temperature = temperature
        self.rounds = int(rounds)
        self.known = known
        self.lin_mode = lin_mode
        self.pdhg = pdhg or {}
        self.seed = seed
        self.info: Dict = {}

    # ------------------------------------------------------------------ metric
    def label_metric(self, fold: Fold) -> Tuple[np.ndarray, Dict]:
        if self.distance.lower() in ("uniform", "potts"):
            k = fold.k
            D = np.ones((k, k))
            np.fill_diagonal(D, 0.0)
            return D, {"lsd": "none"}
        if self.metric_mode == "lsd":
            D, res = fold.metric(self.distance, self.alpha, self.lin_mode, self.lsd_method)
            return D, {"lsd": res.as_dict()}
        return fold.semimetric(self.distance, self.alpha, self.lin_mode), {"lsd": "skipped"}

    # ---------------------------------------------------------------- predict
    def predict(self, fold: Fold) -> np.ndarray:
        t0 = time.perf_counter()
        rng = np.random.default_rng(self.seed + 1000 * fold.fold_id)
        D, minfo = self.label_metric(fold)
        U = fold.unknown
        unary_u = assignment_costs(fold, U, self.costs)
        if self.known == "split":
            red, free_nodes = split_reduced(fold, unary_u, D)
        else:
            red, free_nodes, _ = choice_reduced(fold, unary_u, D)

        # ---- LP relaxation -------------------------------------------------
        lp = self.lp
        if lp == "auto":
            lp = "transport" if transport_lp_size(red) <= self.transport_max_vars else "tree"
        tree = None
        uni = tree_from_metric_if_uniform(D)
        if lp == "tree":
            tree = uni if uni is not None else best_tree(D, self.tree_samples, rng)
        backend = self.backend
        if backend == "auto":
            if lp == "tree":
                size = len(red.fu) * 2 * (tree.n_nodes - 1)
                backend = "highs" if size <= 300_000 else "pdhg"
            else:
                backend = "highs"
        rel = solve_relaxation(red, D, method=lp, tree=tree, backend=backend,
                               transport_max_vars=self.transport_max_vars, pdhg_kwargs=self.pdhg)
        X = rel.X

        # ---- integral solution: Chekuri et al. tree selection + rounding ----
        if tree is None:
            cands = [uni] if uni is not None else [rescale_to_fit(T, D) for T in sample_trees(D, self.tree_samples, rng)]
            costs = [emd_tree_cost(red, X, T) for T in cands]
            round_tree = cands[int(np.argmin(costs))]
        else:
            round_tree = tree
        best_lab, best_e = None, np.inf
        for _ in range(10):
            lab = round_on_tree(X, round_tree, rng, red.allowed)
            e = reduced_energy(red, lab, D)
            if e < best_e:
                best_lab, best_e = lab, e
        self.labels_ = np.full(fold.n, -1, dtype=np.int64)
        self.labels_[free_nodes] = best_lab

        # ---- scores --------------------------------------------------------
        marg = row_normalize(fold.Y_train) if self.known == "split" else np.zeros((fold.n, fold.k))
        if self.known != "split":
            marg[red.fixed] = 0.0
            marg[red.fixed, red.fixed_label[red.fixed]] = 1.0
        marg[free_nodes] = X
        S = np.zeros((fold.n, fold.k))
        if self.score == "lp":
            S[free_nodes] = X
        elif self.score in ("local", "lp+local"):
            loc = local_scores(fold, U, unary_u, marg, D, self.temperature)
            S[U] = loc if self.score == "local" else np.maximum(loc, marg[U])
        elif self.score == "frequency":
            cnt = np.zeros((len(free_nodes), fold.k))
            for r in range(self.rounds):
                if tree is not None or uni is not None:
                    T = tree if tree is not None else uni
                else:
                    T = frt_tree(D, rng)
                lab = round_on_tree(X, T, rng, red.allowed)
                cnt[np.arange(len(free_nodes)), lab] += 1
            S[free_nodes] = cnt / self.rounds
        else:
            raise ValueError(f"unknown score mode {self.score!r}")
        S[fold.train] = fold.Y_train[fold.train]
        self.info = {
            **minfo,
            "lp": rel.method,
            "backend": backend,
            "lp_objective": rel.objective,
            "rounded_energy": best_e,
            "n_free": int(red.n_free),
            "n_free_edges": int(len(red.fu)),
            "seconds": time.perf_counter() - t0,
            **{f"lp_{k}": v for k, v in rel.info.items() if isinstance(v, (int, float, str))},
        }
        return np.clip(S, 0.0, 1.0)
