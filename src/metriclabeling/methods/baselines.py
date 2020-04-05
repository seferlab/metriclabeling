"""Comparison methods of Section 2.6.

* :class:`Majority`       Schwikowski et al. (2000)
* :class:`Neighborhood`   Hishigaki et al. (2001), radius 2, chi-square
* :class:`FunctionalFlow` Nabieva et al. (2005)
* :class:`MRF`            Lee et al. (2006), diffusion-kernel logistic regression
                          with correlated functions
* :class:`GenMultiCut`    Vazquez et al. (2003) / Karaoz et al. (2004), solved
                          50 times with perturbed weights as in Nabieva et al.
* :class:`SemiMAP`        Kumar & Koller (2009): r-HST mixture + hierarchical
                          graph cuts for semimetric labeling
"""

from __future__ import annotations

import logging
import time
from typing import Dict, Optional

import numpy as np
import scipy.sparse as sp

from ..folds import Fold
from ..graphcut import alpha_expansion, hierarchical_moves, reduced_energy
from ..hst import frt_tree, rescale_to_fit, star_tree
from ..lp import round_on_tree, solve_relaxation
from .common import assignment_costs, local_scores, perturbed_adj, row_normalize, split_reduced

log = logging.getLogger(__name__)


def _scores_from_counts(fold: Fold, U: np.ndarray, freq: np.ndarray, unary_u: np.ndarray, D: np.ndarray,
                        score: str, temperature: Optional[float]) -> np.ndarray:
    """Scores of the multi-run methods: label ``frequency`` over runs (as in the
    paper), ``local`` min-marginals around the frequency solution (the same
    scoring as :class:`MetricLabeling`), or ``frequency+local`` (max)."""
    S = np.zeros((fold.n, fold.k))
    if score == "frequency":
        S[U] = freq
    elif score in ("local", "frequency+local", "lp+local"):
        marg = row_normalize(fold.Y_train)
        marg[U] = freq
        loc = local_scores(fold, U, unary_u, marg, D, temperature)
        S[U] = loc if score == "local" else np.maximum(loc, freq)
    else:
        raise ValueError(f"unknown score mode {score!r}")
    S[fold.train] = fold.Y_train[fold.train]
    return S


# --------------------------------------------------------------------------- #
class Majority:
    """Score of function ``i`` = weighted fraction of annotated neighbours having ``i``."""

    name = "majority"

    def __init__(self, **unused) -> None:
        self.info: Dict = {}

    def predict(self, fold: Fold) -> np.ndarray:
        A = fold.adj.tocsr()
        known = fold.known.astype(float)
        num = np.asarray(A @ fold.Y_train.astype(float))
        den = np.asarray(A @ known).ravel()
        S = np.divide(num, den[:, None], out=np.zeros_like(num), where=den[:, None] > 0)
        S[fold.train] = fold.Y_train[fold.train]
        return S


# --------------------------------------------------------------------------- #
class Neighborhood:
    """Chi-square over-representation of each function within radius ``r``.

    ``chi2 = (n - e)^2 / e`` for ``n > e`` (0 otherwise), where ``n`` counts
    annotated proteins with the function within distance ``r`` and ``e`` is the
    expected count under the global frequency.  Mapped to [0, 1) by
    ``chi2 / (chi2 + 3.84)`` (0.5 at the 5% critical value).
    """

    name = "neighborhood"

    def __init__(self, radius: int = 2, **unused) -> None:
        self.radius = int(radius)
        self.info: Dict = {}

    def predict(self, fold: Fold) -> np.ndarray:
        A = (fold.adj > 0).astype(np.float64).tocsr()
        R = A.copy()
        P = A.copy()
        for _ in range(self.radius - 1):
            P = (P @ A)
            P.data[:] = 1.0
            R = R + P
        R = (R > 0).astype(np.float64).tolil()
        R.setdiag(0)
        R = R.tocsr()
        R.eliminate_zeros()
        known = fold.known
        Yk = fold.Y_train.astype(float)
        n_obs = np.asarray(R @ Yk)
        N = np.asarray(R @ known.astype(float)).ravel()
        pi = Yk[known].mean(axis=0)
        e = N[:, None] * pi[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            chi = np.where((n_obs > e) & (e > 0), (n_obs - e) ** 2 / e, 0.0)
        S = chi / (chi + 3.84)
        S[fold.train] = fold.Y_train[fold.train]
        return S


# --------------------------------------------------------------------------- #
class FunctionalFlow:
    """FunctionalFlow (Nabieva et al. 2005) with ``d`` propagation steps."""

    name = "functional_flow"

    def __init__(self, steps: int = 6, **unused) -> None:
        self.steps = int(steps)
        self.info: Dict = {}

    def predict(self, fold: Fold) -> np.ndarray:
        A = sp.triu(fold.adj, k=1).tocoo()
        n, k = fold.n, fold.k
        src = np.r_[A.row, A.col]
        dst = np.r_[A.col, A.row]
        w = np.r_[A.data, A.data].astype(float)
        W = np.asarray(fold.adj.sum(axis=1)).ravel()
        frac = (w / np.maximum(W[src], 1e-12))[:, None]
        Din = sp.csr_matrix((np.ones(len(dst)), (dst, np.arange(len(dst)))), shape=(n, len(dst)))
        Dout = sp.csr_matrix((np.ones(len(src)), (src, np.arange(len(src)))), shape=(n, len(src)))
        R = np.where(fold.Y_train, np.inf, 0.0)
        score = np.zeros((n, k))
        for _ in range(self.steps):
            Rs, Rd = R[src], R[dst]
            with np.errstate(invalid="ignore"):
                g = np.where(Rs >= Rd, np.minimum(w[:, None], frac * Rs), 0.0)
            g = np.nan_to_num(g, nan=0.0, posinf=0.0)
            inflow = np.asarray(Din @ g)
            outflow = np.asarray(Dout @ g)
            finite = np.isfinite(R)
            R = np.where(finite, R + inflow - outflow, R)
            score += inflow
        score[fold.train] = 0.0
        mx = score.max()
        S = score / mx if mx > 0 else score
        S[fold.train] = fold.Y_train[fold.train]
        return S


# --------------------------------------------------------------------------- #
class MRF:
    """Kernel logistic regression MRF of Lee et al. (2006).

    For every function ``i`` a logistic regression is trained on the
    training proteins with features

    * ``K1 = sum_{v != u} K(u,v) Y_v^i``   (kernel-weighted neighbours with ``i``)
    * ``K0 = sum_{v != u} K(u,v) (1 - Y_v^i)``  (annotated neighbours without ``i``)
    * ``K^j`` for the ``n_correlated`` functions ``j`` most correlated with
      ``i`` in the training labels,

    where ``K = exp(-beta L)`` is the diffusion kernel of the (normalised)
    graph Laplacian.
    """

    name = "mrf"

    def __init__(self, beta: float = 1.0, n_correlated: int = 5, normalized: bool = True,
                 C: float = 1.0, **unused) -> None:
        self.beta = float(beta)
        self.n_correlated = int(n_correlated)
        self.normalized = bool(normalized)
        self.C = float(C)
        self.info: Dict = {}

    def _kernel(self, fold: Fold):
        key = ("diffusion", self.beta, self.normalized)
        if key not in fold.shared:
            A = fold.adj.toarray()
            d = A.sum(axis=1)
            if self.normalized:
                dis = 1.0 / np.sqrt(np.maximum(d, 1e-12))
                L = np.eye(len(d)) - dis[:, None] * A * dis[None, :]
            else:
                L = np.diag(d) - A
            evals, evecs = np.linalg.eigh(L)
            ex = np.exp(-self.beta * evals)
            diag = (evecs ** 2) @ ex
            fold.shared[key] = (evecs, ex, diag)
        return fold.shared[key]

    def predict(self, fold: Fold) -> np.ndarray:
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler

        t0 = time.perf_counter()
        V, ex, diag = self._kernel(fold)
        known = fold.known.astype(float)
        Y = fold.Y_train.astype(float)
        Z = np.c_[Y, known]
        KZ = V @ (ex[:, None] * (V.T @ Z)) - diag[:, None] * Z  # exclude self
        KY, Ka = KZ[:, :-1], KZ[:, -1]
        tr = fold.train
        Yt = Y[tr]
        C = np.corrcoef(Yt.T) if Yt.shape[0] > 1 else np.zeros((fold.k, fold.k))
        C = np.nan_to_num(C)
        np.fill_diagonal(C, -np.inf)
        S = np.zeros((fold.n, fold.k))
        for i in range(fold.k):
            y = Yt[:, i]
            if y.min() == y.max():
                S[:, i] = y.mean()
                continue
            corr = np.argsort(-C[i])[: self.n_correlated]
            F = np.c_[KY[:, i], Ka - KY[:, i], KY[:, corr]]
            sc = StandardScaler().fit(F[tr])
            clf = LogisticRegression(C=self.C, max_iter=1000)
            clf.fit(sc.transform(F[tr]), y)
            S[:, i] = clf.predict_proba(sc.transform(F))[:, 1]
        S[fold.train] = fold.Y_train[fold.train]
        self.info = {"seconds": time.perf_counter() - t0}
        return S


# --------------------------------------------------------------------------- #
class GenMultiCut:
    """Generalized multiway cut: Metric Labeling with the uniform (Potts) metric.

    Solved ``runs`` times with edge weights perturbed by
    ``U(-w_max 1e-5, w_max 1e-5)``; the score of a label is the fraction of
    runs assigning it (Nabieva et al. 2005).  ``solver``: ``expansion``
    (alpha-expansion graph cuts) or ``lp`` (exact LP relaxation of the Potts
    model + Kleinberg-Tardos rounding).
    """

    name = "genmulticut"

    def __init__(self, runs: int = 50, solver: str = "expansion", perturb: float = 1e-5, seed: int = 0,
                 score: str = "frequency", temperature: Optional[float] = None, **unused) -> None:
        self.runs = int(runs)
        self.score = score
        self.temperature = temperature
        self.solver = solver
        self.perturb = float(perturb)
        self.seed = seed
        self.info: Dict = {}

    def predict(self, fold: Fold) -> np.ndarray:
        t0 = time.perf_counter()
        rng = np.random.default_rng(self.seed + 1000 * fold.fold_id)
        k = fold.k
        D = np.ones((k, k))
        np.fill_diagonal(D, 0.0)
        U = fold.unknown
        zero = np.zeros((len(U), k))
        cnt = np.zeros((len(U), k))
        tree = star_tree(k)
        for r in range(self.runs):
            adj = perturbed_adj(fold.adj, rng, self.perturb)
            red, _ = split_reduced(fold, zero, D, adj=adj)
            if self.solver == "lp":
                rel = solve_relaxation(red, D, method="tree", tree=tree, backend="highs")
                lab = round_on_tree(rel.X, tree, rng)
            else:
                lab = alpha_expansion(red, D, rng=rng)
            cnt[np.arange(len(U)), lab] += 1
        S = _scores_from_counts(fold, U, cnt / self.runs, np.zeros((len(U), k)), D, self.score, self.temperature)
        self.info = {"seconds": time.perf_counter() - t0}
        return S


# --------------------------------------------------------------------------- #
class SemiMAP:
    """Semimetric MAP estimation of Kumar & Koller (UAI 2009).

    The semimetric is approximated by a mixture of r-HST (FRT) metrics, each
    tree labeling problem is solved with hierarchical graph cuts, and the
    labeling with the lowest energy under the *original* semimetric is kept.
    As in the paper this is repeated ``runs`` times with perturbed edge
    weights (and freshly sampled trees); scores are label frequencies.
    """

    name = "semimap"

    def __init__(self, distance: str = "lca", alpha: float = 0.0, costs: str = "nonuniform",
                 runs: int = 50, trees: int = 2, perturb: float = 1e-5, lin_mode: str = "reciprocal",
                 seed: int = 0, score: str = "frequency", temperature: Optional[float] = None, **unused) -> None:
        self.score = score
        self.temperature = temperature
        self.distance = distance
        self.alpha = float(alpha)
        self.costs = costs
        self.runs = int(runs)
        self.trees = int(trees)
        self.perturb = float(perturb)
        self.lin_mode = lin_mode
        self.seed = seed
        self.info: Dict = {}

    def predict(self, fold: Fold) -> np.ndarray:
        t0 = time.perf_counter()
        rng = np.random.default_rng(self.seed + 1000 * fold.fold_id)
        S_semi = fold.semimetric(self.distance, self.alpha, self.lin_mode)
        U = fold.unknown
        unary_u = assignment_costs(fold, U, self.costs)
        cnt = np.zeros((len(U), fold.k))
        energies = []
        for r in range(self.runs):
            adj = perturbed_adj(fold.adj, rng, self.perturb)
            red_s, _ = split_reduced(fold, unary_u, S_semi, adj=adj)
            best, best_e = None, np.inf
            for _ in range(self.trees):
                T = rescale_to_fit(frt_tree(S_semi, rng), S_semi)
                red_t, _ = split_reduced(fold, unary_u, T.label_distance(), adj=adj)
                lab = hierarchical_moves(red_t, T)
                e = reduced_energy(red_s, lab, S_semi)
                if e < best_e:
                    best, best_e = lab, e
            energies.append(best_e)
            cnt[np.arange(len(U)), best] += 1
        S = _scores_from_counts(fold, U, cnt / self.runs, unary_u, S_semi, self.score, self.temperature)
        self.info = {"seconds": time.perf_counter() - t0, "mean_energy": float(np.mean(energies))}
        return S
