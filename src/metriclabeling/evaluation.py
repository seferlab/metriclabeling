"""Evaluation (Section 2.7): each (protein, GO term) pair of a test protein is
one example; ROC points at score thresholds 0, 0.05, ..., 1; averages over the
five folds."""

from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np

THRESHOLDS = np.round(np.arange(0.0, 1.0 + 1e-9, 0.05), 2)


def roc_points(y_true: np.ndarray, scores: np.ndarray, thresholds: Sequence[float] = THRESHOLDS):
    """(FPR, TPR) at each threshold (prediction = score >= threshold; for the
    threshold 0 every pair is predicted positive)."""
    y = np.asarray(y_true, dtype=bool).ravel()
    s = np.asarray(scores, dtype=float).ravel()
    P, N = y.sum(), (~y).sum()
    fpr, tpr = [], []
    for t in thresholds:
        pred = s >= t - 1e-12 if t > 0 else np.ones_like(y)
        tpr.append((pred & y).sum() / max(P, 1))
        fpr.append((pred & ~y).sum() / max(N, 1))
    return np.array(fpr), np.array(tpr)


def auc_from_points(fpr: np.ndarray, tpr: np.ndarray) -> float:
    """Trapezoidal area under a ROC polyline (adds the (0,0) and (1,1) corners)."""
    x = np.r_[0.0, fpr, 1.0]
    y = np.r_[0.0, tpr, 1.0]
    o = np.lexsort((y, x))
    return float(np.trapezoid(y[o], x[o]) if hasattr(np, "trapezoid") else np.trapz(y[o], x[o]))


def exact_auc(y_true: np.ndarray, scores: np.ndarray) -> float:
    """Exact ROC AUC (Mann-Whitney) over all pairs; ties count 1/2."""
    from scipy.stats import rankdata

    y = np.asarray(y_true, dtype=bool).ravel()
    s = np.asarray(scores, dtype=float).ravel()
    P, N = y.sum(), (~y).sum()
    if P == 0 or N == 0:
        return float("nan")
    r = rankdata(s)
    return float((r[y].sum() - P * (P + 1) / 2) / (P * N))


def evaluate_fold(Y_true: np.ndarray, S: np.ndarray, test: np.ndarray) -> Dict:
    y = Y_true[test]
    s = S[test]
    fpr, tpr = roc_points(y, s)
    return {
        "fpr": fpr.tolist(),
        "tpr": tpr.tolist(),
        "auc_grid": auc_from_points(fpr, tpr),
        "auc": exact_auc(y, s),
        "n_pos": int(y.sum()),
        "n_neg": int((~y).sum()),
    }


def aggregate(fold_results: List[Dict]) -> Dict:
    """Average ROC points and AUCs over folds."""
    fpr = np.mean([r["fpr"] for r in fold_results], axis=0)
    tpr = np.mean([r["tpr"] for r in fold_results], axis=0)
    return {
        "fpr": fpr.tolist(),
        "tpr": tpr.tolist(),
        "auc_grid": auc_from_points(fpr, tpr),
        "auc": float(np.mean([r["auc"] for r in fold_results])),
        "auc_std": float(np.std([r["auc"] for r in fold_results])),
    }
