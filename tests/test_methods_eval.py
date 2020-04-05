import numpy as np
import pytest

from metriclabeling.evaluation import aggregate, auc_from_points, evaluate_fold, exact_auc, roc_points
from metriclabeling.folds import make_folds, misannotate, remove_edges
from metriclabeling.methods import make_method


def test_auc_matches_sklearn(rng):
    sk = pytest.importorskip("sklearn.metrics")
    y = rng.random(500) < 0.2
    s = np.round(rng.random(500) + 0.3 * y, 2)
    assert abs(exact_auc(y, s) - sk.roc_auc_score(y, s)) < 1e-12
    fpr, tpr = roc_points(y, s)
    assert fpr[0] == 1 and tpr[0] == 1
    assert 0.5 < auc_from_points(fpr, tpr) <= 1


SPECS = [
    {"method": "majority"},
    {"method": "neighborhood"},
    {"method": "functional_flow"},
    {"method": "mrf"},
    {"method": "genmulticut", "runs": 3},
    {"method": "semimap", "distance": "lca", "alpha": 0.3, "runs": 3},
    {"method": "metric_labeling", "distance": "lca", "alpha": 0.3},
    {"method": "metric_labeling", "distance": "lin", "alpha": 0.1, "lp": "tree", "backend": "pdhg"},
    {"method": "metric_labeling", "distance": "sp", "costs": "uniform", "lp": "transport"},
    {"method": "metric_labeling", "distance": "lca", "known": "choice", "score": "frequency", "rounds": 5},
]


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: "-".join(str(v) for v in s.values()))
def test_methods_run(synth, spec):
    fold = make_folds(synth, 5, seed=0)[0]
    S = make_method(dict(spec)).predict(fold)
    assert S.shape == (synth.n, synth.k)
    assert np.all((S >= 0) & (S <= 1))
    r = evaluate_fold(synth.Y, S, fold.test)
    assert r["auc"] > 0.55


def test_metric_beats_potts_on_synthetic(synth):
    folds = make_folds(synth, 5, seed=1)
    res = {}
    for name, spec in {"go": {"method": "metric_labeling", "distance": "sp"},
                       "potts": {"method": "metric_labeling", "distance": "uniform"}}.items():
        res[name] = aggregate([evaluate_fold(synth.Y, make_method(dict(spec)).predict(f), f.test) for f in folds])["auc"]
    assert res["go"] > res["potts"]


def test_perturbations(synth):
    A2 = remove_edges(synth.adj, 0.3, seed=0)
    assert abs(A2.nnz / synth.adj.nnz - 0.7) < 0.1
    assert (A2 != A2.T).nnz == 0
    f = make_folds(synth, 5)[0]
    g = misannotate(f, 0.5)
    changed = (g.Y_train[f.train] != f.Y_train[f.train]).any(axis=1).mean()
    assert 0.3 < changed <= 0.5 + 1e-9
    assert np.array_equal(g.Y_train[f.test], f.Y_train[f.test])
