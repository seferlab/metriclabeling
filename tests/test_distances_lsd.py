import numpy as np
import pytest

from metriclabeling.distances import (combine, d_kb, d_lca, d_lin, d_sp, is_metric, term_probabilities,
                                      triangle_violations)
from metriclabeling.lsd import lsd, lsd_objective, max_triangle_violation


def _check_dissimilarity(D):
    assert np.allclose(D, D.T)
    assert np.allclose(np.diag(D), 0)
    assert (D >= 0).all()


def test_go_distances(synth):
    go = synth.go
    Dsp = d_sp(go)
    _check_dissimilarity(Dsp)
    assert is_metric(Dsp)  # shortest-path distance is a metric
    assert Dsp.max() <= 1.0 + 1e-12
    Dl = d_lca(go)
    _check_dissimilarity(Dl)
    assert Dl.max() <= 1.0 + 1e-12
    Dlin = d_lin(go, term_probabilities(go, synth.labeled))
    _check_dissimilarity(Dlin)
    Dkb = d_kb(synth.adj, synth.Y)
    _check_dissimilarity(Dkb)
    C = combine(Dl, Dkb, 0.3)
    assert np.allclose(C, 0.7 * Dl + 0.3 * Dkb)


def _random_semimetric(rng, k):
    X = rng.random((k, k))
    S = (X + X.T) / 2
    np.fill_diagonal(S, 0)
    return S


@pytest.mark.parametrize("k", [3, 7, 12])
def test_lsd_triangle_fixing_is_projection(rng, k):
    S = _random_semimetric(rng, k)
    res = lsd(S, tol=1e-10)
    M = res.M
    assert max_triangle_violation(M) < 1e-6
    _check_dissimilarity(M)
    # optimality: no other metric we can construct is closer (shortest-path closure, scaled copies)
    from scipy.sparse.csgraph import shortest_path

    closure = shortest_path(S, directed=False)
    assert res.objective <= lsd_objective(S, closure) + 1e-9
    for t in (0.9, 1.0, 1.1):
        assert res.objective <= lsd_objective(S, t * M) + 1e-9
    # KKT check against an exact QP solver when available
    pytest.importorskip("cvxpy")
    ref = lsd(S, method="cvxpy")
    assert abs(ref.objective - res.objective) < 1e-5
    assert np.abs(ref.M - M).max() < 1e-3


def test_lsd_cutting_plane_matches(rng):
    pytest.importorskip("cvxpy")
    S = _random_semimetric(rng, 10)
    a = lsd(S)
    b = lsd(S, method="cutting_plane", solver="cvxpy")
    assert abs(a.objective - b.objective) < 1e-5


def test_lsd_keeps_metrics(rng):
    X = rng.random((8, 2))
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    res = lsd(D)
    assert res.method == "none" and np.allclose(res.M, D)
    assert triangle_violations(D)[0] == 0
