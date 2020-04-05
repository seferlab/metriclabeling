import itertools

import numpy as np
import pytest

from metriclabeling.graphcut import alpha_expansion, binary_cut, hierarchical_moves, reduced_energy
from metriclabeling.hst import best_tree, frt_tree, star_tree
from metriclabeling.lp import (LabelingProblem, build_transport_lp, reduce_problem, round_on_tree, solve_ilp,
                               solve_relaxation)


def _euclid(rng, k, dim=2):
    X = rng.random((k, dim))
    return np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))


def test_frt_dominates_and_compresses(rng):
    for _ in range(20):
        k = int(rng.integers(2, 30))
        D = _euclid(rng, k)
        T = frt_tree(D, rng)
        Dt = T.label_distance()
        assert np.all(Dt >= D - 1e-9)
        assert T.n_nodes <= 2 * k
        assert all(len(c) != 1 for c in T.children)  # no unary internal nodes
    T = star_tree(4)
    assert np.allclose(T.label_distance(), 1 - np.eye(4))


def test_best_tree_rescaled(rng):
    D = _euclid(rng, 15)
    T = best_tree(D, 8, rng)
    iu = np.triu_indices(15, 1)
    ratio = T.label_distance()[iu] / D[iu]
    assert 0.2 < np.median(ratio) < 5


def _problem(rng, n, k, p=0.4):
    D = _euclid(rng, k)
    U = rng.random((n, k))
    for u in range(n):
        r = rng.random()
        if r < 0.3:
            U[u] = np.inf
            U[u, rng.integers(k)] = 0
        elif r < 0.45:
            m = rng.random(k) < 0.5
            m[rng.integers(k)] = True
            U[u, ~m] = np.inf
    eu, ev = np.nonzero(np.triu(rng.random((n, n)) < p, 1))
    return LabelingProblem(U, eu, ev, rng.random(len(eu)) + 0.2, D)


def test_ilp_lp_bruteforce(rng):
    from scipy.optimize import linprog

    for _ in range(6):
        P = _problem(rng, 6, 3)
        feas = [l for l in itertools.product(range(3), repeat=6)
                if np.all(np.isfinite(P.unary[np.arange(6), list(l)]))]
        best = min(P.energy(np.array(l)) for l in feas)
        lab = solve_ilp(P)
        assert abs(P.energy(lab) - best) < 1e-8
        red = reduce_problem(P)
        assert abs(reduced_energy(red, lab[red.free], P.D) - best) < 1e-8
        c, A, b, _ = build_transport_lp(red, P.D, np.arange(red.n_free), np.arange(len(red.fu)))
        lp = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs").fun + red.const
        assert lp <= best + 1e-8
        ae = alpha_expansion(red, P.D, rng=rng)
        off = P.D[~np.eye(3, dtype=bool)]
        bound = 2 * off.max() / off[off > 0].min()  # Boykov et al. (2001) bound for metrics
        assert reduced_energy(red, ae, P.D) <= bound * best + 1e-8


def test_tree_lp_equals_transport_on_tree_metric(rng):
    from scipy.optimize import linprog

    for _ in range(4):
        P = _problem(rng, 12, 6, 0.35)
        T = frt_tree(P.D, rng)
        Dt = T.label_distance()
        red = reduce_problem(P, Dt)
        c, A, b, _ = build_transport_lp(red, Dt, np.arange(red.n_free), np.arange(len(red.fu)))
        v = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs").fun + red.const
        r_tree = solve_relaxation(red, Dt, method="tree", tree=T)
        r_pdhg = solve_relaxation(red, Dt, method="tree", tree=T, backend="pdhg",
                                  pdhg_kwargs={"tol": 1e-7, "max_iter": 100000})
        assert abs(r_tree.objective - v) < 1e-6
        assert abs(r_pdhg.objective - v) < 1e-3 * max(1, abs(v))
        lab = round_on_tree(r_tree.X, T, rng, red.allowed)
        assert np.all(red.allowed[np.arange(red.n_free), lab])
        hm = hierarchical_moves(red, T)
        assert np.all(red.allowed[np.arange(red.n_free), hm])
        assert reduced_energy(red, hm, Dt) >= v - 1e-8


def test_binary_cut_bruteforce(rng):
    for _ in range(25):
        n = 7
        E = int(rng.integers(3, 12))
        eu, ev = rng.integers(0, n, E), rng.integers(0, n, E)
        keep = eu != ev
        eu, ev = eu[keep], ev[keep]
        U0, U1 = rng.random(n) * 3 - 1, rng.random(n) * 3 - 1
        A, Dd, B = rng.random(len(eu)), rng.random(len(eu)), rng.random(len(eu))
        C = A + Dd - B + rng.random(len(eu))

        def en(b):
            return np.sum(np.where(b, U1, U0)) + np.sum(
                np.select([~b[eu] & ~b[ev], ~b[eu] & b[ev], b[eu] & ~b[ev]], [A, B, C], Dd))

        best = min(en(np.array(bb, bool)) for bb in itertools.product([0, 1], repeat=n))
        assert abs(en(binary_cut(U0, U1, eu, ev, A, B, C, Dd)) - best) < 1e-9


@pytest.mark.parametrize("backend", ["gurobi", "cplex"])
def test_commercial_backends(rng, backend):
    pytest.importorskip({"gurobi": "gurobipy", "cplex": "cplex"}[backend])
    P = _problem(rng, 5, 3)
    red = reduce_problem(P)
    try:
        r = solve_relaxation(red, P.D, method="transport", backend=backend)
    except Exception as exc:  # licence problems
        pytest.skip(str(exc))
    r0 = solve_relaxation(red, P.D, method="transport")
    assert abs(r.objective - r0.objective) < 1e-6
