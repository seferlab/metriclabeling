"""Least-Squared-Distortion (LSD) metric approximation (Section 2.3, Eqs. 12-15).

Given a semimetric ``S`` on ``k`` items, find the metric ``M`` minimising

    sum_{i<j} (S_ij - M_ij)^2   s.t.   M_ij <= M_il + M_lj  for all triangles,

i.e. the Euclidean projection of ``S`` onto the metric cone.  The problem is a
strictly convex QP with ``3 * C(k, 3)`` constraints (352k for k = 90,
13.4M for k = 300).

Solvers
-------
``triangle_fixing`` (default)
    Exact Hildreth / Dykstra cyclic projection (the "triangle fixing"
    algorithm of Brickell, Dhillon, Sra & Tropp, SIAM J. Matrix Anal. 2008),
    JIT-compiled with numba.  Dual-feasible coordinate ascent; converges to the
    unique optimum.  Runs in seconds for k = 300.
``cvxpy``
    The full QP handed to CVXPY (Clarabel / OSQP / GUROBI / CPLEX ...).
    Practical up to k ~ 60.
``cutting_plane``
    Constraint generation: solve the QP with a subset of triangle constraints,
    add the violated ones, repeat (the "iterative approach" mentioned in
    Section 2.6).  Backends: ``cvxpy``, ``gurobi`` (gurobipy) or ``cplex``.
``gurobi`` / ``cplex``
    Shorthands for ``cutting_plane`` with that backend (the paper used CPLEX).
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Optional

import numpy as np

log = logging.getLogger(__name__)

try:  # numba is a hard dependency but keep the module importable without it
    from numba import njit
except Exception:  # pragma: no cover
    def njit(*args, **kwargs):
        def wrap(f):
            return f

        if args and callable(args[0]):
            return args[0]
        return wrap


@dataclass
class LSDResult:
    M: np.ndarray
    objective: float  # sum_{i<j} (S_ij - M_ij)^2
    violations_before: int
    max_violation_after: float
    iterations: int
    seconds: float
    method: str

    def as_dict(self) -> dict:
        return {
            "objective": self.objective,
            "violations_before": self.violations_before,
            "max_violation_after": self.max_violation_after,
            "iterations": self.iterations,
            "seconds": self.seconds,
            "method": self.method,
        }


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def lsd_objective(S: np.ndarray, M: np.ndarray) -> float:
    iu = np.triu_indices(S.shape[0], 1)
    return float(np.sum((S[iu] - M[iu]) ** 2))


def max_triangle_violation(M: np.ndarray) -> float:
    k = M.shape[0]
    worst = 0.0
    for l in range(k):
        V = M - (M[:, l][:, None] + M[l, :][None, :])
        V[l, :] = -np.inf
        V[:, l] = -np.inf
        np.fill_diagonal(V, -np.inf)
        worst = max(worst, float(V.max()) if k > 2 else 0.0)
    return max(worst, 0.0)


def violated_triangles(M: np.ndarray, tol: float, max_n: Optional[int] = None):
    """Return arrays (i, j, l, v): constraint ``M_ij <= M_il + M_lj`` violated by ``v``."""
    k = M.shape[0]
    I, J, L, V = [], [], [], []
    iu_mask = np.triu(np.ones((k, k), dtype=bool), 1)
    for l in range(k):
        Vl = M - (M[:, l][:, None] + M[l, :][None, :])
        mask = iu_mask & (Vl > tol)
        mask[l, :] = False
        mask[:, l] = False
        ii, jj = np.nonzero(mask)
        if len(ii):
            I.append(ii)
            J.append(jj)
            L.append(np.full(len(ii), l))
            V.append(Vl[ii, jj])
    if not I:
        e = np.zeros(0, dtype=np.int64)
        return e, e, e, np.zeros(0)
    I, J, L, V = (np.concatenate(x) for x in (I, J, L, V))
    if max_n is not None and len(V) > max_n:
        top = np.argpartition(-V, max_n - 1)[:max_n]
        I, J, L, V = I[top], J[top], L[top], V[top]
    return I.astype(np.int64), J.astype(np.int64), L.astype(np.int64), V


def _pair_index(i, j, k):
    i, j = np.minimum(i, j), np.maximum(i, j)
    return i * k - i * (i + 1) // 2 + (j - i - 1)


# --------------------------------------------------------------------------- #
# triangle fixing (Hildreth / Dykstra), numba kernel
# --------------------------------------------------------------------------- #
@njit(cache=True)
def _tf_sweep(M, Z):  # pragma: no cover - compiled
    k = M.shape[0]
    t = 0
    maxv = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            for l in range(j + 1, k):
                a = M[i, j]
                b = M[i, l]
                c = M[j, l]
                # constraint: a <= b + c
                v = a - b - c
                if v > maxv:
                    maxv = v
                d = v / 3.0
                z = Z[t]
                if d < -z:
                    d = -z
                if d != 0.0:
                    Z[t] = z + d
                    a -= d
                    b += d
                    c += d
                # constraint: b <= a + c
                v = b - a - c
                if v > maxv:
                    maxv = v
                d = v / 3.0
                z = Z[t + 1]
                if d < -z:
                    d = -z
                if d != 0.0:
                    Z[t + 1] = z + d
                    b -= d
                    a += d
                    c += d
                # constraint: c <= a + b
                v = c - a - b
                if v > maxv:
                    maxv = v
                d = v / 3.0
                z = Z[t + 2]
                if d < -z:
                    d = -z
                if d != 0.0:
                    Z[t + 2] = z + d
                    c -= d
                    a += d
                    b += d
                M[i, j] = a
                M[j, i] = a
                M[i, l] = b
                M[l, i] = b
                M[j, l] = c
                M[l, j] = c
                t += 3
    return maxv


def _triangle_fixing(S: np.ndarray, tol: float, max_iter: int, verbose: bool):
    k = S.shape[0]
    M = S.astype(np.float64).copy()
    ntri = k * (k - 1) * (k - 2) // 6
    Z = np.zeros(3 * ntri, dtype=np.float64)
    it = 0
    for it in range(1, max_iter + 1):
        prev = M.copy()
        _tf_sweep(M, Z)
        change = float(np.max(np.abs(M - prev)))
        if verbose and (it % 50 == 0 or it == 1):
            log.info("LSD sweep %d: change %.3e", it, change)
        if change < tol:
            break
    return M, it


# --------------------------------------------------------------------------- #
# QP backends for the full problem / cutting planes
# --------------------------------------------------------------------------- #
def _constraint_rows(I, J, L, k):
    """Rows of A in A m <= 0 for constraints m_IJ - m_IL - m_JL <= 0."""
    n = len(I)
    r = np.repeat(np.arange(n), 3)
    c = np.stack([_pair_index(I, J, k), _pair_index(I, L, k), _pair_index(J, L, k)], axis=1).ravel()
    v = np.tile(np.array([1.0, -1.0, -1.0]), n)
    return r, c, v


def _all_triangles(k):
    tri = np.array([(i, j, l) for i in range(k) for j in range(i + 1, k) for l in range(j + 1, k)], dtype=np.int64)
    if len(tri) == 0:
        e = np.zeros(0, dtype=np.int64)
        return e, e, e
    i, j, l = tri[:, 0], tri[:, 1], tri[:, 2]
    # three rotations: long side (i,j), (i,l), (j,l)
    I = np.concatenate([i, i, j])
    J = np.concatenate([j, l, l])
    L = np.concatenate([l, j, i])
    return I, J, L


def _solve_qp(s: np.ndarray, A, backend: str, solver: Optional[str], verbose: bool) -> np.ndarray:
    """min ||m - s||^2  s.t.  A m <= 0, m >= 0."""
    P = len(s)
    if backend == "cvxpy":
        import cvxpy as cp

        m = cp.Variable(P)
        cons = [m >= 0]
        if A.shape[0]:
            cons.append(A @ m <= 0)
        prob = cp.Problem(cp.Minimize(cp.sum_squares(m - s)), cons)
        kwargs = {}
        if solver is None:
            solver = "CLARABEL" if "CLARABEL" in cp.installed_solvers() else None
        prob.solve(solver=solver, verbose=verbose, **kwargs)
        return np.asarray(m.value, dtype=float)
    if backend == "gurobi":
        import gurobipy as gp

        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", int(verbose))
            env.start()
            with gp.Model(env=env) as model:
                m = model.addMVar(P, lb=0.0)
                model.setObjective(m @ m - 2 * s @ m, gp.GRB.MINIMIZE)
                if A.shape[0]:
                    model.addMConstr(A, m, "<", np.zeros(A.shape[0]))
                model.optimize()
                return np.asarray(m.X, dtype=float)
    if backend == "cplex":
        import cplex

        c = cplex.Cplex()
        if not verbose:
            c.set_log_stream(None)
            c.set_results_stream(None)
            c.set_warning_stream(None)
        c.objective.set_sense(c.objective.sense.minimize)
        c.variables.add(obj=list(-2.0 * s), lb=[0.0] * P)
        c.objective.set_quadratic([2.0] * P)  # 0.5 * x'Qx with Q = 2I
        if A.shape[0]:
            A = A.tocsr()
            lin = [
                cplex.SparsePair(ind=A.indices[A.indptr[r] : A.indptr[r + 1]].tolist(),
                                 val=A.data[A.indptr[r] : A.indptr[r + 1]].tolist())
                for r in range(A.shape[0])
            ]
            c.linear_constraints.add(lin_expr=lin, senses="L" * A.shape[0], rhs=[0.0] * A.shape[0])
        c.solve()
        return np.asarray(c.solution.get_values(), dtype=float)
    raise ValueError(f"unknown QP backend {backend!r}")


def _vec_to_mat(m: np.ndarray, k: int) -> np.ndarray:
    M = np.zeros((k, k))
    iu = np.triu_indices(k, 1)
    M[iu] = m
    return M + M.T


def _full_qp(S, backend, solver, verbose):
    import scipy.sparse as sp

    k = S.shape[0]
    iu = np.triu_indices(k, 1)
    I, J, L = _all_triangles(k)
    r, c, v = _constraint_rows(I, J, L, k)
    A = sp.csr_matrix((v, (r, c)), shape=(len(I), len(iu[0])))
    m = _solve_qp(S[iu], A, backend, solver, verbose)
    return _vec_to_mat(m, k), 1


def _cutting_plane(S, backend, solver, tol, max_iter, max_add, verbose):
    import scipy.sparse as sp

    k = S.shape[0]
    iu = np.triu_indices(k, 1)
    s = S[iu]
    P = len(s)
    rows_r, rows_c, rows_v = [], [], []
    seen = set()
    nrows = 0
    M = S.copy()
    it = 0
    for it in range(1, max_iter + 1):
        I, J, L, V = violated_triangles(M, tol, max_add)
        if len(I) == 0:
            break
        keys = _pair_index(I, J, k) * k + L
        new = [q for q, key in enumerate(keys.tolist()) if key not in seen]
        if not new:
            break
        seen.update(keys[new].tolist())
        r, c, v = _constraint_rows(I[new], J[new], L[new], k)
        rows_r.append(r + nrows)
        rows_c.append(c)
        rows_v.append(v)
        nrows += len(new)
        A = sp.csr_matrix(
            (np.concatenate(rows_v), (np.concatenate(rows_r), np.concatenate(rows_c))), shape=(nrows, P)
        )
        m = _solve_qp(s, A, backend, solver, verbose)
        M = _vec_to_mat(m, k)
        if verbose:
            log.info("cutting plane %d: %d constraints, max violation %.2e", it, nrows, V.max())
    return M, it


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def lsd(
    S: np.ndarray,
    method: str = "triangle_fixing",
    tol: float = 1e-8,
    max_iter: int = 20000,
    solver: Optional[str] = None,
    max_add: int = 200000,
    skip_if_metric: bool = True,
    verbose: bool = False,
) -> LSDResult:
    """Closest metric (in least squares) to a symmetric dissimilarity matrix ``S``.

    Parameters
    ----------
    method: ``triangle_fixing`` | ``cvxpy`` | ``cutting_plane`` | ``gurobi`` | ``cplex``
    solver: for ``cvxpy``/``cutting_plane``: the CVXPY solver name, or the
        cutting-plane backend (``cvxpy``/``gurobi``/``cplex``).
    tol: convergence tolerance (max entry change per sweep for triangle fixing,
        max allowed violation for cutting planes).
    """
    from .distances import triangle_violations

    S = np.asarray(S, dtype=np.float64)
    S = 0.5 * (S + S.T)
    np.fill_diagonal(S, 0.0)
    k = S.shape[0]
    t0 = time.perf_counter()
    nviol, _ = triangle_violations(S, tol=1e-12)
    if (nviol == 0 and skip_if_metric) or k < 3:
        return LSDResult(S.copy(), 0.0, nviol, 0.0, 0, time.perf_counter() - t0, "none")

    method = method.lower()
    if method == "triangle_fixing":
        M, it = _triangle_fixing(S, tol, max_iter, verbose)
    elif method == "cvxpy":
        M, it = _full_qp(S, "cvxpy", solver, verbose)
    elif method in ("cutting_plane", "gurobi", "cplex"):
        if method in ("gurobi", "cplex"):
            backend, qp_solver = method, None
        elif solver is not None and solver.lower() in ("cvxpy", "gurobi", "cplex"):
            backend, qp_solver = solver.lower(), None
        else:
            backend, qp_solver = "cvxpy", solver
        M, it = _cutting_plane(S, backend, qp_solver, max(tol, 1e-9), 1000, max_add, verbose)
    else:
        raise ValueError(f"unknown LSD method {method!r}")
    M = 0.5 * (M + M.T)
    np.fill_diagonal(M, 0.0)
    M = np.maximum(M, 0.0)
    return LSDResult(
        M=M,
        objective=lsd_objective(S, M),
        violations_before=nviol,
        max_violation_after=max_triangle_violation(M),
        iterations=it,
        seconds=time.perf_counter() - t0,
        method=method,
    )
