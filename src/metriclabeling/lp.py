"""Metric Labeling: problem reduction, LP relaxations, solvers and rounding.

Problem (Eq. 5)::

    Q(f) = sum_p c(p, f(p)) + sum_{(p,q) in E} w(p,q) d(f(p), f(q))

LP relaxations implemented
--------------------------
``transport``
    The Chekuri et al. (2005) LP of Eqs. (6)-(11): node variables x(u,i) and
    edge variables x(u,i,v,j).  Exact but ``k^2`` variables per edge.
``tree``
    The same LP when ``d`` is a tree (HST) metric.  The earth-mover cost on a
    tree is ``sum_t len(t) |x_u(T_t) - x_v(T_t)|``; we write it as a min-cost
    flow on the tree, ``O(k)`` variables per edge (Kleinberg & Tardos 1999).
    Solving it on a well-chosen FRT tree and rounding gives the
    ``O(log k)``-approximation used in the paper.

Nodes whose label is fixed (finite cost for one label only) are eliminated
first: their edges become linear terms of the free neighbours, so the LP only
contains the proteins to be predicted.

Solvers: HiGHS (via SciPy, default), Gurobi (``gurobipy``), CPLEX (``cplex``)
and, for the tree LP, a matrix-free primal-dual hybrid gradient method
(``pdhg``) that scales to the yeast network.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

from .hst import HST

log = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# problem containers
# --------------------------------------------------------------------------- #
@dataclass
class LabelingProblem:
    """A Metric Labeling instance.

    ``unary`` is ``(n, k)``; ``np.inf`` marks forbidden labels.  Edges are
    undirected (each listed once) with non-negative weights.
    """

    unary: np.ndarray
    eu: np.ndarray
    ev: np.ndarray
    ew: np.ndarray
    D: np.ndarray

    @property
    def n(self) -> int:
        return self.unary.shape[0]

    @property
    def k(self) -> int:
        return self.unary.shape[1]

    def energy(self, labels: np.ndarray, D: Optional[np.ndarray] = None) -> float:
        D = self.D if D is None else D
        un = self.unary[np.arange(self.n), labels]
        return float(un.sum() + np.sum(self.ew * D[labels[self.eu], labels[self.ev]]))

    def fractional_energy(self, X: np.ndarray, D: Optional[np.ndarray] = None) -> float:
        """Energy of a fractional labeling under the *independent* coupling (upper bound on LP cost)."""
        D = self.D if D is None else D
        U = np.where(np.isfinite(self.unary), self.unary, 0.0)
        un = float(np.sum(U * X))
        pair = np.einsum("ei,ij,ej->e", X[self.eu], D, X[self.ev])
        return un + float(np.sum(self.ew * pair))


@dataclass
class Reduced:
    """Instance restricted to the free nodes (fixed nodes folded into ``lin``)."""

    n: int
    k: int
    free: np.ndarray  # indices (into the original problem) of free nodes
    fixed: np.ndarray
    fixed_label: np.ndarray  # label of each fixed node (len n, -1 for free)
    lin: np.ndarray  # (n_free, k) linear cost incl. fixed neighbours; inf = forbidden
    fu: np.ndarray  # free-free edges in local indices
    fv: np.ndarray
    fw: np.ndarray
    const: float

    @property
    def n_free(self) -> int:
        return len(self.free)

    @property
    def allowed(self) -> np.ndarray:
        return np.isfinite(self.lin)


def reduce_problem(P: LabelingProblem, D: Optional[np.ndarray] = None) -> Reduced:
    """Eliminate nodes with a single allowed label."""
    D = P.D if D is None else D
    allowed = np.isfinite(P.unary)
    na = allowed.sum(axis=1)
    if np.any(na == 0):
        raise ValueError("every node needs at least one allowed label")
    is_fixed = na == 1
    fixed = np.flatnonzero(is_fixed)
    free = np.flatnonzero(~is_fixed)
    fixed_label = -np.ones(P.n, dtype=np.int64)
    fixed_label[fixed] = np.argmax(allowed[fixed], axis=1)
    pos = -np.ones(P.n, dtype=np.int64)
    pos[free] = np.arange(len(free))
    lin = P.unary[free].copy()
    const = float(P.unary[fixed, fixed_label[fixed]].sum())

    fu_, fv_ = is_fixed[P.eu], is_fixed[P.ev]
    both = fu_ & fv_
    const += float(np.sum(P.ew[both] * D[fixed_label[P.eu[both]], fixed_label[P.ev[both]]]))
    # free u, fixed v
    m = ~fu_ & fv_
    np.add.at(lin, pos[P.eu[m]], P.ew[m][:, None] * D[fixed_label[P.ev[m]]])
    m = fu_ & ~fv_
    np.add.at(lin, pos[P.ev[m]], P.ew[m][:, None] * D[fixed_label[P.eu[m]]])
    m = ~fu_ & ~fv_
    return Reduced(
        n=P.n,
        k=P.k,
        free=free,
        fixed=fixed,
        fixed_label=fixed_label,
        lin=lin,
        fu=pos[P.eu[m]],
        fv=pos[P.ev[m]],
        fw=P.ew[m].astype(float),
        const=const,
    )


# --------------------------------------------------------------------------- #
# generic LP solve:  min c'x  s.t.  A x = b,  x >= 0
# --------------------------------------------------------------------------- #
def solve_lp(c: np.ndarray, A: sp.csr_matrix, b: np.ndarray, backend: str = "highs", verbose: bool = False,
             time_limit: Optional[float] = None) -> np.ndarray:
    backend = backend.lower()
    if backend.startswith("highs"):
        from scipy.optimize import linprog

        method = {"highs": "highs", "highs-ipm": "highs-ipm", "highs-ds": "highs-ds"}.get(backend, "highs")
        opts = {"disp": verbose, "presolve": True}
        if time_limit:
            opts["time_limit"] = time_limit
        res = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method=method, options=opts)
        if res.status != 0:
            raise RuntimeError(f"HiGHS failed: {res.message}")
        return np.asarray(res.x)
    if backend == "gurobi":
        import gurobipy as gp

        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", int(verbose))
            env.start()
            with gp.Model(env=env) as m:
                x = m.addMVar(len(c), lb=0.0, obj=c)
                m.addMConstr(A, x, "=", b)
                if time_limit:
                    m.Params.TimeLimit = time_limit
                m.optimize()
                return np.asarray(x.X)
    if backend == "cplex":
        import cplex

        cp_ = cplex.Cplex()
        if not verbose:
            cp_.set_log_stream(None)
            cp_.set_results_stream(None)
            cp_.set_warning_stream(None)
        cp_.objective.set_sense(cp_.objective.sense.minimize)
        cp_.variables.add(obj=c.tolist(), lb=[0.0] * len(c))
        A = A.tocsr()
        rows = [
            cplex.SparsePair(ind=A.indices[A.indptr[r]: A.indptr[r + 1]].tolist(),
                             val=A.data[A.indptr[r]: A.indptr[r + 1]].tolist())
            for r in range(A.shape[0])
        ]
        cp_.linear_constraints.add(lin_expr=rows, senses="E" * A.shape[0], rhs=b.tolist())
        if time_limit:
            cp_.parameters.timelimit.set(time_limit)
        cp_.solve()
        return np.asarray(cp_.solution.get_values())
    raise ValueError(f"unknown LP backend {backend!r}")


# --------------------------------------------------------------------------- #
# LP builders
# --------------------------------------------------------------------------- #
def _x_index(allowed: np.ndarray) -> Tuple[np.ndarray, int]:
    xid = -np.ones(allowed.shape, dtype=np.int64)
    nz = np.flatnonzero(allowed.ravel())
    xid.ravel()[nz] = np.arange(len(nz))
    return xid, len(nz)


def _simplex_rows(xid: np.ndarray):
    r, c = np.nonzero(xid >= 0)
    return r, xid[r, c]


def build_transport_lp(red: Reduced, D: np.ndarray, nodes: np.ndarray, edges: np.ndarray):
    """Chekuri et al. LP (Eqs. 6-11) restricted to ``nodes`` / ``edges`` (local ids)."""
    allowed = red.allowed[nodes]
    loc = {int(u): i for i, u in enumerate(nodes)}
    xid, nx = _x_index(allowed)
    cost = [np.where(allowed, red.lin[nodes], 0.0)[allowed]]
    rows, cols, vals = [], [], []
    r, cidx = _simplex_rows(xid)
    rows.append(r)
    cols.append(cidx)
    vals.append(np.ones(len(r)))
    nrow = len(nodes)
    nvar = nx
    for e in edges:
        u, v, w = loc[int(red.fu[e])], loc[int(red.fv[e])], red.fw[e]
        Su = np.flatnonzero(allowed[u])
        Sv = np.flatnonzero(allowed[v])
        a, bb = len(Su), len(Sv)
        fid = nvar + np.arange(a * bb).reshape(a, bb)
        cost.append((w * D[np.ix_(Su, Sv)]).ravel())
        # sum_j f(i,j) - x(u,i) = 0
        rr = nrow + np.repeat(np.arange(a), bb)
        rows += [rr, nrow + np.arange(a)]
        cols += [fid.ravel(), xid[u, Su]]
        vals += [np.ones(a * bb), -np.ones(a)]
        nrow += a
        # sum_i f(i,j) - x(v,j) = 0
        rr = nrow + np.tile(np.arange(bb), a)
        rows += [rr, nrow + np.arange(bb)]
        cols += [fid.ravel(), xid[v, Sv]]
        vals += [np.ones(a * bb), -np.ones(bb)]
        nrow += bb
        nvar += a * bb
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nrow, nvar))
    b = np.zeros(nrow)
    b[: len(nodes)] = 1.0
    return np.concatenate(cost), A, b, xid


def build_tree_lp(red: Reduced, tree: HST, nodes: np.ndarray, edges: np.ndarray):
    """Earth-mover LP on a tree metric written as a min-cost flow per edge."""
    allowed = red.allowed[nodes]
    loc = -np.ones(red.n_free, dtype=np.int64)
    loc[nodes] = np.arange(len(nodes))
    xid, nx = _x_index(allowed)
    src, dst, ln = tree.arcs()
    nA, nN = len(src), tree.n_nodes
    E = len(edges)
    U = loc[red.fu[edges]]
    V = loc[red.fv[edges]]
    W = red.fw[edges]

    cost = np.concatenate([red.lin[nodes][allowed], (W[:, None] * ln[None, :]).ravel()])
    rows, cols, vals = [], [], []
    r, cidx = _simplex_rows(xid)
    rows.append(r)
    cols.append(cidx)
    vals.append(np.ones(len(r)))
    base_row = len(nodes)
    # flow conservation: for edge e and tree node h:
    #   sum_{a out of h} g - sum_{a into h} g - x(u, label h) + x(v, label h) = 0
    e_idx = np.repeat(np.arange(E), nA)
    a_idx = np.tile(np.arange(nA), E)
    gid = nx + e_idx * nA + a_idx
    rows += [base_row + e_idx * nN + src[a_idx], base_row + e_idx * nN + dst[a_idx]]
    cols += [gid, gid]
    vals += [np.ones(E * nA), -np.ones(E * nA)]
    # x terms at leaves
    leaf = tree.leaf
    for sign, X in ((-1.0, U), (1.0, V)):
        xe = xid[X]  # (E, k)
        ok = xe >= 0
        ee, ii = np.nonzero(ok)
        rows.append(base_row + ee * nN + leaf[ii])
        cols.append(xe[ee, ii])
        vals.append(np.full(len(ee), sign))
    nrow = base_row + E * nN
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(nrow, nx + E * nA))
    b = np.zeros(nrow)
    b[: len(nodes)] = 1.0
    return cost, A, b, xid


# --------------------------------------------------------------------------- #
# PDHG for the tree LP
# --------------------------------------------------------------------------- #
def project_simplex_masked(V: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Row-wise Euclidean projection onto the simplex restricted to ``mask``."""
    big = np.where(mask, V, -1e18)
    U = -np.sort(-big, axis=1)
    css = np.cumsum(np.where(U > -1e17, U, 0.0), axis=1) - 1.0
    ind = np.arange(1, V.shape[1] + 1)[None, :]
    cond = (U - css / ind > 0) & (U > -1e17)
    rho = cond.sum(axis=1)
    rho = np.maximum(rho, 1)
    theta = css[np.arange(V.shape[0]), rho - 1] / rho
    out = np.maximum(V - theta[:, None], 0.0)
    out[~mask] = 0.0
    return out


def pdhg_tree_lp(
    lin: np.ndarray,
    fu: np.ndarray,
    fv: np.ndarray,
    fw: np.ndarray,
    tree: HST,
    max_iter: int = 20000,
    tol: float = 1e-4,
    check_every: int = 50,
    X0: Optional[np.ndarray] = None,
    verbose: bool = False,
) -> Tuple[np.ndarray, Dict]:
    """Diagonally-preconditioned Chambolle-Pock for

        min_x  <L, x> + sum_e w_e sum_t len_t |(x_u - x_v) M_t|,   x_u in simplex(allowed_u)

    Returns the fractional labeling and convergence information.
    """
    n, k = lin.shape
    mask = np.isfinite(lin)
    L = np.where(mask, lin, 0.0)
    M = tree.incidence()[:, 1:]  # drop root (k, T)
    ln = tree.length[1:]
    T = M.shape[1]
    E = len(fu)
    if E == 0:
        X = np.zeros_like(L)
        best = np.where(mask, lin, np.inf).argmin(axis=1)
        X[np.arange(n), best] = 1.0
        return X, {"iterations": 0, "gap": 0.0}
    B = sp.csr_matrix(
        (np.r_[np.ones(E), -np.ones(E)], (np.r_[fu, fv], np.r_[np.arange(E), np.arange(E)])), shape=(n, E)
    )
    Bt = B.T.tocsr()
    box = fw[:, None] * ln[None, :]  # (E, T)
    deg = np.asarray(abs(B).sum(axis=1)).ravel()
    depth = M.sum(axis=1)  # number of non-root ancestors per label
    tau = 1.0 / np.maximum(deg[:, None] * depth.max(), 1.0)  # (n,1) per node
    sigma = 1.0 / (2.0 * np.maximum(M.sum(axis=0), 1.0))[None, :]  # (1, T)

    X = X0.copy() if X0 is not None else project_simplex_masked(-L, mask)
    Xbar = X.copy()
    Pd = np.zeros((E, T))
    info = {"iterations": 0, "gap": np.inf}
    t0 = time.perf_counter()
    for it in range(1, max_iter + 1):
        KX = (Bt @ Xbar) @ M  # (E, T)
        Pd = np.clip(Pd + sigma * KX, -box, box)
        G = B @ (Pd @ M.T)  # (n, k)
        Xn = project_simplex_masked(X - tau * (L + G), mask)
        Xbar = 2 * Xn - X
        X = Xn
        if it % check_every == 0 or it == max_iter:
            primal = float(np.sum(L * X) + np.sum(box * np.abs((Bt @ X) @ M)))
            dual_lin = L + B @ (Pd @ M.T)
            dual = float(np.sum(np.where(mask, dual_lin, np.inf).min(axis=1)))
            gap = (primal - dual) / max(1.0, abs(primal))
            info = {"iterations": it, "gap": gap, "primal": primal, "dual": dual,
                    "seconds": time.perf_counter() - t0}
            if verbose:
                log.info("pdhg it %d primal %.6f dual %.6f gap %.2e", it, primal, dual, gap)
            if gap < tol:
                break
    return X, info


# --------------------------------------------------------------------------- #
# relaxation driver
# --------------------------------------------------------------------------- #
@dataclass
class Relaxation:
    X: np.ndarray  # (n_free, k) fractional labeling of the free nodes
    method: str
    objective: float
    info: Dict = field(default_factory=dict)


def _components(red: Reduced) -> List[Tuple[np.ndarray, np.ndarray]]:
    nf = red.n_free
    if nf == 0:
        return []
    A = sp.coo_matrix((np.ones(len(red.fu)), (red.fu, red.fv)), shape=(nf, nf))
    ncomp, lab = connected_components(A, directed=False)
    order = np.argsort(lab, kind="stable")
    bounds = np.searchsorted(lab[order], np.arange(ncomp + 1))
    edge_comp = lab[red.fu]
    eorder = np.argsort(edge_comp, kind="stable")
    ebounds = np.searchsorted(edge_comp[eorder], np.arange(ncomp + 1))
    return [(order[bounds[c]: bounds[c + 1]], eorder[ebounds[c]: ebounds[c + 1]]) for c in range(ncomp)]


def _isolated_solution(lin_rows: np.ndarray) -> np.ndarray:
    """Optimal LP solution of nodes without free neighbours: split ties evenly."""
    best = np.min(lin_rows, axis=1, keepdims=True)
    tie = np.isclose(lin_rows, best, rtol=1e-12, atol=1e-12) & np.isfinite(lin_rows)
    return tie / tie.sum(axis=1, keepdims=True)


def tree_lp_objective(red: Reduced, X: np.ndarray, tree: HST) -> float:
    M = tree.incidence()[:, 1:]
    ln = tree.length[1:]
    lin = np.where(np.isfinite(red.lin), red.lin, 0.0)
    sep = np.abs((X[red.fu] - X[red.fv]) @ M) @ ln
    return float(np.sum(lin * X) + np.sum(red.fw * sep) + red.const)


def transport_lp_size(red: Reduced) -> int:
    na = red.allowed.sum(axis=1)
    return int(np.sum(na[red.fu] * na[red.fv]) + na.sum())


def solve_relaxation(
    red: Reduced,
    D: np.ndarray,
    method: str = "auto",
    tree: Optional[HST] = None,
    backend: str = "highs",
    transport_max_vars: int = 3_000_000,
    batch_vars: int = 400_000,
    pdhg_kwargs: Optional[Dict] = None,
    verbose: bool = False,
) -> Relaxation:
    """Solve the LP relaxation on the free nodes.

    method: ``transport`` (exact LP on ``D``), ``tree`` (LP on the tree metric
    of ``tree``) or ``auto`` (transport if it has at most
    ``transport_max_vars`` variables, otherwise tree).
    backend: ``highs`` / ``gurobi`` / ``cplex`` / ``pdhg`` (tree only).
    """
    t0 = time.perf_counter()
    if method == "auto":
        method = "transport" if transport_lp_size(red) <= transport_max_vars else "tree"
    if method == "tree" and tree is None:
        raise ValueError("tree LP requires an HST")
    X = np.zeros((red.n_free, red.k))
    if red.n_free == 0:
        return Relaxation(X, method, red.const, {"seconds": 0.0})

    if method == "tree" and backend == "pdhg":
        X, info = pdhg_tree_lp(red.lin, red.fu, red.fv, red.fw, tree, verbose=verbose, **(pdhg_kwargs or {}))
        obj = tree_lp_objective(red, X, tree)
        info.update(seconds=time.perf_counter() - t0, backend="pdhg")
        return Relaxation(X, method, obj, info)

    comps = _components(red)
    # isolated nodes in closed form
    singles = np.array([c[0][0] for c in comps if len(c[0]) == 1], dtype=np.int64)
    if len(singles):
        X[singles] = _isolated_solution(red.lin[singles])
    multi = [c for c in comps if len(c[0]) > 1]
    na = red.allowed.sum(axis=1)
    # batch small components together
    batches: List[Tuple[List[np.ndarray], List[np.ndarray]]] = []
    cur_n, cur_e, cur_size = [], [], 0
    for nodes, edges in multi:
        if method == "transport":
            size = int(np.sum(na[red.fu[edges]] * na[red.fv[edges]]) + na[nodes].sum())
        else:
            size = int(len(edges) * 2 * (tree.n_nodes - 1) + na[nodes].sum())
        if cur_n and cur_size + size > batch_vars:
            batches.append((cur_n, cur_e))
            cur_n, cur_e, cur_size = [], [], 0
        cur_n.append(nodes)
        cur_e.append(edges)
        cur_size += size
    if cur_n:
        batches.append((cur_n, cur_e))
    nvars = 0
    obj = red.const
    if len(singles):
        obj += float(np.sum(np.where(np.isfinite(red.lin[singles]), red.lin[singles], 0.0) * X[singles]))
    for bn, be in batches:
        nodes = np.concatenate(bn)
        edges = np.concatenate(be)
        if method == "transport":
            c, A, b, xid = build_transport_lp(red, D, nodes, edges)
        else:
            c, A, b, xid = build_tree_lp(red, tree, nodes, edges)
        nvars += len(c)
        sol = solve_lp(c, A, b, backend=backend, verbose=verbose)
        obj += float(c @ sol)
        Xb = np.zeros((len(nodes), red.k))
        ok = xid >= 0
        Xb[ok] = sol[xid[ok]]
        X[nodes] = Xb
    X = np.clip(X, 0.0, None)
    X /= np.maximum(X.sum(axis=1, keepdims=True), 1e-12)
    info = {"seconds": time.perf_counter() - t0, "backend": backend, "lp_vars": nvars,
            "components": len(comps), "batches": len(batches)}
    return Relaxation(X, method, obj, info)


# --------------------------------------------------------------------------- #
# rounding on an HST (Kleinberg-Tardos / Chekuri et al.)
# --------------------------------------------------------------------------- #
def round_on_tree(X: np.ndarray, tree: HST, rng: Optional[np.random.Generator] = None,
                  allowed: Optional[np.ndarray] = None) -> np.ndarray:
    """Hierarchical randomized rounding of a fractional labeling on an HST.

    Starting at the root, the nodes currently assigned to a tree node ``t``
    are split among its children with the Kleinberg-Tardos uniform-metric
    rounding: repeatedly pick a random child ``c`` and threshold ``alpha``
    and send every still-unassigned node with ``y_u(c) / y_u(t) >= alpha``
    to ``c``.  ``y_u(t)`` is the fractional mass of ``u`` under ``t``.
    """
    rng = np.random.default_rng(rng)
    n, k = X.shape
    Minc = tree.incidence()  # (k, N)
    Y = X @ Minc  # (n, N)
    assign = np.zeros(n, dtype=np.int64)  # current tree node (root = 0)
    for t in range(tree.n_nodes):
        ch = tree.children[t]
        if not ch:
            continue
        members = np.flatnonzero(assign == t)
        if len(members) == 0:
            continue
        if len(ch) == 1:
            assign[members] = ch[0]
            continue
        yt = Y[members, t]
        R = Y[np.ix_(members, ch)] / np.maximum(yt[:, None], 1e-15)
        R = np.where(yt[:, None] > 1e-15, R, 1.0 / len(ch))
        if allowed is not None:  # never route a node into a subtree without allowed labels
            okc = (allowed[members].astype(float) @ Minc[:, ch]) > 0
            R = np.where(okc, R, 0.0)
        todo = np.ones(len(members), dtype=bool)
        guard = 0
        while todo.any():
            c = rng.integers(len(ch))
            a = rng.uniform(1e-12, 1.0)
            hit = todo & (R[:, c] >= a)
            assign[members[hit]] = ch[c]
            todo &= ~hit
            guard += 1
            if guard > 100000:  # numerical safety: send leftovers to the heaviest child
                rest = members[todo]
                assign[rest] = np.asarray(ch)[np.argmax(R[todo], axis=1)]
                break
    labels = tree.label_of_leaf[assign]
    if np.any(labels < 0):
        raise RuntimeError("rounding ended at a non-leaf node")
    return labels


def emd_tree_cost(red: Reduced, X: np.ndarray, tree: HST) -> float:
    """Fractional cost of ``X`` measured with the tree metric (Chekuri's tree-selection criterion)."""
    return tree_lp_objective(red, X, tree)


def full_labels(red: Reduced, free_labels: np.ndarray) -> np.ndarray:
    lab = red.fixed_label.copy()
    lab[red.free] = free_labels
    return lab


def full_marginals(red: Reduced, X: np.ndarray) -> np.ndarray:
    out = np.zeros((red.n, red.k))
    out[red.fixed, red.fixed_label[red.fixed]] = 1.0
    out[red.free] = X
    return out


def solve_ilp(P: LabelingProblem, backend: str = "highs", time_limit: Optional[float] = None) -> np.ndarray:
    """Exact Metric Labeling via the integer program (6)-(11).  Small instances only."""
    red = reduce_problem(P)
    nodes = np.arange(red.n_free)
    edges = np.arange(len(red.fu))
    c, A, b, xid = build_transport_lp(red, P.D, nodes, edges)
    if backend.startswith("highs"):
        from scipy.optimize import milp, LinearConstraint, Bounds

        res = milp(c, constraints=LinearConstraint(A, b, b), integrality=np.ones(len(c)),
                   bounds=Bounds(0, 1), options={"time_limit": time_limit} if time_limit else None)
        if res.x is None:
            raise RuntimeError(res.message)
        sol = res.x
    elif backend == "gurobi":
        import gurobipy as gp

        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            with gp.Model(env=env) as m:
                x = m.addMVar(len(c), vtype=gp.GRB.BINARY, obj=c)
                m.addMConstr(A, x, "=", b)
                m.optimize()
                sol = x.X
    else:
        raise ValueError(backend)
    Xf = np.zeros((red.n_free, red.k))
    ok = xid >= 0
    Xf[ok] = sol[xid[ok]]
    return full_labels(red, Xf.argmax(axis=1))
