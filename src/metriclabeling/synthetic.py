"""Synthetic GO-like benchmark used for tests and as an end-to-end sanity check.

A random ontology DAG is generated, ``k`` mutually non-nested terms are used as
labels, and a protein network is planted so that interacting proteins carry
labels that are *close in the ontology* (edge probability decays with the GO
shortest-path distance between their labels).  On such data Metric Labeling
with a GO distance should beat the uniform-metric (multiway cut) model.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sp

from .dataset import Dataset, build_go_info, label_matrix, largest_component
from .io.obo import Ontology

log = logging.getLogger(__name__)


def random_ontology(depth: int = 4, branching: int = 3, extra_parent_prob: float = 0.15,
                    rng: Optional[np.random.Generator] = None) -> Ontology:
    rng = np.random.default_rng(rng)
    ids, parents, level = ["GO:0008150"], [[]], [0]
    frontier = [0]
    by_level = {0: [0]}
    for d in range(1, depth + 1):
        nxt = []
        for p in frontier:
            for _ in range(branching):
                i = len(ids)
                ids.append(f"GO:{9000000 + i:07d}")
                ps = [p]
                if d > 1 and rng.random() < extra_parent_prob:
                    q = int(rng.choice(by_level[d - 1]))
                    if q != p:
                        ps.append(q)
                parents.append(ps)
                level.append(d)
                nxt.append(i)
        by_level[d] = nxt
        frontier = nxt
    return Ontology(ids, [f"term {i}" for i in range(len(ids))], ["biological_process"] * len(ids), parents)


def synthetic_dataset(n: int = 400, k: int = 20, avg_degree: float = 8.0, decay: float = 1.2,
                      second_label_prob: float = 0.3, noise_edges: float = 0.15, seed: int = 0) -> Dataset:
    rng = np.random.default_rng(seed)
    onto = random_ontology(rng=rng)
    leaves = [i for i in range(len(onto)) if not onto.children[i]]
    selected = sorted(rng.choice(leaves, size=k, replace=False).tolist())
    # GO shortest paths between selected terms
    from scipy.sparse.csgraph import shortest_path

    Dgo = shortest_path(onto.undirected_adjacency(), unweighted=True, directed=False, indices=selected)[:, selected]
    home = rng.integers(k, size=n)
    labels: List[List[int]] = []
    for p in range(n):
        ls = [int(home[p])]
        if rng.random() < second_label_prob:
            w = np.exp(-decay * Dgo[home[p]])
            w[home[p]] = 0
            ls.append(int(rng.choice(k, p=w / w.sum())))
        labels.append(ls)
    # planted network
    P = np.exp(-decay * Dgo[home][:, home])
    np.fill_diagonal(P, 0)
    P *= (avg_degree * (1 - noise_edges) * n / 2) / P[np.triu_indices(n, 1)].sum()
    U = rng.random((n, n)) < P
    m_noise = int(noise_edges * avg_degree * n / 2)
    U[rng.integers(n, size=m_noise), rng.integers(n, size=m_noise)] = True
    U = np.triu(U, 1)
    U = U | U.T
    A = sp.csr_matrix(U.astype(float))
    lcc = largest_component(A)
    A = A[lcc][:, lcc].tocsr()
    protein_terms = [[selected[j] for j in labels[p]] for p in lcc]
    # hide the labels of 10% of proteins entirely (unannotated proteins)
    unannot = rng.random(len(lcc)) < 0.1
    protein_terms = [[] if unannot[i] else t for i, t in enumerate(protein_terms)]
    Y = label_matrix(onto, selected, protein_terms)
    go = build_go_info(onto, selected, protein_terms)
    return Dataset(f"synthetic-{n}-{k}", [f"P{i}" for i in lcc], A, [onto.ids[t] for t in selected],
                   [onto.names[t] for t in selected], Y, go, {"species": "synthetic", "go_diameter": go.go_diameter})


DEFAULT_METHODS = [
    {"label": "Majority", "method": "majority"},
    {"label": "Neighbour", "method": "neighborhood"},
    {"label": "FunctionalFlow", "method": "functional_flow"},
    {"label": "MRF", "method": "mrf"},
    {"label": "GenMulticut", "method": "genmulticut", "runs": 10},
    {"label": "Potts ML", "method": "metric_labeling", "distance": "uniform"},
    {"label": "SP ML", "method": "metric_labeling", "distance": "sp"},
    {"label": "LCA LSD", "method": "metric_labeling", "distance": "lca"},
    {"label": "Lin LSD", "method": "metric_labeling", "distance": "lin"},
    {"label": "LCA 0.3 LSD", "method": "metric_labeling", "distance": "lca", "alpha": 0.3},
    {"label": "LCA 0.3 SemiMAP", "method": "semimap", "distance": "lca", "alpha": 0.3, "runs": 10},
]


def run_synthetic_benchmark(out: Path, seed: int = 0, n: int = 400, k: int = 20,
                            methods: Optional[List[Dict]] = None, folds: int = 5) -> Dict[str, float]:
    import warnings

    from .evaluation import aggregate, evaluate_fold
    from .folds import make_folds
    from .methods import make_method
    from .plotting import plot_roc

    warnings.filterwarnings("ignore")
    ds = synthetic_dataset(n=n, k=k, seed=seed)
    res: Dict[str, List[Dict]] = {}
    for fold in make_folds(ds, folds, seed=seed):
        for spec in methods or DEFAULT_METHODS:
            m = make_method({**spec, "seed": seed})
            S = m.predict(fold)
            res.setdefault(spec["label"], []).append(evaluate_fold(ds.Y, S, fold.test))
    agg = {lab: aggregate(r) for lab, r in res.items()}
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "curves.json").write_text(json.dumps({lab: {"fpr": a["fpr"], "tpr": a["tpr"]} for lab, a in agg.items()}))
    summary = {lab: round(a["auc"], 4) for lab, a in agg.items()}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    plot_roc({lab: a for lab, a in agg.items()}, out / "synthetic_roc.png", f"Synthetic benchmark (n={ds.n}, k={k})")
    for lab, v in summary.items():
        print(f"{lab:18s} AUC {v:.3f}")
    return summary
