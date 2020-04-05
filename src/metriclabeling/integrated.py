"""Integrated S. cerevisiae network (Section 2.5).

Sources (each a binary network over canonical SGD ids):

* ``ppi``          BioGRID physical interactions
* ``complex``      co-membership in the TAP complexes of Gavin et al. (2002)
                   and Ho et al. (2002).  Complexes are reconstructed from the
                   BioGRID *Affinity Capture-MS* records of those papers
                   (bait + its preys) and expanded to all member pairs.  A
                   file with one complex per line can be given instead.
* ``expression``   |Pearson r| > 0.8 in the Gasch et al. (2000) compendium
* ``localization`` same GFP localisation in Huh et al. (2003), read from the
                   GAF cellular-component annotations citing PMID:14562095
                   (or from a two-column ``gene<TAB>location`` file).
                   Locations with more than ``max_location_size`` proteins
                   (e.g. *cytoplasm*) are skipped, otherwise they would add
                   millions of uninformative edges.

Edges are combined with noisy-or, ``w(u,v) = 1 - prod_i (1 - r_i)``, where
the reliability ``r_i`` of source ``i`` is the fraction of its edges (between
two training proteins) that connect proteins sharing a function.  Because
``r_i`` uses training labels only, weights are recomputed for every fold.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp

from .dataset import (
    REPO_ROOT,
    Dataset,
    build_go_info,
    edges_to_adjacency,
    label_matrix,
    node_annotations,
    largest_component,
    read_network_and_annotations,
    select_terms,
    species_files,
)
from .io import IdMapper, read_biogrid, read_gaf
from .io.files import open_text

log = logging.getLogger(__name__)

GAVIN_2002 = "11805826"
HO_2002 = "11805837"
HUH_2003 = "PMID:14562095"


def complexes_from_biogrid(biogrid_path: Path, mapper: IdMapper,
                           pubmeds: Iterable[str] = (GAVIN_2002, HO_2002)) -> List[List[str]]:
    inter = read_biogrid(biogrid_path, physical_only=True, systems=["Affinity Capture-MS"], pubmeds=list(pubmeds))
    groups: Dict[Tuple[str, str], set] = defaultdict(set)
    for it in inter:
        a = mapper.map_any(it.names_a)
        b = mapper.map_any(it.names_b)
        if a is None or b is None:
            continue
        key = (it.pubmed, a)
        groups[key].add(a)
        groups[key].add(b)
    return [sorted(g) for g in groups.values() if len(g) > 1]


def read_complex_file(path: Path, mapper: IdMapper) -> List[List[str]]:
    out = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            mem = [mapper.map(x) for x in line.replace(",", " ").split()]
            mem = sorted({m for m in mem if m})
            if len(mem) > 1:
                out.append(mem)
    return out


def clique_edges(groups: Iterable[List[str]]) -> List[Tuple[str, str]]:
    edges = set()
    for g in groups:
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                edges.add((g[i], g[j]))
    return sorted(edges)


def _to_float(x: str) -> float:
    try:
        return float(x)
    except ValueError:
        return float("nan")


def expression_edges(path: Path, mapper: IdMapper, threshold: float = 0.8) -> List[Tuple[str, str]]:
    """Pairs of genes with |Pearson correlation| > threshold (missing values mean-imputed)."""
    names, rows = [], []
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        start = 2 if len(header) > 1 and header[1].upper() == "NAME" else 1
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if not c or not c[0] or c[0].upper() in ("EWEIGHT",):
                continue
            g = mapper.map(c[0])
            if g is None:
                continue
            vals = np.array([_to_float(x) for x in c[start:]], dtype=float)
            if np.sum(np.isfinite(vals)) < 10:
                continue
            names.append(g)
            rows.append(vals)
    X = np.vstack(rows).astype(np.float32)
    mu = np.nanmean(X, axis=1, keepdims=True)
    X = np.where(np.isfinite(X), X - mu, 0.0)
    X /= np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    edges = []
    for s in range(0, len(names), 1000):
        C = X[s : s + 1000] @ X.T
        ii, jj = np.nonzero(np.abs(C) > threshold)
        for i, j in zip(ii + s, jj):
            if i < j and names[i] != names[j]:
                edges.append((names[i], names[j]))
    return edges


def localization_groups(gaf_path: Path, max_size: int = 200,
                        location_file: Optional[Path] = None, mapper: Optional[IdMapper] = None) -> List[List[str]]:
    groups: Dict[str, set] = defaultdict(set)
    if location_file is not None:
        with open_text(location_file) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 2 and not line.startswith("#"):
                    g = mapper.map(c[0]) if mapper else c[0]
                    if g:
                        for loc in c[1].split(","):
                            groups[loc.strip()].add(g)
    else:
        for r in read_gaf(gaf_path, aspects=("C",), exclude_evidence=(), references=[HUH_2003]):
            groups[r.go_id].add(r.canonical_id)
    out = []
    for loc, mem in groups.items():
        if 1 < len(mem) <= max_size:
            out.append(sorted(mem))
        elif len(mem) > max_size:
            log.info("localization: skipping %s (%d proteins)", loc, len(mem))
    return out


def build_integrated_dataset(
    species: str = "yeast",
    n_terms: int = 90,
    source: str = "original",
    data_dir: Path = REPO_ROOT / "data",
    min_count: int = 3,
    max_frac: float = 0.05,
    expression_threshold: float = 0.8,
    max_location_size: int = 200,
    complex_file: Optional[Path] = None,
    location_file: Optional[Path] = None,
    terms_file: Optional[Path] = None,
) -> Dataset:
    if species != "yeast":
        raise ValueError("the integrated network is defined for S. cerevisiae only")
    files = species_files(species, source, data_dir)
    onto, diam, mapper, ppi_edges, bp, _ = read_network_and_annotations(species, source, data_dir)

    sources: Dict[str, List[Tuple[str, str]]] = {"ppi": ppi_edges}
    groups = read_complex_file(complex_file, mapper) if complex_file else complexes_from_biogrid(files["biogrid"], mapper)
    sources["complex"] = clique_edges(groups)
    if files.get("expression") and Path(files["expression"]).exists():
        sources["expression"] = expression_edges(files["expression"], mapper, expression_threshold)
    else:
        log.warning("expression file not found; integrated network built without it")
    sources["localization"] = clique_edges(localization_groups(files["gaf"], max_location_size, location_file, mapper))
    for k, e in sources.items():
        log.info("integrated source %s: %d edges", k, len(e))

    all_edges = [e for es in sources.values() for e in es]
    nodes_all, A_all = edges_to_adjacency(all_edges)
    lcc = largest_component(A_all)
    nodes = [nodes_all[i] for i in lcc]
    A = A_all[lcc][:, lcc].tocsr()
    A.data[:] = 1.0
    mats = {}
    for k, e in sources.items():
        _, M = edges_to_adjacency(e, nodes=nodes)
        M.data[:] = 1.0
        mats[k] = M.tocsr()

    protein_terms = [node_annotations(p, bp) for p in nodes]
    from collections import Counter

    counts: Counter = Counter()
    n_ann = 0
    for ts in protein_terms:
        if ts:
            n_ann += 1
            counts.update(onto.propagate(ts).tolist())
    if terms_file is not None:
        wanted = [l.split()[0] for l in Path(terms_file).read_text().splitlines() if l.strip() and not l.startswith("#")]
        selected = [t for t in (onto.resolve(w) for w in wanted) if t is not None]
    else:
        selected = select_terms(onto, counts, n_terms, n_ann, min_count=min_count, max_frac=max_frac)
    Y = label_matrix(onto, selected, protein_terms)
    go = build_go_info(onto, selected, protein_terms, go_diameter=diam)
    meta = {
        "species": species,
        "source": source,
        "network": "integrated",
        "go_diameter": diam,
        "sources": mats,
        **{f"edges_{k}": int(v.nnz // 2) for k, v in mats.items()},
    }
    return Dataset(f"{species}-integrated-{len(selected)}", nodes, A, [onto.ids[t] for t in selected],
                   [onto.names[t] for t in selected], Y, go, meta)


def source_reliability(M: sp.csr_matrix, Y_train: np.ndarray, known: np.ndarray) -> float:
    """Fraction of the source's edges between two training proteins that join proteins sharing a function."""
    T = sp.triu(M, k=1).tocoo()
    m = known[T.row] & known[T.col]
    if not m.any():
        return 0.0
    share = (Y_train[T.row[m]] & Y_train[T.col[m]]).any(axis=1)
    return float(share.mean())


def noisy_or_adjacency(sources: Dict[str, sp.csr_matrix], Y_train: np.ndarray, known: np.ndarray):
    """Noisy-or combination of the sources with per-fold reliabilities."""
    rel = {k: source_reliability(M, Y_train, known) for k, M in sources.items()}
    n = next(iter(sources.values())).shape[0]
    logq = sp.csr_matrix((n, n))
    for k, M in sources.items():
        r = min(max(rel[k], 1e-6), 1 - 1e-9)
        logq = logq + M * np.log1p(-r)
    W = logq.tocsr()
    W.data = -np.expm1(W.data)  # 1 - prod(1 - r_i)
    W.eliminate_zeros()
    return W, rel
