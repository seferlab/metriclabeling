"""Building the protein networks + GO label sets used in the experiments.

A :class:`Dataset` bundles

* a weighted, undirected protein network (largest connected component),
* the ``k`` selected GO Biological-Process terms and the binary label matrix
  ``Y`` (true-path propagated, so a protein annotated to a descendant of a
  selected term carries that term),
* :class:`GOInfo`, the GO-structure quantities needed by the distances of
  Section 2.4 (shortest paths, upward distances, ancestor sets, per-protein
  propagated annotations for Lin's information content).
"""

from __future__ import annotations

import logging
import pickle
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components, shortest_path

from .io import IdMapper, read_biogrid, read_gaf, read_obo
from .io.obo import Ontology

log = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]


# --------------------------------------------------------------------------- #
# containers
# --------------------------------------------------------------------------- #
@dataclass
class GOInfo:
    """GO-structure quantities restricted to what the distances need."""

    term_ids: List[str]
    term_names: List[str]
    sp_len: np.ndarray  # (k,k) undirected shortest-path lengths in the GO DAG
    go_diameter: float  # diameter of the undirected GO (sub-)graph
    anc_ids: List[str]  # union of ancestors of the selected terms (incl. themselves)
    up: np.ndarray  # (k, |A|) upward distance term -> ancestor (inf if not an ancestor)
    anc_depth: np.ndarray  # (|A|,) shortest root -> ancestor distance
    P: sp.csr_matrix  # (n, |A|) protein is annotated (propagated) to ancestor a
    has_bp: np.ndarray  # (n,) protein has at least one (non-IEA) BP annotation

    @property
    def is_ancestor(self) -> np.ndarray:
        return np.isfinite(self.up)


@dataclass
class Dataset:
    name: str
    proteins: List[str]
    adj: sp.csr_matrix
    terms: List[str]
    term_names: List[str]
    Y: np.ndarray  # (n,k) bool
    go: GOInfo
    meta: Dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.proteins)

    @property
    def k(self) -> int:
        return len(self.terms)

    @property
    def labeled(self) -> np.ndarray:
        """Indices of proteins with at least one selected label."""
        return np.flatnonzero(self.Y.any(axis=1))

    @property
    def edges(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        A = sp.triu(self.adj, k=1).tocoo()
        return A.row.astype(np.int64), A.col.astype(np.int64), A.data.astype(np.float64)

    def summary(self) -> Dict:
        u, v, _ = self.edges
        return {
            "name": self.name,
            "proteins_in_lcc": self.n,
            "interactions_in_lcc": int(len(u)),
            "labeled_proteins": int(len(self.labeled)),
            "terms": self.k,
            "annotations": int(self.Y.sum()),
            **{k: v for k, v in self.meta.items() if isinstance(v, (int, float, str))},
        }

    def save(self, path: Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as fh:
            pickle.dump(self, fh, protocol=4)

    @staticmethod
    def load(path: Path) -> "Dataset":
        with open(path, "rb") as fh:
            return pickle.load(fh)

    def with_adjacency(self, adj: sp.csr_matrix, name: Optional[str] = None) -> "Dataset":
        return Dataset(
            name or self.name, self.proteins, adj.tocsr(), self.terms, self.term_names, self.Y, self.go, dict(self.meta)
        )


# --------------------------------------------------------------------------- #
# graph helpers
# --------------------------------------------------------------------------- #
def edges_to_adjacency(
    edges: Iterable[Tuple[str, str]],
    weights: Optional[Iterable[float]] = None,
    nodes: Optional[Sequence[str]] = None,
    combine: str = "max",
) -> Tuple[List[str], sp.csr_matrix]:
    """Build a symmetric adjacency matrix; duplicate edges are merged."""
    edges = list(edges)
    if weights is None:
        weights = [1.0] * len(edges)
    weights = list(weights)
    if nodes is None:
        nodes = sorted({x for e in edges for x in e})
    index = {p: i for i, p in enumerate(nodes)}
    best: Dict[Tuple[int, int], float] = {}
    for (a, b), w in zip(edges, weights):
        if a == b or a not in index or b not in index:
            continue
        i, j = index[a], index[b]
        if i > j:
            i, j = j, i
        if combine == "max":
            best[(i, j)] = max(best.get((i, j), 0.0), float(w))
        else:
            best[(i, j)] = best.get((i, j), 0.0) + float(w)
    n = len(nodes)
    if not best:
        return list(nodes), sp.csr_matrix((n, n))
    ij = np.array(list(best.keys()), dtype=np.int64)
    w = np.array(list(best.values()), dtype=np.float64)
    A = sp.coo_matrix((np.r_[w, w], (np.r_[ij[:, 0], ij[:, 1]], np.r_[ij[:, 1], ij[:, 0]])), shape=(n, n))
    return list(nodes), A.tocsr()


def largest_component(adj: sp.csr_matrix) -> np.ndarray:
    """Indices of the nodes of the largest connected component."""
    if adj.shape[0] == 0:
        return np.zeros(0, dtype=np.int64)
    _, comp = connected_components(adj, directed=False)
    big = np.bincount(comp).argmax()
    return np.flatnonzero(comp == big)


def graph_diameter(adj: sp.csr_matrix) -> float:
    """Exact diameter of a connected, unweighted graph (chunked BFS)."""
    n = adj.shape[0]
    if n <= 1:
        return 0.0
    try:
        import networkx as nx

        G = nx.from_scipy_sparse_array(adj)
        if nx.is_connected(G):
            return float(nx.diameter(G, usebounds=True))
    except Exception:  # pragma: no cover - fallback
        pass
    best = 0.0
    for start in range(0, n, 512):
        idx = np.arange(start, min(n, start + 512))
        D = shortest_path(adj, unweighted=True, directed=False, indices=idx)
        D[~np.isfinite(D)] = 0
        best = max(best, float(D.max()))
    return best


# --------------------------------------------------------------------------- #
# GO term selection
# --------------------------------------------------------------------------- #
def _descendants(onto: Ontology, i: int, cache: Dict[int, set]) -> set:
    if i in cache:
        return cache[i]
    seen = {i}
    stack = [i]
    while stack:
        u = stack.pop()
        for c in onto.children[u]:
            if c not in seen:
                seen.add(c)
                stack.append(c)
    cache[i] = seen
    return seen


def select_terms(
    onto: Ontology,
    counts: Dict[int, int],
    k: int,
    n_annotated: int,
    min_count: int = 3,
    max_frac: float = 0.05,
    max_count: Optional[int] = None,
) -> List[int]:
    """Greedy selection of the ``k`` most frequent, mutually non-nested terms.

    The paper selects "those annotations that are seen more than others and
    also that are not parent of each other".  Very general terms (e.g.
    *cellular process*) would block almost the whole ontology, so terms
    annotating more than ``max_frac`` of the annotated proteins (or more than
    ``max_count``) are skipped.  The root is never selected.
    """
    if max_count is None:
        max_count = max(int(round(max_frac * n_annotated)), 2 * min_count)
    root = onto.root
    cand = sorted(
        (t for t, c in counts.items() if t != root and min_count <= c <= max_count),
        key=lambda t: (-counts[t], onto.ids[t]),
    )
    chosen: List[int] = []
    blocked: set = set()
    desc_cache: Dict[int, set] = {}
    for t in cand:
        if t in blocked:
            continue
        chosen.append(t)
        blocked.update(onto.ancestors(t).tolist())
        blocked.update(_descendants(onto, t, desc_cache))
        if len(chosen) == k:
            break
    if len(chosen) < k:
        log.warning("only %d non-nested terms satisfy the frequency filter (requested %d)", len(chosen), k)
    return chosen


# --------------------------------------------------------------------------- #
# GO info
# --------------------------------------------------------------------------- #
def build_go_info(
    onto: Ontology,
    selected: Sequence[int],
    protein_terms: Sequence[Sequence[int]],
    go_diameter: Optional[float] = None,
) -> GOInfo:
    """Pre-compute the GO quantities needed by :mod:`metriclabeling.distances`."""
    k = len(selected)
    root = onto.root
    # ancestors of the selected terms
    anc_union = sorted(set().union(*[set(onto.ancestors(t).tolist()) for t in selected]))
    a_index = {a: j for j, a in enumerate(anc_union)}
    up = np.full((k, len(anc_union)), np.inf)
    for i, t in enumerate(selected):
        for a, d in onto.up_distances(t).items():
            up[i, a_index[a]] = d
    anc_depth = np.array([onto.up_distances(a)[root] for a in anc_union], dtype=float)

    # undirected shortest paths among selected terms
    U = onto.undirected_adjacency()
    D = shortest_path(U, unweighted=True, directed=False, indices=np.asarray(selected))
    sp_len = D[:, np.asarray(selected)]
    if go_diameter is None:
        go_diameter = graph_diameter(U)

    # propagated annotations restricted to the ancestor set
    rows, cols = [], []
    has_bp = np.zeros(len(protein_terms), dtype=bool)
    for p, ts in enumerate(protein_terms):
        if len(ts) == 0:
            continue
        has_bp[p] = True
        for a in onto.propagate(ts).tolist():
            j = a_index.get(a)
            if j is not None:
                rows.append(p)
                cols.append(j)
    P = sp.csr_matrix(
        (np.ones(len(rows), dtype=bool), (rows, cols)), shape=(len(protein_terms), len(anc_union)), dtype=bool
    )
    return GOInfo(
        term_ids=[onto.ids[t] for t in selected],
        term_names=[onto.names[t] for t in selected],
        sp_len=sp_len,
        go_diameter=float(go_diameter),
        anc_ids=[onto.ids[a] for a in anc_union],
        up=up,
        anc_depth=anc_depth,
        P=P,
        has_bp=has_bp,
    )


def label_matrix(onto: Ontology, selected: Sequence[int], protein_terms: Sequence[Sequence[int]]) -> np.ndarray:
    sel = {t: j for j, t in enumerate(selected)}
    Y = np.zeros((len(protein_terms), len(selected)), dtype=bool)
    for p, ts in enumerate(protein_terms):
        if len(ts) == 0:
            continue
        for a in onto.propagate(ts).tolist():
            j = sel.get(a)
            if j is not None:
                Y[p, j] = True
    return Y


# --------------------------------------------------------------------------- #
# config / paths
# --------------------------------------------------------------------------- #
def load_species_config(path: Optional[Path] = None) -> Dict:
    import yaml

    path = Path(path) if path else REPO_ROOT / "configs" / "species.yaml"
    with open(path) as fh:
        return yaml.safe_load(fh)


def raw_dir(data_dir: Path, source: str) -> Path:
    return Path(data_dir) / "raw" / source


def species_files(species: str, source: str, data_dir: Path, cfg: Optional[Dict] = None) -> Dict[str, Path]:
    cfg = cfg or load_species_config()
    sp_cfg = cfg["species"][species][source]
    src_cfg = cfg["sources"][source]
    base = raw_dir(data_dir, source)
    out = {
        "go": base / src_cfg["go"],
        "gaf": base / sp_cfg["gaf"],
        "biogrid": base / sp_cfg["biogrid"],
    }
    if "biogrid_member" in sp_cfg:
        out["biogrid_member"] = sp_cfg["biogrid_member"]
    if src_cfg.get("expression"):
        out["expression"] = base / src_cfg["expression"]
    return out


# --------------------------------------------------------------------------- #
# main builder
# --------------------------------------------------------------------------- #
class _Cache:
    onto: Dict[str, Ontology] = {}
    diam: Dict[str, float] = {}


def load_bp_ontology(go_path: Path) -> Tuple[Ontology, float]:
    key = str(go_path)
    if key not in _Cache.onto:
        onto = read_obo(go_path, relations=("is_a", "part_of"), namespace="biological_process")
        _Cache.onto[key] = onto
    onto = _Cache.onto[key]
    if key not in _Cache.diam:
        diam_file = Path(str(go_path) + ".bp_diameter.txt")
        if diam_file.exists():
            _Cache.diam[key] = float(diam_file.read_text().strip())
        else:
            d = graph_diameter(onto.undirected_adjacency())
            try:
                diam_file.write_text(str(d))
            except OSError:
                pass
            _Cache.diam[key] = d
    return onto, _Cache.diam[key]


def node_id(mapper: IdMapper, names: Sequence[str]) -> str:
    """Canonical node id of a BioGRID interactor.

    One GAF object -> its canonical id; several (e.g. one symbol with several
    UniProt entries) -> the ids joined by ``|`` (annotations are merged);
    none -> ``BIOGRID:<first identifier>``.
    """
    hits = mapper.map_all(names)
    if not hits:
        return f"BIOGRID:{names[0]}"
    return "|".join(hits)


def node_annotations(node: str, bp: Dict[str, set]) -> List[int]:
    if node.startswith("BIOGRID:"):
        return []
    out: set = set()
    for c in node.split("|"):
        out |= bp.get(c, set())
    return sorted(out)


def read_network_and_annotations(
    species: str,
    source: str = "original",
    data_dir: Path = REPO_ROOT / "data",
    physical_only: bool = True,
    cross_species: bool = False,
):
    """Return (onto, go_diameter, mapper, edges, bp_annotations, interactions).

    Interactions are restricted to the species' NCBI taxa (``taxon`` in
    ``configs/species.yaml``; both partners) unless ``cross_species`` is set,
    in which case one partner may come from another organism if its name maps
    to the species' GAF (as orthologous symbols often do).
    """
    cfg = load_species_config()
    files = species_files(species, source, data_dir, cfg)
    onto, diam = load_bp_ontology(files["go"])
    records = read_gaf(files["gaf"], aspects=None, exclude_evidence=(), exclude_not=False)
    mapper = IdMapper.from_records(records)
    bp: Dict[str, set] = defaultdict(set)
    for r in records:
        if r.aspect != "P" or r.evidence == "IEA" or r.is_not:
            continue
        t = onto.resolve(r.go_id)
        if t is not None:
            bp[r.canonical_id].add(t)

    biogrid_path = files["biogrid"]
    if "biogrid_member" in files and str(biogrid_path).endswith(".zip"):
        import tempfile
        import zipfile

        with zipfile.ZipFile(biogrid_path) as zf:
            member = [n for n in zf.namelist() if files["biogrid_member"] in n][0]
            tmp = Path(tempfile.mkdtemp()) / Path(member).name
            tmp.write_bytes(zf.read(member))
            biogrid_path = tmp
    inter = read_biogrid(biogrid_path, physical_only=physical_only)
    taxa = {str(t) for t in cfg["species"][species].get("taxon", [])}
    if not taxa:
        counts = Counter(it.taxon_a for it in inter)
        taxa = {counts.most_common(1)[0][0]} if counts else set()
    edges = []
    for it in inter:
        ina, inb = it.taxon_a in taxa, it.taxon_b in taxa
        if not (ina and inb):
            if not cross_species or not (ina or inb):
                continue
        a = node_id(mapper, it.names_a)
        b = node_id(mapper, it.names_b)
        if cross_species and not (ina and inb) and (a.startswith("BIOGRID:") or b.startswith("BIOGRID:")):
            continue
        if a != b:
            edges.append((a, b))
    return onto, diam, mapper, edges, bp, inter


def build_dataset(
    species: str,
    n_terms: int = 90,
    source: str = "original",
    data_dir: Path = REPO_ROOT / "data",
    min_count: int = 3,
    max_frac: float = 0.05,
    terms_file: Optional[Path] = None,
    network: str = "ppi",
    cross_species: Optional[bool] = None,
) -> Dataset:
    """Build the PPI dataset of one species (Section 2.5).

    Physical BioGRID interactions, duplicates merged, largest connected
    component, non-IEA BP annotations, ``n_terms`` selected GO terms.
    """
    if cross_species is None:  # species default from configs/species.yaml
        cross_species = bool(load_species_config()["species"][species].get("cross_species", False))
    onto, diam, mapper, edges, bp, inter = read_network_and_annotations(species, source, data_dir,
                                                                        cross_species=cross_species)
    nodes_all, A_all = edges_to_adjacency(edges)
    n_all = len(nodes_all)
    m_all = int(A_all.nnz // 2)
    lcc = largest_component(A_all)
    nodes = [nodes_all[i] for i in lcc]
    A = A_all[lcc][:, lcc].tocsr()

    protein_terms = [node_annotations(p, bp) for p in nodes]
    if terms_file is not None:
        wanted = [l.split()[0] for l in Path(terms_file).read_text().splitlines() if l.strip() and not l.startswith("#")]
        selected = [onto.resolve(t) for t in wanted]
        selected = [t for t in selected if t is not None]
    else:
        counts: Counter = Counter()
        n_ann = 0
        for ts in protein_terms:
            if ts:
                n_ann += 1
                counts.update(onto.propagate(ts).tolist())
        selected = select_terms(onto, counts, n_terms, n_ann, min_count=min_count, max_frac=max_frac)

    Y = label_matrix(onto, selected, protein_terms)
    go = build_go_info(onto, selected, protein_terms, go_diameter=diam)
    meta = {
        "species": species,
        "source": source,
        "network": network,
        "cross_species": cross_species,
        "proteins_total": n_all,
        "interactions_total": m_all,
        "annotated_proteins_total": int(sum(1 for p in nodes_all if node_annotations(p, bp))),
        "mapped_proteins_total": int(sum(1 for p in nodes_all if not p.startswith("BIOGRID:"))),
        "annotated_proteins_lcc": int(go.has_bp.sum()),
        "go_terms_bp": len(onto),
        "go_diameter": diam,
    }
    ds = Dataset(
        name=f"{species}-{network}-{len(selected)}",
        proteins=nodes,
        adj=A,
        terms=[onto.ids[t] for t in selected],
        term_names=[onto.names[t] for t in selected],
        Y=Y,
        go=go,
        meta=meta,
    )
    log.info("built %s: %s", ds.name, ds.summary())
    return ds


def processed_path(data_dir: Path, source: str, species: str, network: str, n_terms: int,
                   options: Optional[Dict] = None) -> Path:
    suffix = ""
    if options:
        import hashlib
        import json

        blob = json.dumps({k: str(v) for k, v in sorted(options.items())})
        suffix = "_" + hashlib.md5(blob.encode()).hexdigest()[:8]
    return Path(data_dir) / "processed" / source / f"{species}_{network}_{n_terms}{suffix}.pkl"


def get_dataset(
    species: str,
    n_terms: int = 90,
    source: str = "original",
    data_dir: Path = REPO_ROOT / "data",
    network: str = "ppi",
    rebuild: bool = False,
    **kwargs,
) -> Dataset:
    """Load a cached dataset, building (and caching) it if necessary."""
    path = processed_path(data_dir, source, species, network, n_terms, kwargs)
    if path.exists() and not rebuild:
        ds = Dataset.load(path)
        if network != "ppi":
            return ds
        want = kwargs.get("cross_species")
        if want is None:
            want = bool(load_species_config()["species"][species].get("cross_species", False))
        if bool(ds.meta.get("cross_species", False)) == bool(want):
            return ds
        log.info("cached %s was built with other network options; rebuilding", path)
    if network == "ppi":
        ds = build_dataset(species, n_terms, source, data_dir, **kwargs)
    elif network == "integrated":
        from .integrated import build_integrated_dataset

        ds = build_integrated_dataset(species, n_terms, source, data_dir, **kwargs)
    else:
        raise ValueError(f"unknown network type {network!r}")
    ds.save(path)
    return ds
