"""Minimal, dependency-free OBO 1.2/1.4 reader for the Gene Ontology.

Only what the paper needs is kept: term ids, names, namespaces, ``is_a`` and
``relationship`` edges (``part_of`` by default), alternative ids and obsolete
flags.  The resulting :class:`Ontology` is a DAG over integer indices with
helpers for ancestors, upward BFS distances and undirected shortest paths.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp

from .files import PathLike, open_text

BP_ROOT = "GO:0008150"
MF_ROOT = "GO:0003674"
CC_ROOT = "GO:0005575"
NAMESPACE_ALIASES = {
    "P": "biological_process",
    "BP": "biological_process",
    "F": "molecular_function",
    "MF": "molecular_function",
    "C": "cellular_component",
    "CC": "cellular_component",
}


@dataclass
class _RawTerm:
    id: str
    name: str = ""
    namespace: str = ""
    is_obsolete: bool = False
    alt_ids: List[str] = field(default_factory=list)
    parents: List[Tuple[str, str]] = field(default_factory=list)  # (parent id, relation)


def _parse_obo(path: PathLike) -> Tuple[List[_RawTerm], str]:
    terms: List[_RawTerm] = []
    default_ns = ""
    cur: Optional[_RawTerm] = None
    in_term = False
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("!"):
                continue
            if line.startswith("["):
                if cur is not None:
                    terms.append(cur)
                cur = None
                in_term = line == "[Term]"
                continue
            if ":" not in line:
                continue
            key, _, val = line.partition(":")
            val = val.strip()
            # strip trailing comments / qualifiers
            if " ! " in val:
                val = val.split(" ! ", 1)[0].strip()
            if not in_term:
                if key == "default-namespace":
                    default_ns = val
                continue
            if key == "id":
                cur = _RawTerm(id=val)
            elif cur is None:
                continue
            elif key == "name":
                cur.name = val
            elif key == "namespace":
                cur.namespace = val
            elif key == "is_obsolete":
                cur.is_obsolete = val.lower().startswith("true")
            elif key == "alt_id":
                cur.alt_ids.append(val)
            elif key == "is_a":
                cur.parents.append((val.split()[0], "is_a"))
            elif key == "relationship":
                parts = val.split()
                if len(parts) >= 2:
                    cur.parents.append((parts[1], parts[0]))
    if cur is not None:
        terms.append(cur)
    return terms, default_ns


class Ontology:
    """A DAG of ontology terms over integer indices.

    ``parents[i]`` lists the indices of the direct parents of term ``i``
    (restricted to the relations that were requested when building it).
    """

    def __init__(
        self,
        ids: Sequence[str],
        names: Sequence[str],
        namespaces: Sequence[str],
        parents: Sequence[Sequence[int]],
        alt_ids: Optional[Dict[str, str]] = None,
    ) -> None:
        self.ids: List[str] = list(ids)
        self.names: List[str] = list(names)
        self.namespaces: List[str] = list(namespaces)
        self.parents: List[List[int]] = [sorted(set(p)) for p in parents]
        self.index: Dict[str, int] = {t: i for i, t in enumerate(self.ids)}
        self.alt_ids: Dict[str, str] = dict(alt_ids or {})
        n = len(self.ids)
        children: List[List[int]] = [[] for _ in range(n)]
        for c, ps in enumerate(self.parents):
            for p in ps:
                children[p].append(c)
        self.children: List[List[int]] = children
        self._anc_cache: Dict[int, np.ndarray] = {}
        self._up_cache: Dict[int, Dict[int, int]] = {}

    # ------------------------------------------------------------------ basics
    def __len__(self) -> int:
        return len(self.ids)

    def resolve(self, go_id: str) -> Optional[int]:
        """Index of ``go_id`` (following ``alt_id`` links), or ``None``."""
        if go_id in self.index:
            return self.index[go_id]
        primary = self.alt_ids.get(go_id)
        if primary is not None:
            return self.index.get(primary)
        return None

    @property
    def roots(self) -> List[int]:
        return [i for i, p in enumerate(self.parents) if not p]

    @property
    def root(self) -> int:
        roots = self.roots
        if len(roots) != 1:
            raise ValueError(f"ontology has {len(roots)} roots; restrict to one namespace first")
        return roots[0]

    # ------------------------------------------------------------- traversals
    def ancestors(self, i: int) -> np.ndarray:
        """Sorted indices of all ancestors of ``i`` *including* ``i``."""
        cached = self._anc_cache.get(i)
        if cached is not None:
            return cached
        seen = {i}
        stack = [i]
        while stack:
            u = stack.pop()
            for p in self.parents[u]:
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        arr = np.fromiter(sorted(seen), dtype=np.int64)
        self._anc_cache[i] = arr
        return arr

    def up_distances(self, i: int) -> Dict[int, int]:
        """Shortest *upward* path length from ``i`` to each of its ancestors."""
        cached = self._up_cache.get(i)
        if cached is not None:
            return cached
        dist = {i: 0}
        dq = deque([i])
        while dq:
            u = dq.popleft()
            for p in self.parents[u]:
                if p not in dist:
                    dist[p] = dist[u] + 1
                    dq.append(p)
        self._up_cache[i] = dist
        return dist

    def propagate(self, terms: Iterable[int]) -> np.ndarray:
        """Union of the ancestor sets (true-path rule)."""
        out: set = set()
        for t in terms:
            out.update(self.ancestors(int(t)).tolist())
        return np.fromiter(sorted(out), dtype=np.int64)

    def depth(self, i: int) -> int:
        """Shortest path length from the root down to ``i``."""
        return self.up_distances(i)[self.root]

    def undirected_adjacency(self) -> sp.csr_matrix:
        rows, cols = [], []
        for c, ps in enumerate(self.parents):
            for p in ps:
                rows += [c, p]
                cols += [p, c]
        n = len(self)
        data = np.ones(len(rows), dtype=np.float64)
        A = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
        A.data[:] = 1.0
        return A

    # ------------------------------------------------------------ restriction
    def restrict(self, namespace: str) -> "Ontology":
        """Sub-ontology for one namespace (edges outside it are dropped)."""
        namespace = NAMESPACE_ALIASES.get(namespace, namespace)
        keep = [i for i, ns in enumerate(self.namespaces) if ns == namespace]
        remap = {old: new for new, old in enumerate(keep)}
        parents = [[remap[p] for p in self.parents[i] if p in remap] for i in keep]
        ids = [self.ids[i] for i in keep]
        alt = {a: p for a, p in self.alt_ids.items() if p in set(ids)}
        return Ontology(
            ids,
            [self.names[i] for i in keep],
            [self.namespaces[i] for i in keep],
            parents,
            alt,
        )


def read_obo(
    path: PathLike,
    relations: Sequence[str] = ("is_a", "part_of"),
    namespace: Optional[str] = None,
) -> Ontology:
    """Read an OBO file.

    Parameters
    ----------
    relations:
        Edge types to keep (``is_a`` is always a relation name here).  The
        paper's GO distances use the DAG structure; we use ``is_a`` and
        ``part_of`` by default, as the original code did.
    namespace:
        Optionally restrict to one namespace (``"P"``/``"biological_process"``).
    """
    raw, default_ns = _parse_obo(path)
    raw = [t for t in raw if not t.is_obsolete]
    ids = [t.id for t in raw]
    index = {t: i for i, t in enumerate(ids)}
    rel = set(relations)
    parents = []
    for t in raw:
        ps = [index[p] for p, r in t.parents if r in rel and p in index]
        parents.append(ps)
    alt: Dict[str, str] = {}
    for t in raw:
        for a in t.alt_ids:
            alt[a] = t.id
    onto = Ontology(
        ids,
        [t.name for t in raw],
        [t.namespace or default_ns for t in raw],
        parents,
        alt,
    )
    if namespace is not None:
        onto = onto.restrict(namespace)
    return onto
