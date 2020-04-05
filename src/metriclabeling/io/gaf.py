"""Reader for GO Annotation Files (GAF 1.0 / 2.x) and an identifier mapper.

The paper uses non-IEA annotations of the Biological Process aspect.  BioGRID
uses systematic names / official symbols while GAF files use database object
ids, so :class:`IdMapper` maps every identifier, symbol and synonym that a GAF
file lists for a gene product to one canonical id.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .files import PathLike, open_text


@dataclass(frozen=True)
class GafRecord:
    db: str
    object_id: str
    symbol: str
    qualifier: str
    go_id: str
    reference: str
    evidence: str
    aspect: str
    synonyms: Tuple[str, ...]
    taxon: str

    @property
    def canonical_id(self) -> str:
        return f"{self.db}:{self.object_id}"

    @property
    def is_not(self) -> bool:
        return "NOT" in self.qualifier.upper().split("|")


def read_gaf(
    path: PathLike,
    aspects: Optional[Sequence[str]] = ("P",),
    exclude_evidence: Iterable[str] = ("IEA",),
    exclude_not: bool = True,
    taxon: Optional[str] = None,
    references: Optional[Iterable[str]] = None,
) -> List[GafRecord]:
    """Parse a GAF file.

    Parameters
    ----------
    aspects: keep only these aspects (``P``/``F``/``C``); ``None`` keeps all.
    exclude_evidence: evidence codes to drop (the paper drops ``IEA``).
    exclude_not: drop annotations with a ``NOT`` qualifier.
    taxon: keep only rows whose taxon column contains ``taxon:<id>``.
    references: if given, keep only rows citing one of these references
        (e.g. ``{"PMID:14562095"}`` to extract Huh et al. 2003 localisations).
    """
    aspects_set = set(aspects) if aspects is not None else None
    excl = {e.upper() for e in exclude_evidence}
    refs = set(references) if references is not None else None
    out: List[GafRecord] = []
    with open_text(path) as fh:
        for line in fh:
            if not line or line[0] == "!":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                continue
            aspect = cols[8].strip()
            if aspects_set is not None and aspect not in aspects_set:
                continue
            evidence = cols[6].strip().upper()
            if evidence in excl:
                continue
            rec = GafRecord(
                db=cols[0].strip(),
                object_id=cols[1].strip(),
                symbol=cols[2].strip(),
                qualifier=cols[3].strip(),
                go_id=cols[4].strip(),
                reference=cols[5].strip(),
                evidence=evidence,
                aspect=aspect,
                synonyms=tuple(s.strip() for s in cols[10].split("|") if s.strip()),
                taxon=cols[12].strip(),
            )
            if exclude_not and rec.is_not:
                continue
            if taxon is not None and f"taxon:{taxon}" not in rec.taxon.split("|")[0]:
                continue
            if refs is not None and not (set(rec.reference.split("|")) & refs):
                continue
            out.append(rec)
    return out


class IdMapper:
    """Map gene identifiers / symbols / synonyms to canonical GAF ids.

    Priority: object id > symbol > synonym.  A name that maps to several
    canonical ids at its best priority level is considered ambiguous and is
    not mapped.  Matching is case-insensitive.
    """

    def __init__(self) -> None:
        self._levels: List[Dict[str, Set[str]]] = [defaultdict(set) for _ in range(3)]

    @classmethod
    def from_records(cls, records: Iterable[GafRecord]) -> "IdMapper":
        m = cls()
        for r in records:
            m.add(r.canonical_id, r.object_id, r.symbol, r.synonyms)
        return m

    def add(self, canonical: str, object_id: str, symbol: str, synonyms: Iterable[str]) -> None:
        self._levels[0][object_id.upper()].add(canonical)
        self._levels[0][canonical.upper()].add(canonical)
        if symbol:
            self._levels[1][symbol.upper()].add(canonical)
        for s in synonyms:
            self._levels[2][s.upper()].add(canonical)

    def map(self, name: str) -> Optional[str]:
        key = name.strip().upper()
        if not key or key in ("N/A", "-"):
            return None
        for level in self._levels:
            hits = level.get(key)
            if hits:
                if len(hits) == 1:
                    return next(iter(hits))
                return None  # ambiguous at the best level
        return None

    def map_all(self, names: Iterable[str], max_hits: int = 8) -> List[str]:
        """All canonical ids matched by the first name that matches anything.

        Unlike :meth:`map`, several hits are accepted at the id and symbol
        levels (e.g. one gene symbol with several UniProt entries in GOA
        files); synonyms must still be unambiguous.
        """
        for n in names:
            key = n.strip().upper()
            if not key or key in ("N/A", "-"):
                continue
            for lvl, level in enumerate(self._levels):
                hits = level.get(key)
                if hits:
                    if lvl == 2 and len(hits) > 1:
                        break
                    if len(hits) > max_hits:
                        break
                    return sorted(hits)
        return []

    def map_any(self, names: Iterable[str]) -> Optional[str]:
        """Map the first name (in order) that resolves unambiguously."""
        for n in names:
            hit = self.map(n)
            if hit is not None:
                return hit
        return None
