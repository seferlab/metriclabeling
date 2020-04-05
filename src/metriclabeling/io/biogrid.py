"""Reader for BioGRID interaction files (``tab`` v1, ``tab2`` and ``tab3``).

The 2011 data used by the paper (BioGRID 3.0.65, ``*.tab.txt``) has no
"experimental system type" column, so physical/genetic is decided from the
experimental-system name using :data:`PHYSICAL_SYSTEMS`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, List, Optional, Sequence, Tuple

from .files import PathLike, open_text

PHYSICAL_SYSTEMS = frozenset(
    {
        "affinity capture-luminescence",
        "affinity capture-ms",
        "affinity capture-rna",
        "affinity capture-western",
        "biochemical activity",
        "co-crystal structure",
        "co-fractionation",
        "co-localization",
        "co-purification",
        "far western",
        "fret",
        "pca",
        "protein-peptide",
        "protein-rna",
        "proximity label-ms",
        "reconstituted complex",
        "two-hybrid",
        "cross-linking-ms (xl-ms)",
    }
)

GENETIC_SYSTEMS = frozenset(
    {
        "dosage growth defect",
        "dosage lethality",
        "dosage rescue",
        "negative genetic",
        "phenotypic enhancement",
        "phenotypic suppression",
        "positive genetic",
        "synthetic growth defect",
        "synthetic haploinsufficiency",
        "synthetic lethality",
        "synthetic rescue",
    }
)


@dataclass
class Interaction:
    names_a: Tuple[str, ...]  # identifiers in priority order (systematic, symbol, aliases, entrez)
    names_b: Tuple[str, ...]
    system: str
    system_type: str  # "physical" or "genetic"
    pubmed: str
    taxon_a: str
    taxon_b: str


def _split_aliases(s: str) -> List[str]:
    if not s or s in ("-", "N/A"):
        return []
    return [a for a in s.split("|") if a and a not in ("-", "N/A")]


def _system_type(system: str, explicit: Optional[str]) -> str:
    if explicit:
        return explicit.strip().lower()
    s = system.strip().lower()
    if s in PHYSICAL_SYSTEMS:
        return "physical"
    if s in GENETIC_SYSTEMS:
        return "genetic"
    return "unknown"


def iter_biogrid(path: PathLike) -> Iterator[Interaction]:
    """Yield interactions from any BioGRID tab format."""
    with open_text(path) as fh:
        header: Optional[List[str]] = None
        for line in fh:
            line = line.rstrip("\n")
            if header is None:
                if line.startswith("INTERACTOR_A") or line.startswith("#BioGRID Interaction ID"):
                    header = [h.strip().lstrip("#") for h in line.split("\t")]
                    col = {h.lower(): i for i, h in enumerate(header)}
                    v1 = line.startswith("INTERACTOR_A")
                continue
            if not line.strip():
                continue
            c = line.split("\t")
            if v1:
                if len(c) < 11:
                    continue
                names_a = tuple([c[0], c[2]] + _split_aliases(c[4]))
                names_b = tuple([c[1], c[3]] + _split_aliases(c[5]))
                system, pubmed = c[6], c[8]
                ta, tb = c[9], c[10]
                stype = _system_type(system, None)
            else:
                def g(name: str, default: str = "") -> str:
                    i = col.get(name.lower())
                    return c[i] if i is not None and i < len(c) else default

                names_a = tuple(
                    x
                    for x in [g("Systematic Name Interactor A"), g("Official Symbol Interactor A")]
                    + _split_aliases(g("Synonyms Interactor A"))
                    + [g("Entrez Gene Interactor A")]
                    if x and x != "-"
                )
                names_b = tuple(
                    x
                    for x in [g("Systematic Name Interactor B"), g("Official Symbol Interactor B")]
                    + _split_aliases(g("Synonyms Interactor B"))
                    + [g("Entrez Gene Interactor B")]
                    if x and x != "-"
                )
                system = g("Experimental System")
                stype = _system_type(system, g("Experimental System Type") or None)
                pubmed = g("Pubmed ID") or g("Publication Source")
                ta = g("Organism Interactor A") or g("Organism ID Interactor A")
                tb = g("Organism Interactor B") or g("Organism ID Interactor B")
            yield Interaction(names_a, names_b, system, stype, pubmed, ta.strip(), tb.strip())


def read_biogrid(
    path: PathLike,
    physical_only: bool = True,
    taxon: Optional[str] = None,
    systems: Optional[Sequence[str]] = None,
    pubmeds: Optional[Sequence[str]] = None,
) -> List[Interaction]:
    """Read BioGRID interactions with optional filters.

    Parameters
    ----------
    physical_only: keep only physical interactions (as in the paper).
    taxon: keep interactions where *both* partners have this NCBI taxon id.
    systems: keep only these experimental systems (case-insensitive).
    pubmeds: keep only interactions from these PubMed ids.
    """
    sysset = {s.lower() for s in systems} if systems else None
    pmset = {str(p).replace("PUBMED:", "") for p in pubmeds} if pubmeds else None
    out = []
    for it in iter_biogrid(path):
        if physical_only and it.system_type != "physical":
            continue
        if taxon is not None and (it.taxon_a != str(taxon) or it.taxon_b != str(taxon)):
            continue
        if sysset is not None and it.system.lower() not in sysset:
            continue
        if pmset is not None:
            ids = {p.strip().replace("PUBMED:", "") for p in it.pubmed.replace(";", "|").split("|")}
            if not ids & pmset:
                continue
        out.append(it)
    return out
