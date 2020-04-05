"""Readers for the raw data used in the paper (GO OBO, GAF, BioGRID, expression)."""

from .files import open_text
from .obo import Ontology, read_obo
from .gaf import GafRecord, IdMapper, read_gaf
from .biogrid import read_biogrid, PHYSICAL_SYSTEMS, GENETIC_SYSTEMS

__all__ = [
    "open_text",
    "Ontology",
    "read_obo",
    "GafRecord",
    "IdMapper",
    "read_gaf",
    "read_biogrid",
    "PHYSICAL_SYSTEMS",
    "GENETIC_SYSTEMS",
]
