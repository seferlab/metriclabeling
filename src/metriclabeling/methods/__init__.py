"""Prediction methods.  Every method exposes ``predict(fold) -> (n, k) scores in [0, 1]``."""

from __future__ import annotations

from typing import Dict

from .baselines import MRF, FunctionalFlow, GenMultiCut, Majority, Neighborhood, SemiMAP
from .metric_labeling import MetricLabeling

REGISTRY = {
    "majority": Majority,
    "neighborhood": Neighborhood,
    "neighbour": Neighborhood,
    "functional_flow": FunctionalFlow,
    "functionalflow": FunctionalFlow,
    "mrf": MRF,
    "genmulticut": GenMultiCut,
    "multicut": GenMultiCut,
    "semimap": SemiMAP,
    "metric_labeling": MetricLabeling,
    "ml": MetricLabeling,
}


def make_method(spec: Dict):
    """Instantiate a method from a config dict ``{"method": name, **params}``."""
    spec = dict(spec)
    name = spec.pop("method").lower()
    spec.pop("name", None)
    spec.pop("label", None)
    if name not in REGISTRY:
        raise KeyError(f"unknown method {name!r}; available: {sorted(REGISTRY)}")
    return REGISTRY[name](**spec)


__all__ = ["REGISTRY", "make_method", "Majority", "Neighborhood", "FunctionalFlow", "MRF", "GenMultiCut",
           "SemiMAP", "MetricLabeling"]
