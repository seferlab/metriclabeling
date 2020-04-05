"""Replication of "Metric Labeling and Semi-metric Embedding for Protein
Annotation Prediction" (Emre Sefer and Carl Kingsford).

Main entry points
-----------------
* :func:`metriclabeling.dataset.get_dataset`   - networks + GO labels
* :mod:`metriclabeling.distances`              - d_SP, d_LCA, d_Lin, d_KB (Sec. 2.4)
* :func:`metriclabeling.lsd.lsd`               - least-squared-distortion metric (Sec. 2.3)
* :class:`metriclabeling.methods.MetricLabeling` - the method (Secs. 2.1-2.2)
* :mod:`metriclabeling.methods.baselines`      - comparison methods (Sec. 2.6)
* :func:`metriclabeling.experiment.run_experiment` - 5-fold CV experiments (Sec. 2.7)
"""

__version__ = "0.1.0"
