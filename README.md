# Metric Labeling and Semi-metric Embedding for Protein Annotation Prediction — replication

A from-scratch, open-source Python implementation of

> Emre Sefer and Carl Kingsford. *Metric Labeling and Semi-metric Embedding for
> Protein Annotation Prediction.* (Early version: RECOMB 2011.)

It covers every component of the paper: the GO distances, the least-squared-
distortion (LSD) metric approximation, the Metric Labeling LP and its
`O(log k)` rounding, the five baselines plus the semimetric MAP algorithm of
Kumar & Koller, the integrated yeast network, five-fold cross-validation, and
config files that regenerate Figures 1–3, Table 1, Table 2 and the robustness
experiments.

No commercial solver is needed. HiGHS (via SciPy), a numba triangle-fixing
solver for LSD, a PDHG first-order LP solver and PyMaxflow graph cuts cover
everything. Gurobi and CPLEX are picked up as optional backends when they are
installed.

---

## Quick start

```bash
git clone <this repo> && cd metriclabeling-replication
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"            # add ".[gurobi]" / ".[cplex]" for commercial solvers

pytest -q                          # unit tests (~1 min)
mlabel synthetic                   # end-to-end check on a planted GO-like benchmark
mlabel download                    # 2011 data snapshot of the paper (~45 MB)
mlabel run configs/experiments/smoke.yaml --jobs 4   # S. pombe, all methods, ~1 min
```

Everything else:

```bash
mlabel table1                                        # Table 1 (network statistics)
mlabel run configs/experiments/fig1_yeast_ppi.yaml --jobs 8
mlabel run configs/experiments/fig1_yeast_integrated.yaml
mlabel run configs/experiments/fig2_terms.yaml
mlabel run configs/experiments/fig3a_distances.yaml configs/experiments/fig3b_alpha.yaml
mlabel run configs/experiments/table2.yaml --jobs 8  # prints the Table 2 layout
mlabel run configs/experiments/robustness_edges.yaml configs/experiments/robustness_misannotation.yaml
mlabel lsd --species yeast --distances sp,lca,lin --alphas 0 0.1 0.3   # LSD statistics
# or simply
JOBS=8 bash scripts/reproduce_all.sh
```

Runs can be resumed. Each (dataset, fold, method) result is appended to
`results/<experiment>/runs.jsonl`, and finished entries are skipped on
re-runs. `--folds 0,1`, `--methods "LCA 0.3 LSD"`, `--species pombe` and
`--runs 10` restrict a run. Outputs are described in
[`results/README.md`](results/README.md).

---

## Repository layout

```
configs/species.yaml          data sources (2011 snapshot + current releases)
configs/experiments/*.yaml    one file per figure / table / ablation
src/metriclabeling/
  io/                         OBO, GAF, BioGRID readers
  dataset.py                  networks, term selection, GO pre-computations
  integrated.py               integrated yeast network
  distances.py                d_SP, d_LCA, d_Lin, d_KB, Eq. 16
  lsd.py                      LSD metric approximation
  hst.py                      FRT trees, tree selection
  lp.py                       LP relaxations, solvers, rounding, ILP
  graphcut.py                 min-cut, alpha-expansion, hierarchical moves
  methods/                    Metric Labeling + baselines
  folds.py evaluation.py      CV, perturbations, ROC/AUC
  experiment.py plotting.py   runner, figures, Table 2
  cli.py                      `mlabel` command
tests/                        unit + integration tests
```

## Citation

```bibtex
@inproceedings{sefer2011metric,
  title     = {Metric Labeling and Semi-metric Embedding for Protein Annotation Prediction},
  author    = {Sefer, Emre and Kingsford, Carl},
  booktitle = {Research in Computational Molecular Biology (RECOMB)},
  year      = {2011}
}
```

(Fill in the volume/pages, or cite the extended journal version.)
