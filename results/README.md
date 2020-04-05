# Results

Each experiment writes `results/<name>/`:

| file | content |
|---|---|
| `config.json` | the configuration that was run |
| `runs.jsonl` | one JSON record per (dataset, perturbation level, fold, method): AUCs, ROC points, run time, solver info |
| `summary.csv` | fold-averaged exact AUC (`auc`), its std, and the AUC of the 21-threshold ROC polyline (`auc_grid`) |
| `curves.json` | fold-averaged ROC points at score thresholds 0, 0.05, ..., 1 |
| `figures/` | ROC figures (PNG + PDF); AUC-vs-noise plots for robustness runs |
| `table2_auc.{csv,md}` | Table 2 layout (only for `table2`) |

`results/reference/` holds the outputs produced while building this repository
(see the README for what was run).
