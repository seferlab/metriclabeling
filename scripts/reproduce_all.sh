#!/usr/bin/env bash
# Reproduce every table and figure of the paper.
#   JOBS=8 bash scripts/reproduce_all.sh
# Runs are resumable: re-running skips (dataset, fold, method) results that exist.
set -euo pipefail
JOBS=${JOBS:-4}
cd "$(dirname "$0")/.."

mlabel download --source original
mlabel build --species all --terms 90 150 300
mlabel table1
mlabel synthetic
for cfg in fig1_yeast_ppi fig1_yeast_integrated fig2_terms fig3a_distances fig3b_alpha \
           table2 robustness_edges robustness_misannotation ablation; do
  echo "=== $cfg"
  mlabel run "configs/experiments/${cfg}.yaml" --jobs "$JOBS" -v
done
echo "Figures: results/*/figures/  Table 2: results/table2/table2_auc.md"
