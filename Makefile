JOBS ?= 4
PY ?= python3

.PHONY: install test download build smoke synthetic paper fig1 fig2 fig3 table2 robustness ablation

install:
	$(PY) -m pip install -e ".[dev]"

test:
	$(PY) -m pytest -q

download:
	mlabel download --source original

build: download
	mlabel build --species all --terms 90 150 300
	mlabel table1

synthetic:
	mlabel synthetic

smoke: download
	mlabel run configs/experiments/smoke.yaml --jobs $(JOBS)

fig1: download
	mlabel run configs/experiments/fig1_yeast_ppi.yaml configs/experiments/fig1_yeast_integrated.yaml --jobs $(JOBS)

fig2: download
	mlabel run configs/experiments/fig2_terms.yaml --jobs $(JOBS)

fig3: download
	mlabel run configs/experiments/fig3a_distances.yaml configs/experiments/fig3b_alpha.yaml --jobs $(JOBS)

table2: download
	mlabel run configs/experiments/table2.yaml --jobs $(JOBS)

robustness: download
	mlabel run configs/experiments/robustness_edges.yaml configs/experiments/robustness_misannotation.yaml --jobs $(JOBS)

ablation: download
	mlabel run configs/experiments/ablation.yaml --jobs $(JOBS)

paper: build fig1 fig2 fig3 table2 robustness
