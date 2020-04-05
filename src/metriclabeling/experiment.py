"""YAML-driven experiment runner with resumable, optionally parallel execution.

A config looks like::

    name: fig1_yeast_ppi
    source: original          # data snapshot
    folds: 5
    seed: 0
    datasets:
      - {species: yeast, n_terms: 90, network: ppi}
    perturbation: {type: edge_removal, levels: [0.0, 0.1, 0.2]}   # optional
    methods:
      - {label: Majority, method: majority}
      - {label: "LCA 0.3 LSD", method: metric_labeling, distance: lca, alpha: 0.3}

Every (dataset, perturbation level, fold, method) result is appended to
``<out>/<name>/runs.jsonl`` so interrupted runs resume where they stopped.
"""

from __future__ import annotations

import json
import logging
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np

from .dataset import REPO_ROOT, get_dataset
from .evaluation import aggregate, evaluate_fold
from .folds import make_folds, misannotate, remove_edges
from .methods import make_method

log = logging.getLogger(__name__)


def load_config(path: Path) -> Dict:
    import yaml

    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg.setdefault("name", Path(path).stem)
    cfg.setdefault("source", "original")
    cfg.setdefault("folds", 5)
    cfg.setdefault("seed", 0)
    cfg.setdefault("perturbation", None)
    return cfg


def _dataset_key(d: Dict) -> str:
    return f"{d['species']}-{d.get('network', 'ppi')}-{d.get('n_terms', 90)}"


def _levels(cfg: Dict) -> List[float]:
    p = cfg.get("perturbation")
    return [0.0] if not p else [float(x) for x in p["levels"]]


def _jsonable(v):
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in v.items()}
    return v


def run_fold_task(cfg: Dict, dspec: Dict, level: float, fold_id: int, methods: List[Dict],
                  data_dir: str) -> List[Dict]:
    """Run ``methods`` on one fold; returns one record per method."""
    import warnings

    warnings.filterwarnings("ignore")
    ds_kwargs = {k: v for k, v in dspec.items() if k not in ("species", "n_terms", "network")}
    ds = get_dataset(dspec["species"], dspec.get("n_terms", 90), cfg["source"], Path(data_dir),
                     network=dspec.get("network", "ppi"), **ds_kwargs)
    pert = (cfg.get("perturbation") or {}).get("type")
    adj = ds.adj
    if pert == "edge_removal" and level > 0:
        adj = remove_edges(ds.adj, level, seed=cfg["seed"] + 17)
    folds = make_folds(ds, cfg["folds"], seed=cfg["seed"], adj=adj)
    fold = folds[fold_id]
    if pert == "misannotation" and level > 0:
        fold = misannotate(fold, level, seed=cfg["seed"])
    out = []
    for spec in methods:
        label = spec.get("label") or spec.get("name") or spec["method"]
        t0 = time.perf_counter()
        rec = {"dataset": _dataset_key(dspec), "species": dspec["species"], "n_terms": dspec.get("n_terms", 90),
               "network": dspec.get("network", "ppi"), "level": level, "fold": fold_id, "method": label}
        try:
            m = make_method({**spec, "seed": spec.get("seed", cfg["seed"])})
            S = m.predict(fold)
            ev = evaluate_fold(ds.Y, S, fold.test)
            rec.update(ev)
            rec["seconds"] = time.perf_counter() - t0
            rec["info"] = {k: _jsonable(v) for k, v in getattr(m, "info", {}).items()
                           if isinstance(v, (int, float, str, dict, np.floating, np.integer))}
        except Exception as exc:  # keep going; record the failure
            rec["error"] = f"{type(exc).__name__}: {exc}"
            rec["traceback"] = traceback.format_exc()
            log.error("%s fold %d %s failed: %s", rec["dataset"], fold_id, label, exc)
        out.append(rec)
        log.info("%s level=%.2f fold=%d %-28s auc=%s (%.1fs)", rec["dataset"], level, fold_id, label,
                 f"{rec.get('auc', float('nan')):.3f}", time.perf_counter() - t0)
    return out


def _done(runs_file: Path) -> set:
    done = set()
    if runs_file.exists():
        for line in runs_file.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "error" not in r:
                done.add((r["dataset"], float(r["level"]), int(r["fold"]), r["method"]))
    return done


def run_experiment(cfg: Dict, out_dir: Path = REPO_ROOT / "results", data_dir: Path = REPO_ROOT / "data",
                   jobs: int = 1, folds: Optional[Iterable[int]] = None, methods: Optional[List[str]] = None) -> Path:
    out = Path(out_dir) / cfg["name"]
    out.mkdir(parents=True, exist_ok=True)
    (out / "config.json").write_text(json.dumps(cfg, indent=2))
    runs_file = out / "runs.jsonl"
    done = _done(runs_file)
    fold_ids = list(folds) if folds is not None else list(range(cfg["folds"]))
    mspecs = [m for m in cfg["methods"] if methods is None or (m.get("label") or m["method"]) in methods]

    # make sure datasets exist before forking workers
    for d in cfg["datasets"]:
        kw = {k: v for k, v in d.items() if k not in ("species", "n_terms", "network")}
        get_dataset(d["species"], d.get("n_terms", 90), cfg["source"], Path(data_dir), network=d.get("network", "ppi"), **kw)

    tasks = []
    for d in cfg["datasets"]:
        for level in _levels(cfg):
            for f in fold_ids:
                todo = [m for m in mspecs if (_dataset_key(d), level, f, m.get("label") or m["method"]) not in done]
                if todo:
                    tasks.append((cfg, d, level, f, todo, str(data_dir)))
    log.info("%s: %d fold tasks to run", cfg["name"], len(tasks))

    def _write(records):
        with open(runs_file, "a") as fh:
            for r in records:
                fh.write(json.dumps(r) + "\n")

    if jobs <= 1:
        for t in tasks:
            _write(run_fold_task(*t))
    else:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = [ex.submit(run_fold_task, *t) for t in tasks]
            for fu in as_completed(futs):
                _write(fu.result())
    summarize(out)
    return out


def load_runs(out: Path):
    import pandas as pd

    rows = [json.loads(l) for l in (Path(out) / "runs.jsonl").read_text().splitlines() if l.strip()]
    return pd.DataFrame(rows)


def summarize(out: Path) -> None:
    """Write ``summary.csv`` (mean AUCs) and ``curves.json`` (fold-averaged ROC points)."""
    import pandas as pd

    df = load_runs(out)
    if "error" in df.columns:
        df = df[df["error"].isna()]
    if df.empty:
        return
    # keep the latest record per key
    df = df.drop_duplicates(subset=["dataset", "level", "fold", "method"], keep="last")
    curves = {}
    rows = []
    for (dset, level, method), g in df.groupby(["dataset", "level", "method"], sort=False):
        agg = aggregate(g.to_dict("records"))
        curves.setdefault(dset, {}).setdefault(str(level), {})[method] = {"fpr": agg["fpr"], "tpr": agg["tpr"]}
        rows.append({"dataset": dset, "species": g["species"].iloc[0], "n_terms": g["n_terms"].iloc[0],
                     "network": g["network"].iloc[0], "level": level, "method": method,
                     "auc": agg["auc"], "auc_std": agg["auc_std"], "auc_grid": agg["auc_grid"],
                     "folds": len(g), "seconds": g["seconds"].mean()})
    pd.DataFrame(rows).to_csv(Path(out) / "summary.csv", index=False)
    (Path(out) / "curves.json").write_text(json.dumps(curves))
