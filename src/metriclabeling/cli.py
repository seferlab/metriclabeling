"""Command-line interface: ``mlabel <command> ...``."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from .dataset import REPO_ROOT

SPECIES = ["yeast", "pombe", "fly", "worm", "arabidopsis", "mouse", "human"]


def _common(p):
    p.add_argument("--data-dir", default=str(REPO_ROOT / "data"))
    p.add_argument("--source", default="original", choices=["original", "current"])
    p.add_argument("-v", "--verbose", action="store_true")


def cmd_download(a):
    from .download import download

    path = download(a.source, Path(a.data_dir), force=a.force)
    print(f"data in {path}")


def cmd_build(a):
    from .dataset import get_dataset

    species = SPECIES if a.species == "all" else a.species.split(",")
    for s in species:
        for t in a.terms:
            ds = get_dataset(s, t, a.source, Path(a.data_dir), network=a.network, rebuild=a.rebuild)
            print(json.dumps(ds.summary()))


def cmd_table1(a):
    import pandas as pd

    from .dataset import get_dataset

    rows = []
    for s in (SPECIES if a.species == "all" else a.species.split(",")):
        ds = get_dataset(s, 90, a.source, Path(a.data_dir))
        m = ds.meta
        rows.append({"species": s, "proteins": m["proteins_total"], "proteins_in_lcc": ds.n,
                     "annotated_proteins": m["annotated_proteins_total"], "interactions": m["interactions_total"],
                     "interactions_in_lcc": int(ds.adj.nnz // 2)})
    df = pd.DataFrame(rows)
    out = Path(a.out) / "table1.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(df.to_string(index=False))


def cmd_run(a):
    from .experiment import load_config, run_experiment
    from .plotting import plot_experiment, plot_robustness

    for cfg_path in a.config:
        cfg = load_config(Path(cfg_path))
        if a.species:
            cfg["datasets"] = [d for d in cfg["datasets"] if d["species"] in a.species.split(",")]
        if a.runs is not None:
            for m in cfg["methods"]:
                if "runs" in m or m["method"] in ("genmulticut", "semimap"):
                    m["runs"] = a.runs
        folds = [int(x) for x in a.folds.split(",")] if a.folds else None
        methods = a.methods.split(",") if a.methods else None
        out = run_experiment(cfg, Path(a.out), Path(a.data_dir), jobs=a.jobs, folds=folds, methods=methods)
        if not a.no_plot:
            if cfg.get("perturbation"):
                plot_robustness(out)
            plot_experiment(out, combine=cfg.get("plot") == "combined")
            if cfg.get("table2"):
                from .plotting import table2

                print(table2(out).to_string(index=False))
        print(f"results in {out}")


def cmd_plot(a):
    from .plotting import plot_experiment, plot_robustness, table2

    for d in a.result_dirs:
        cfg = json.loads((Path(d) / "config.json").read_text())
        if cfg.get("perturbation"):
            plot_robustness(Path(d))
        paths = plot_experiment(Path(d), combine=cfg.get("plot") == "combined")
        if cfg.get("table2"):
            print(table2(Path(d)).to_string(index=False))
        print("\n".join(map(str, paths)))


def cmd_table2(a):
    from .plotting import table2

    print(table2(Path(a.result_dir), a.metric).to_string(index=False))


def cmd_lsd(a):
    """LSD statistics for the paper's distances on one dataset (fold 0)."""
    import numpy as np

    from .dataset import get_dataset
    from .folds import make_folds

    ds = get_dataset(a.species, a.terms, a.source, Path(a.data_dir))
    fold = make_folds(ds, 5, seed=0)[0]
    print(f"{ds.name}: k={ds.k}")
    for name in a.distances.split(","):
        for alpha in a.alphas:
            S = fold.semimetric(name, alpha)
            M, res = fold.metric(name, alpha, lsd_method=a.method)
            iu = np.triu_indices(ds.k, 1)
            rel = np.sqrt(res.objective / max(np.sum(S[iu] ** 2), 1e-12))
            print(f"  {name:4s} alpha={alpha:.1f} violated triangles={res.violations_before:7d} "
                  f"LSD objective={res.objective:.4f} rel.change={rel:.4f} "
                  f"max violation after={res.max_violation_after:.1e} time={res.seconds:.2f}s")


def cmd_synthetic(a):
    from .synthetic import run_synthetic_benchmark

    run_synthetic_benchmark(Path(a.out), seed=a.seed, n=a.n, k=a.k)


def main(argv=None):
    p = argparse.ArgumentParser(prog="mlabel", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("download", help="download raw data")
    _common(s)
    s.add_argument("--force", action="store_true")
    s.set_defaults(func=cmd_download)

    s = sub.add_parser("build", help="build/cache datasets")
    _common(s)
    s.add_argument("--species", default="all")
    s.add_argument("--terms", type=int, nargs="+", default=[90])
    s.add_argument("--network", default="ppi", choices=["ppi", "integrated"])
    s.add_argument("--rebuild", action="store_true")
    s.set_defaults(func=cmd_build)

    s = sub.add_parser("table1", help="network statistics (Table 1)")
    _common(s)
    s.add_argument("--species", default="all")
    s.add_argument("--out", default=str(REPO_ROOT / "results"))
    s.set_defaults(func=cmd_table1)

    s = sub.add_parser("run", help="run experiment config(s)")
    _common(s)
    s.add_argument("config", nargs="+")
    s.add_argument("--out", default=str(REPO_ROOT / "results"))
    s.add_argument("--jobs", type=int, default=1)
    s.add_argument("--folds", default=None, help="comma-separated fold ids (default: all)")
    s.add_argument("--methods", default=None, help="comma-separated method labels to run")
    s.add_argument("--species", default=None, help="restrict to these species")
    s.add_argument("--runs", type=int, default=None, help="override #runs of GenMultiCut/SemiMAP")
    s.add_argument("--no-plot", action="store_true")
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("plot", help="re-draw figures of result directories")
    s.add_argument("result_dirs", nargs="+")
    s.add_argument("-v", "--verbose", action="store_true")
    s.set_defaults(func=cmd_plot)

    s = sub.add_parser("table2", help="Table 2 from a result directory")
    s.add_argument("result_dir")
    s.add_argument("--metric", default="auc", choices=["auc", "auc_grid"])
    s.add_argument("-v", "--verbose", action="store_true")
    s.set_defaults(func=cmd_table2)

    s = sub.add_parser("lsd", help="LSD statistics for the GO distances")
    _common(s)
    s.add_argument("--species", default="yeast")
    s.add_argument("--terms", type=int, default=90)
    s.add_argument("--distances", default="sp,lca,lin")
    s.add_argument("--alphas", type=float, nargs="+", default=[0.0, 0.1, 0.3])
    s.add_argument("--method", default="triangle_fixing")
    s.set_defaults(func=cmd_lsd)

    s = sub.add_parser("synthetic", help="end-to-end check on a synthetic GO-like benchmark")
    s.add_argument("--out", default=str(REPO_ROOT / "results" / "synthetic"))
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--n", type=int, default=400)
    s.add_argument("--k", type=int, default=20)
    s.add_argument("-v", "--verbose", action="store_true")
    s.set_defaults(func=cmd_synthetic)

    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO if a.verbose else logging.WARNING,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    a.func(a)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
