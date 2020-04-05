"""Figures (ROC curves as in Figs. 1-3, robustness curves) and Table 2."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

# validated categorical order (blue, orange, aqua, yellow, magenta, green, violet, red)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
BASELINE_GRAYS = ["#3d3c39", "#6b6a66", "#8f8e89", "#55544f", "#a3a29c", "#77766f"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h"]
BASELINES = {"majority", "neighbour", "neighborhood", "genmulticut", "functionalflow", "functional_flow", "mrf"}
TEXT = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e4e3df"


def _style(ax, xlabel: str, ylabel: str, title: Optional[str] = None):
    ax.set_xlabel(xlabel, color=TEXT)
    ax.set_ylabel(ylabel, color=TEXT)
    if title:
        ax.set_title(title, color=TEXT, fontsize=11, loc="left")
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED)


def _series_styles(labels: Sequence[str]) -> Dict[str, Dict]:
    """Baselines: recessive grays + distinct markers; Metric Labeling variants: categorical hues."""
    styles, ci, bi = {}, 0, 0
    for i, lab in enumerate(labels):
        key = lab.lower().replace(" ", "")
        if key in BASELINES:
            styles[lab] = dict(color=BASELINE_GRAYS[bi % len(BASELINE_GRAYS)], linestyle="--",
                               marker=MARKERS[i % len(MARKERS)], linewidth=1.5)
            bi += 1
        else:
            styles[lab] = dict(color=SERIES[ci % len(SERIES)], linestyle="-",
                               marker=MARKERS[i % len(MARKERS)], linewidth=2.0)
            ci += 1
    return styles


def plot_roc(curves: Dict[str, Dict], path: Path, title: Optional[str] = None, order: Optional[List[str]] = None):
    """``curves``: label -> {"fpr": [...], "tpr": [...]} ."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = order or list(curves)
    labels = [l for l in labels if l in curves]
    styles = _series_styles(labels)
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    ax.plot([0, 1], [0, 1], color=GRID, linewidth=1, zorder=0)
    for lab in labels:
        c = curves[lab]
        x = np.r_[0.0, c["fpr"], 1.0]
        y = np.r_[0.0, c["tpr"], 1.0]
        o = np.lexsort((y, x))
        ax.plot(x[o], y[o], label=lab, markersize=4.5, markeredgecolor="white", markeredgewidth=0.6, **styles[lab])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.01)
    _style(ax, "False positive rate", "Sensitivity", title)
    leg = ax.legend(loc="lower right", fontsize=8, frameon=False)
    for t in leg.get_texts():
        t.set_color(TEXT)
    fig.tight_layout()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_experiment(out_dir: Path, combine: bool = False, fig_dir: Optional[Path] = None) -> List[Path]:
    """ROC figure(s) for an experiment directory (one per dataset/level, or all
    datasets in one panel when ``combine`` - used for Fig. 2)."""
    out_dir = Path(out_dir)
    curves = json.loads((out_dir / "curves.json").read_text())
    cfg = json.loads((out_dir / "config.json").read_text())
    order = [m.get("label") or m["method"] for m in cfg["methods"]]
    fig_dir = Path(fig_dir) if fig_dir else out_dir / "figures"
    paths = []
    if combine:
        merged, morder = {}, []
        for dset, by_level in curves.items():
            terms = dset.split("-")[-1]
            for lab in order:
                if lab in by_level.get("0.0", {}):
                    name = f"{lab} {terms}"
                    merged[name] = by_level["0.0"][lab]
                    morder.append(name)
        p = fig_dir / f"{cfg['name']}.png"
        plot_roc(merged, p, cfg.get("title"), morder)
        return [p]
    for dset, by_level in curves.items():
        for level, cur in by_level.items():
            suffix = "" if float(level) == 0 else f"_level{float(level):.2f}"
            p = fig_dir / f"{cfg['name']}_{dset}{suffix}.png"
            title = cfg.get("title") or dset
            if float(level) > 0:
                title += f" ({cfg['perturbation']['type'].replace('_', ' ')} {float(level):.0%})"
            plot_roc(cur, p, title, order)
            paths.append(p)
    return paths


def plot_robustness(out_dir: Path, fig_dir: Optional[Path] = None) -> List[Path]:
    """Mean AUC vs perturbation level, one line per method."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    out_dir = Path(out_dir)
    df = pd.read_csv(out_dir / "summary.csv")
    cfg = json.loads((out_dir / "config.json").read_text())
    order = [m.get("label") or m["method"] for m in cfg["methods"]]
    fig_dir = Path(fig_dir) if fig_dir else out_dir / "figures"
    paths = []
    for dset, g in df.groupby("dataset"):
        labels = [l for l in order if l in set(g["method"])]
        styles = _series_styles(labels)
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
        for lab in labels:
            h = g[g["method"] == lab].sort_values("level")
            ax.plot(h["level"] * 100, h["auc"], label=lab, markersize=5, markeredgecolor="white", **styles[lab])
        kind = cfg["perturbation"]["type"].replace("_", " ")
        _style(ax, f"{kind} (%)", "AUC", f"{dset}: robustness to {kind}")
        leg = ax.legend(fontsize=8, frameon=False, loc="best")
        for t in leg.get_texts():
            t.set_color(TEXT)
        fig.tight_layout()
        p = fig_dir / f"{cfg['name']}_{dset}_auc.png"
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=200)
        fig.savefig(p.with_suffix(".pdf"))
        plt.close(fig)
        paths.append(p)
    return paths


SPECIES_ORDER = ["arabidopsis", "fly", "human", "pombe", "mouse", "worm", "yeast"]
SPECIES_ABBR = {"arabidopsis": "A. tha.", "fly": "D. mel.", "human": "H. sap.", "pombe": "S. pom.",
                "mouse": "M. mus.", "worm": "C. ele.", "yeast": "S. cer."}


def table2(out_dir: Path, metric: str = "auc") -> "object":
    """Table 2 layout: one row per (species, #terms), one column per method."""
    import pandas as pd

    out_dir = Path(out_dir)
    df = pd.read_csv(out_dir / "summary.csv")
    df = df[df["level"] == 0]
    cfg = json.loads((out_dir / "config.json").read_text())
    order = [m.get("label") or m["method"] for m in cfg["methods"]]
    tab = df.pivot_table(index=["species", "n_terms"], columns="method", values=metric)
    tab = tab[[c for c in order if c in tab.columns]]
    sp_rank = {s: i for i, s in enumerate(SPECIES_ORDER)}
    tab = tab.reset_index()
    tab["_o"] = tab["species"].map(sp_rank)
    tab = tab.sort_values(["_o", "n_terms"]).drop(columns="_o")
    tab.insert(0, "Org.", tab.pop("species").map(lambda s: SPECIES_ABBR.get(s, s)))
    tab = tab.rename(columns={"n_terms": "Terms"})
    tab.to_csv(out_dir / f"table2_{metric}.csv", index=False, float_format="%.3f")
    with open(out_dir / f"table2_{metric}.md", "w") as fh:
        fh.write(_markdown(tab))
    return tab


def _markdown(df) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append(f"{v:.3f}" if isinstance(v, float) else str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"
