#!/usr/bin/env python3
"""
STEP 4 -- Figures for the slides.

Produces five PNGs in out/figures/. Each one answers a question a judge will
actually ask, and nothing is plotted just because it is plottable.

  fig1_leakage.png        Why grouped CV. THE methods slide.
  fig2_label.png          What "dissemination" is made of.
  fig3_performance.png    Predicted vs observed, plus the baseline comparison.
  fig4_importance.png     Which backbone features drive the prediction.
  fig5_risk_space.png     The two-component risk space -- the product slide.

Usage
-----
    python figures.py                       # uses out/
    python figures.py --out sweep_full_all  # any run directory
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Validated categorical palette (light surface). Three slots only: the
# all-pairs CVD and normal-vision floors fail past three, so anything needing
# more categories gets faceted or folded into "Other" instead of a fourth hue.
C_BLUE, C_ORANGE, C_AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.8,
    "figure.dpi": 140,
})


def _finish(ax, title, sub=None, xlab=None, ylab=None):
    ax.set_title(title, loc="left", pad=24 if sub else 8)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=9,
                color=INK2, va="bottom")
    if xlab:
        ax.set_xlabel(xlab)
    if ylab:
        ax.set_ylabel(ylab)
    ax.set_axisbelow(True)


def fig_leakage(metrics: dict, outdir: Path) -> None:
    """
    The methods slide. A two-bar comparison is the right form: the job is a
    single magnitude contrast between two numbers, and the gap IS the message.
    """
    grouped = metrics["grouped_cv"]["r2"]
    random = metrics["grouped_cv"].get("random_cv_r2_INVALID")
    if random is None:
        print("  (skipping fig1: re-run train.py with --random-cv-demo)")
        return

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    vals = [random, grouped]
    labels = ["Random split\n(wrong)", "Grouped by lineage\n(correct)"]
    bars = ax.bar(labels, vals, width=0.5, color=[C_ORANGE, C_BLUE],
                  zorder=3)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.012, f"{v:.3f}",
                ha="center", va="bottom", fontweight="bold", color=INK)
    ax.annotate("", xy=(1, grouped), xytext=(1, random),
                arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.4))
    ax.text(1.06, (grouped + random) / 2,
            f"inflation\n+{random - grouped:.3f}", color=INK2, fontsize=9,
            va="center")
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_ylim(min(0, min(vals)) - 0.05, max(vals) * 1.25)
    _finish(ax, "Near-duplicate plasmids inflate performance",
            "Identical model, identical features — only the split differs",
            ylab="R² (out-of-fold)")
    fig.tight_layout()
    fig.savefig(outdir / "fig1_leakage.png", bbox_inches="tight")
    plt.close(fig)


def fig_label(df: pd.DataFrame, outdir: Path) -> None:
    """
    Small multiples, not one crowded panel: the components are on different
    scales, and the job is to show each distribution's shape, not to compare
    their magnitudes against each other.
    """
    comps = [c for c in ["group_size", "n_host_genera", "host_rank_max",
                         "n_countries", "n_compartments", "year_span"]
             if c in df.columns]
    lineages = df.drop_duplicates("group")
    n = len(comps)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.7 * nrow))
    axes = np.atleast_1d(axes).ravel()
    nice = {"group_size": "Records per lineage",
            "n_host_genera": "Distinct host genera",
            "host_rank_max": "Host-range rank (ordinal)",
            "n_countries": "Distinct countries",
            "n_compartments": "One-Health compartments",
            "year_span": "Year span of records"}
    for ax, c in zip(axes, comps):
        v = pd.to_numeric(lineages[c], errors="coerce").dropna()
        logx = c in ("group_size",)
        if logx:
            v = v[v > 0]
            bins = np.logspace(0, np.log10(max(v.max(), 10)), 26)
        else:
            bins = min(26, max(6, int(v.nunique())))
        ax.hist(v, bins=bins, color=C_BLUE, zorder=3)
        if logx:
            ax.set_xscale("log")
        ax.set_title(nice.get(c, c), loc="left", fontsize=10,
                     fontweight="normal", color=INK2)
        ax.set_ylabel("lineages", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.set_axisbelow(True)
    for ax in axes[n:]:
        ax.set_visible(False)
    fig.suptitle("What the dissemination label is made of", x=0.01, ha="left",
                 fontsize=12, fontweight="bold")
    fig.text(0.01, 0.955, f"One row per lineage (n = {len(lineages):,}); each "
             "component rank-transformed, then averaged",
             fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(outdir / "fig2_label.png", bbox_inches="tight")
    plt.close(fig)


def fig_performance(oof: pd.DataFrame, metrics: dict, outdir: Path) -> None:
    """Left: predicted vs observed. Right: the model against its baselines."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.4),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    x = oof.dissemination_true.to_numpy()
    y = oof.dissemination_oof.to_numpy()
    ax1.scatter(x, y, s=14, alpha=0.35, color=C_BLUE, linewidths=0, zorder=3)
    lo, hi = 0, 100
    ax1.plot([lo, hi], [lo, hi], color=MUTED, lw=1.2, ls="--", zorder=4)
    rho = spearmanr(x, y).statistic
    r2 = metrics["grouped_cv"]["r2"]
    ax1.text(0.03, 0.97, f"Spearman ρ = {rho:.3f}\nR² = {r2:.3f}",
             transform=ax1.transAxes, va="top", fontsize=10,
             bbox=dict(boxstyle="round,pad=0.45", fc=SURFACE, ec="#e6e5e1"))
    ax1.set_xlim(lo, hi)
    ax1.set_ylim(lo, hi)
    _finish(ax1, "Out-of-fold predictions",
            "Every point predicted by a model that never saw its lineage",
            "Observed dissemination (percentile)", "Predicted")

    names, vals = [], []
    for k, v in metrics["baselines"].items():
        names.append(k.replace("_", " "))
        vals.append(v["r2"])
    names.append("FULL MODEL")
    vals.append(r2)
    colors = [MUTED] * (len(vals) - 1) + [C_BLUE]
    ypos = np.arange(len(vals))
    ax2.barh(ypos, vals, color=colors, height=0.6, zorder=3)
    ax2.set_yticks(ypos, names)
    ax2.invert_yaxis()
    ax2.axvline(0, color=MUTED, lw=1)
    for i, v in enumerate(vals):
        ax2.text(v + (0.006 if v >= 0 else -0.006), i, f"{v:+.3f}",
                 va="center", ha="left" if v >= 0 else "right", fontsize=9,
                 color=INK)
    pad = max(0.05, abs(max(vals)) * 0.35)
    ax2.set_xlim(min(0, min(vals)) - pad, max(vals) + pad)
    _finish(ax2, "Against the baselines",
            "Same folds, same lineage weights", "R² (out-of-fold)")
    fig.tight_layout()
    fig.savefig(outdir / "fig3_performance.png", bbox_inches="tight")
    plt.close(fig)


def fig_importance(imp: pd.DataFrame, outdir: Path, top: int = 18) -> None:
    d = imp.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.4, 0.33 * len(d) + 1.6))
    ax.barh(np.arange(len(d)), d.importance, xerr=d.sd, height=0.65,
            color=C_BLUE, zorder=3,
            error_kw=dict(ecolor=MUTED, lw=1, capsize=2))
    ax.set_yticks(np.arange(len(d)), d.feature)
    ax.tick_params(axis="y", labelsize=9)
    _finish(ax, "What the model is reading",
            "Permutation importance (drop in R² when the feature is shuffled)",
            "Δ R²")
    fig.tight_layout()
    fig.savefig(outdir / "fig4_importance.png", bbox_inches="tight")
    plt.close(fig)


def fig_risk_space(df: pd.DataFrame, outdir: Path) -> None:
    """
    The product slide: risk is a plane, not a line.

    Two categories only (carbapenemase carriers vs the rest) because the
    all-pairs palette floors cap categorical hues at three, and because the
    contrast is the whole point.
    """
    d = df.dropna(subset=["dissemination", "cargo_severity"]).copy()
    if len(d) > 12000:
        d = d.sample(12000, random_state=0)
    carb = d.get("worst_class", pd.Series("", index=d.index)).eq("carbapenem")

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    jit = np.random.default_rng(0).normal(0, 0.006, len(d))
    ax.scatter(d.dissemination[~carb], (d.cargo_severity + jit)[~carb],
               s=9, alpha=0.22, color=C_BLUE, linewidths=0, zorder=3,
               label="other cargo")
    ax.scatter(d.dissemination[carb], (d.cargo_severity + jit)[carb],
               s=14, alpha=0.5, color=C_ORANGE, linewidths=0, zorder=4,
               label="carries a carbapenemase")

    # Iso-risk contours make the multiplication legible.
    xx, yy = np.meshgrid(np.linspace(0, 100, 200), np.linspace(0, 1, 200))
    risk = xx * (0.25 + 0.75 * yy)
    cs = ax.contour(xx, yy, risk, levels=[20, 40, 60, 80],
                    colors=MUTED, linewidths=0.9, linestyles=":", zorder=5)
    ax.clabel(cs, fmt="risk %d", fontsize=8, colors=MUTED)

    ax.set_xlim(0, 100)
    ax.set_ylim(-0.03, 1.05)
    leg = ax.legend(loc="lower left", frameon=True, fontsize=9)
    leg.get_frame().set_edgecolor("#e6e5e1")
    leg.get_frame().set_facecolor(SURFACE)
    _finish(ax, "Risk is a plane, not a line",
            "Same score, different meaning: an empty fast vehicle vs a "
            "stranded carbapenemase",
            "Observed dissemination (percentile)", "Cargo severity")
    fig.tight_layout()
    fig.savefig(outdir / "fig5_risk_space.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=Path("out"))
    args = ap.parse_args()

    run = args.out
    figdir = run / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(run / "dataset.csv.gz", index_col=0, low_memory=False)
    metrics = json.loads((run / "metrics.json").read_text())
    oof = pd.read_csv(run / "oof_predictions.csv", index_col=0)
    imp = pd.read_csv(run / "feature_importance.csv")

    train_rows = df[df.in_training_set] if "in_training_set" in df else df

    fig_leakage(metrics, figdir)
    fig_label(train_rows, figdir)
    fig_performance(oof, metrics, figdir)
    fig_importance(imp, figdir)
    fig_risk_space(df, figdir)

    for p in sorted(figdir.glob("*.png")):
        print(f"  wrote {p}")


if __name__ == "__main__":
    main()
