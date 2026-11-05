#!/usr/bin/env python3
"""
Diagnostic -- are we data-limited or signal-limited?

The question matters for how you spend the remaining hours. If performance
keeps climbing with more lineages, go and get more lineages. If it plateaus,
more data is wasted effort and the next move is better features or a better
label.

Naively comparing the `--subset 1000` run against the full run does NOT answer
it, because train.py scales model capacity to the number of lineages: the two
runs use different models, so any difference confounds capacity with data
volume. This script holds the model FIXED and varies only the number of
training lineages, which is the comparison you actually want.

Subsampling is by LINEAGE, never by plasmid, for the usual reason.

Usage
-----
    python learning_curve.py --dataset sweep_taxonomic_all/dataset.csv.gz \
                             --features sweep_taxonomic_all/feature_columns.txt
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import GroupKFold

from plrisk import config

CAPACITIES = {
    "small":  dict(max_iter=120, max_leaf_nodes=7,  min_samples_leaf=30,
                   l2_regularization=5.0),
    "medium": dict(max_iter=250, max_leaf_nodes=15, min_samples_leaf=20,
                   l2_regularization=2.0),
    "large":  dict(max_iter=400, max_leaf_nodes=31, min_samples_leaf=20,
                   l2_regularization=1.0),
}


def model(cap: str) -> HistGradientBoostingRegressor:
    return HistGradientBoostingRegressor(
        loss="squared_error", learning_rate=0.05, early_stopping=False,
        random_state=config.RANDOM_STATE, **CAPACITIES[cap])


def grouped_cv(X, y, groups, w, cap: str, folds: int = 5):
    oof = np.full(len(y), np.nan)
    for tr, te in GroupKFold(n_splits=folds).split(X, y, groups):
        m = model(cap)
        m.fit(X.iloc[tr], y.iloc[tr], sample_weight=w.iloc[tr])
        oof[te] = m.predict(X.iloc[te])
    return (r2_score(y, oof, sample_weight=w), spearmanr(y, oof).statistic)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", type=Path, required=True)
    ap.add_argument("--features", type=Path, required=True)
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[200, 400, 800, 1600, 3200, 99999])
    ap.add_argument("--repeats", type=int, default=3,
                    help="random lineage subsamples per size, to show spread")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.dataset, index_col=0, low_memory=False)
    feats = [l.strip() for l in open(args.features) if l.strip()]
    tr = df[df.in_training_set].copy()

    all_groups = tr["group"].astype(str).unique()
    print(f"available lineages: {len(all_groups):,}  "
          f"(rows: {len(tr):,})\n")

    header = f"{'lineages':>9} {'rows':>7}  " + "  ".join(
        f"{c:>16}" for c in CAPACITIES)
    print(header)
    print("-" * len(header))

    rows = []
    for size in args.sizes:
        n = min(size, len(all_groups))
        cells = {}
        for cap in CAPACITIES:
            scores = []
            reps = 1 if n == len(all_groups) else args.repeats
            for r in range(reps):
                rng = np.random.default_rng(config.RANDOM_STATE + r)
                pick = (all_groups if n == len(all_groups)
                        else rng.choice(all_groups, n, replace=False))
                sub = tr[tr.group.isin(set(pick))]
                X = sub[feats].astype(float)
                y = sub["dissemination"].astype(float)
                g = sub["group"].astype(str)
                w = 1.0 / g.map(g.value_counts())
                scores.append(grouped_cv(X, y, g, w, cap))
            r2 = np.mean([s[0] for s in scores])
            rho = np.mean([s[1] for s in scores])
            sd = np.std([s[0] for s in scores])
            cells[cap] = (r2, rho, sd)
            rows.append({"n_lineages": n, "capacity": cap,
                         "r2": r2, "spearman": rho, "r2_sd": sd})
        n_rows = len(tr[tr.group.isin(set(all_groups[:n]))])
        print(f"{n:>9,} {n_rows:>7,}  " + "  ".join(
            f"{cells[c][0]:>+7.3f}/{cells[c][1]:>5.3f}" + " " * 3
            for c in CAPACITIES))
        if n == len(all_groups):
            break

    print("\n(each cell is R2 / Spearman rho, grouped 5-fold CV)")
    best = max(rows, key=lambda r: r["r2"])
    print(f"\nbest: {best['capacity']} capacity at {best['n_lineages']:,} "
          f"lineages -> R2 {best['r2']:+.3f}, rho {best['spearman']:.3f}")

    if args.out:
        args.out.write_text(json.dumps(rows, indent=2))
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
