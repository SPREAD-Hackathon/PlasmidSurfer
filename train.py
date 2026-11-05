#!/usr/bin/env python3
"""
STEP 2 -- Train and honestly evaluate the dissemination model.

WHAT IS BEING LEARNED
---------------------
    f(backbone features) -> dissemination potential (0-100)

and nothing else. Cargo severity is not learned, it is looked up (plrisk/cargo.py)
and multiplied in at scoring time.

THE THREE THINGS THAT MAKE THIS EVALUATION CREDIBLE
----------------------------------------------------
1. GROUPED CROSS-VALIDATION. Folds are split by lineage, so a plasmid and its
   40 near-identical siblings are always on the same side. Random CV here
   inflates R^2 from ~0.3 to ~0.9 and means nothing. We print both so you can
   show the gap on a slide -- it is the single most convincing methods point
   you can make in a hackathon room.

2. SAMPLE WEIGHTS OF 1/n_members. Each lineage contributes total weight 1,
   whatever its size. Without this the loss is dominated by whichever epidemic
   plasmid happens to be most sequenced.

3. TWO BASELINES. A mean predictor (floor) and a size-only predictor. If the
   full model cannot beat "how big is the plasmid", you have not built a
   plasmid model, you have built a ruler. Report the comparison either way.

Usage
-----
    python train.py
    python train.py --random-cv-demo        # also print the inflated numbers
    python train.py --dataset out/dataset.csv.gz --out out
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold, KFold

from plrisk import config


def log(m: str) -> None:
    print(m, file=sys.stderr, flush=True)


def load_dataset(path: Path, feature_cols_path: Path):
    if path.suffix == ".parquet":
        df = pd.read_csv(path) if not path.exists() else pd.read_parquet(path)
    else:
        df = pd.read_csv(path, index_col=0, low_memory=False)
    feats = [l.strip() for l in open(feature_cols_path) if l.strip()]
    missing = [c for c in feats if c not in df.columns]
    if missing:
        raise SystemExit(f"Feature columns missing from dataset: {missing[:10]}")
    return df, feats


def make_model(n_groups: int | None = None) -> HistGradientBoostingRegressor:
    """
    Gradient-boosted trees from scikit-learn (no XGBoost/LightGBM dependency,
    so this runs on any laptop in the room).

    !! EARLY STOPPING IS OFF ON PURPOSE -- THIS WAS A REAL BUG !!

    HistGradientBoostingRegressor's built-in early stopping holds out a
    RANDOM `validation_fraction` of the training rows. In this dataset that is
    catastrophic: near-identical siblings from the same lineage land on both
    sides of that internal split, so the internal validation score keeps
    improving while true generalisation gets worse, early stopping never
    fires, and the model memorises. The first version of this file had
    `early_stopping=True` and produced R2 = -0.30 on grouped CV -- worse than
    predicting the mean.

    sklearn has no way to pass a grouped internal validation split, so instead
    early stopping is disabled and capacity is bounded directly, scaled to the
    number of independent LINEAGES (not rows -- rows are inflated by
    near-duplicates). If you add far more lineages, raise these and re-check
    the gap between grouped and random CV.
    """
    # Capacity is scaled to the number of independent LINEAGES, not rows --
    # rows are inflated by near-duplicates and would license a far bigger
    # model than the data supports.
    #
    # These two rungs are measured, not guessed. `learning_curve.py` holds
    # capacity fixed and varies lineage count; a third, larger rung
    # (400 iters / 31 leaves / l2 1.0) was tested and is WORSE at every
    # training size, from 200 lineages (R2 -0.140 vs +0.061) to the full
    # 4,425 (+0.128 vs +0.168). It is not included. Re-run that script
    # before adding capacity back.
    n = n_groups or 1000
    if n < 800:
        max_iter, leaves, leaf_min, l2 = 120, 7, 30, 5.0
    else:
        max_iter, leaves, leaf_min, l2 = 250, 15, 20, 2.0

    return HistGradientBoostingRegressor(
        loss="squared_error",
        max_iter=max_iter,
        learning_rate=0.05,
        max_leaf_nodes=leaves,
        min_samples_leaf=leaf_min,
        l2_regularization=l2,
        early_stopping=False,          # see the docstring -- do not re-enable
        random_state=config.RANDOM_STATE,
    )


def evaluate(X, y, groups, w, n_splits: int, grouped: bool = True,
             n_groups: int | None = None):
    """Run CV and return out-of-fold predictions plus per-fold metrics."""
    splitter = (GroupKFold(n_splits=n_splits) if grouped
                else KFold(n_splits=n_splits, shuffle=True,
                           random_state=config.RANDOM_STATE))
    oof = np.full(len(y), np.nan)
    rows = []
    split_args = (X, y, groups) if grouped else (X, y)
    for k, (tr, te) in enumerate(splitter.split(*split_args), start=1):
        m = make_model(n_groups)
        m.fit(X.iloc[tr], y.iloc[tr], sample_weight=w.iloc[tr])
        p = m.predict(X.iloc[te])
        oof[te] = p
        rows.append({
            "fold": k,
            "n_test": len(te),
            "r2": r2_score(y.iloc[te], p, sample_weight=w.iloc[te]),
            "mae": mean_absolute_error(y.iloc[te], p, sample_weight=w.iloc[te]),
            "spearman": spearmanr(y.iloc[te], p).statistic,
        })
    return oof, pd.DataFrame(rows)


def baseline_scores(X, y, groups, w, n_splits, n_groups=None):
    """Mean predictor and size-only predictor, same folds, same weights."""
    out = {}
    gkf = GroupKFold(n_splits=n_splits)

    for name, cols, est in [
        ("mean_only", None, DummyRegressor(strategy="mean")),
        ("size_only", ["log_size"], make_model(n_groups)),
        ("size_gc_mobility", [c for c in X.columns
                              if c in ("log_size", "gc")
                              or c.startswith("mobility__")], make_model(n_groups)),
    ]:
        preds = np.full(len(y), np.nan)
        for tr, te in gkf.split(X, y, groups):
            Xi = X if cols is None else X[cols]
            e = est.__class__(**est.get_params())
            e.fit(Xi.iloc[tr], y.iloc[tr], sample_weight=w.iloc[tr])
            preds[te] = e.predict(Xi.iloc[te])
        out[name] = {
            "r2": r2_score(y, preds, sample_weight=w),
            "mae": mean_absolute_error(y, preds, sample_weight=w),
            "spearman": spearmanr(y, preds).statistic,
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", type=Path, default=config.OUT_DIR / "dataset.csv.gz")
    ap.add_argument("--features", type=Path,
                    default=config.OUT_DIR / "feature_columns.txt")
    ap.add_argument("--out", type=Path, default=config.OUT_DIR)
    ap.add_argument("--folds", type=int, default=config.N_FOLDS)
    ap.add_argument("--random-cv-demo", action="store_true",
                    help="also report the (invalid) random-split numbers, to "
                         "show the size of the leakage effect")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    df, feats = load_dataset(args.dataset, args.features)
    tr = df[df.in_training_set].copy()
    if not len(tr):
        raise SystemExit("No rows flagged in_training_set. Re-run build_dataset.py.")

    X = tr[feats].astype(float)
    y = tr["dissemination"].astype(float)
    groups = tr["group"].astype(str)

    # Each lineage gets total weight 1, shared among the rows it contributed.
    counts = groups.map(groups.value_counts())
    w = 1.0 / counts

    log("=" * 64)
    log("TRAINING")
    log("=" * 64)
    log(f"  rows (plasmids)     : {len(tr):,}")
    log(f"  lineages (groups)   : {groups.nunique():,}")
    log(f"  features            : {len(feats):,}")
    log(f"  label  mean {y.mean():.1f}  sd {y.std():.1f}  "
        f"range {y.min():.1f}-{y.max():.1f}")
    log(f"  folds               : {args.folds} (GroupKFold on lineage)")

    # ---- honest, grouped CV ----------------------------------------------
    log("\nGrouped cross-validation (the real numbers)")
    n_groups = groups.nunique()
    oof, folds = evaluate(X, y, groups, w, args.folds, grouped=True,
                          n_groups=n_groups)
    log(folds.round(3).to_string(index=False))
    overall = {
        "r2": float(r2_score(y, oof, sample_weight=w)),
        "mae": float(mean_absolute_error(y, oof, sample_weight=w)),
        "spearman": float(spearmanr(y, oof).statistic),
    }
    log(f"\n  OOF  R2 = {overall['r2']:.3f}   MAE = {overall['mae']:.2f}   "
        f"Spearman rho = {overall['spearman']:.3f}")

    # ---- baselines --------------------------------------------------------
    log("\nBaselines (same folds, same weights)")
    bl = baseline_scores(X, y, groups, w, args.folds, n_groups)
    for k, v in bl.items():
        log(f"  {k:<18} R2 = {v['r2']:+.3f}   MAE = {v['mae']:5.2f}   "
            f"rho = {v['spearman']:+.3f}")
    log(f"  {'FULL MODEL':<18} R2 = {overall['r2']:+.3f}   "
        f"MAE = {overall['mae']:5.2f}   rho = {overall['spearman']:+.3f}")

    gain = overall["r2"] - bl["size_only"]["r2"]
    if gain < 0.02:
        log("\n  !! The full model barely beats plasmid size alone. Do not "
            "\n     oversell this. Either the backbone features are not "
            "\n     informative at this sample size, or the label is "
            "\n     dominated by something size-correlated. Say so.")
    else:
        log(f"\n  Full model adds R2 +{gain:.3f} over size alone.")

    # ---- the leakage demonstration ---------------------------------------
    if args.random_cv_demo:
        log("\nRandom (NON-grouped) CV -- THIS IS THE WRONG WAY, shown for "
            "comparison")
        oof_r, folds_r = evaluate(X, y, groups, w, args.folds, grouped=False,
                                  n_groups=n_groups)
        r2r = r2_score(y, oof_r, sample_weight=w)
        log(f"  random-split R2 = {r2r:.3f}  vs grouped R2 = {overall['r2']:.3f}"
            f"   (inflation = {r2r - overall['r2']:+.3f})")
        overall["random_cv_r2_INVALID"] = float(r2r)

    # ---- fit the final model on everything --------------------------------
    log("\nFitting final model on all training rows")
    final = make_model(n_groups)
    final.fit(X, y, sample_weight=w)

    # Permutation importance is the slowest step by far: it refits nothing
    # but re-scores the model once per feature per repeat. With 374 features
    # that is thousands of passes, so both the repeat count and the number of
    # rows scored are capped. The ranking is stable well below the full data;
    # the reported SDs let you see when it is not.
    n_rep = 8 if len(X) <= 3000 else 5
    max_s = min(len(X), 4000)
    log(f"Computing permutation importance "
        f"({n_rep} repeats on {max_s:,} rows; the slowest step)")
    pi = permutation_importance(final, X, y, n_repeats=n_rep,
                                max_samples=max_s,
                                random_state=config.RANDOM_STATE,
                                sample_weight=w, scoring="r2",
                                n_jobs=-1)
    imp = (pd.DataFrame({"feature": X.columns,
                         "importance": pi.importances_mean,
                         "sd": pi.importances_std})
             .sort_values("importance", ascending=False))
    imp.to_csv(args.out / "feature_importance.csv", index=False)
    log("\n  Top 20 features by permutation importance:")
    for _, r in imp.head(20).iterrows():
        log(f"    {r.feature:<42} {r.importance:+.4f} +/- {r.sd:.4f}")

    # ---- persist ----------------------------------------------------------
    with open(args.out / "model.pkl", "wb") as fh:
        pickle.dump({"model": final, "features": feats,
                     "label": "dissemination", "cargo_floor": config.CARGO_FLOOR},
                    fh)
    oof_df = pd.DataFrame({"dissemination_true": y, "dissemination_oof": oof,
                           "group": groups, "weight": w}, index=tr.index)
    oof_df.to_csv(args.out / "oof_predictions.csv")
    with open(args.out / "metrics.json", "w") as fh:
        json.dump({"grouped_cv": overall, "baselines": bl,
                   "n_rows": int(len(tr)), "n_groups": int(groups.nunique()),
                   "n_features": len(feats),
                   "per_fold": folds.to_dict("records")}, fh, indent=2)

    log(f"\nWrote {args.out/'model.pkl'}, {args.out/'metrics.json'}, "
        f"{args.out/'feature_importance.csv'}, {args.out/'oof_predictions.csv'}")


if __name__ == "__main__":
    main()
