#!/usr/bin/env python3
"""
STEP 1 -- Build the modelling dataset.

Reads the PLSDB tables, groups plasmids into lineages, computes the
dissemination label and the cargo severity, assembles the backbone feature
matrix, and writes one tidy table to out/dataset.parquet (plus a CSV copy).

IMPORTANT DESIGN POINT -- read this before you touch --subset:

    Labels are ALWAYS computed on the FULL database, and only then is the
    training set subsetted.

    The label is a lineage-level property: "this family has been seen in 6
    genera across 40 assemblies". If you subset first, that same family looks
    like "2 records, 1 genus" and the signal you are trying to learn is
    destroyed. Subset last, never first.

Usage
-----
    python build_dataset.py                       # all plasmids
    python build_dataset.py --subset 1000         # 1000-plasmid hackathon set
    python build_dataset.py --subset 1000 --focus enterobacterales
    python build_dataset.py --include-stress      # add metal/biocide features
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from plrisk import cargo as cargo_mod
from plrisk import config, features, grouping, labels, metadata as meta_mod

# Genera used by --focus enterobacterales. Matched against the host proxy.
ENTEROBACTERALES = {
    "Escherichia", "Klebsiella", "Enterobacter", "Salmonella", "Citrobacter",
    "Serratia", "Proteus", "Providencia", "Morganella", "Shigella",
    "Raoultella", "Cronobacter", "Hafnia", "Edwardsiella", "Yersinia",
    "Kluyvera", "Leclercia", "Pantoea", "Pluralibacter", "Lelliottia",
}


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def load_inputs() -> dict[str, pd.DataFrame]:
    log("Loading input tables")
    missing = [p for p in (config.TYPING, config.AMR) if not Path(p).exists()]
    if missing:
        raise SystemExit(
            "ERROR: required file(s) not found:\n  "
            + "\n  ".join(str(m) for m in missing)
            + "\nSet the paths in plrisk/config.py or via PLRISK_* env vars."
        )

    typing = pd.read_csv(config.TYPING, low_memory=False)
    log(f"  typing.csv           : {len(typing):>7,} rows")

    amr = pd.read_csv(config.AMR, sep="\t", low_memory=False)
    log(f"  amr.tsv              : {len(amr):>7,} rows")

    pf = None
    if Path(config.PLASMIDFIND).exists():
        pf = pd.read_csv(config.PLASMIDFIND, low_memory=False)
        log(f"  plasmidfinder.csv    : {len(pf):>7,} rows")
    else:
        log("  plasmidfinder.csv    : not found (Inc features will be skipped)")

    return {"typing": typing, "amr": amr, "plasmidfinder": pf}


def subset_for_training(df: pd.DataFrame, n: int | None, focus: str | None,
                        max_per_group: int = config.MAX_MEMBERS_PER_GROUP,
                        seed: int = config.RANDOM_STATE) -> pd.DataFrame:
    """
    Choose which plasmids to train on.

    STRATEGY -- four steps, each of which is worth defending in the talk:

    1. DROP SINGLETON LINEAGES. A group seen exactly once has a dissemination
       label that is mostly sampling accident. They can still be SCORED later;
       this filter applies only to training.
       Note the threshold is 2, not 5 or 10: group size is one third of the
       label, so an aggressive filter here is range restriction on the target
       in disguise. Filtering at >=3 pushes the training label distribution
       into the top quartile and leaves nothing to learn.

    2. CAP MEMBERS PER LINEAGE. PLSDB has lineages with 1,700 near-identical
       records. Without a cap, a handful of epidemic plasmids supply most
       training rows, the loss is dominated by memorising them, and the
       plasmid-level label distribution is dragged to the ceiling. We keep at
       most `max_per_group` members, chosen at random within the lineage.

    3. OPTIONAL TAXONOMIC FOCUS. For an AMR-epidemiology demo, restricting to
       Enterobacterales concentrates the signal: that is where the clinically
       relevant conjugative plasmids are, and it removes cyanobacterial and
       plant-symbiont plasmids whose "dissemination" means something quite
       different. It also makes the host-range label comparable across rows.

    4. SAMPLE WHOLE GROUPS, STRATIFIED BY LABEL DECILE. Sampling groups (never
       individual plasmids) keeps CV folds clean. Stratifying across the label
       range stops a uniform draw from returning mostly mid-range lineages.

    Returns the subset; the caller marks these rows with `in_training_set`.
    """
    rng = np.random.default_rng(seed)
    d = df.copy()

    # --- 1. drop singleton lineages ---------------------------------------
    before = len(d)
    d = d[d.group_size >= config.MIN_GROUP_SIZE_FOR_TRAINING]
    log(f"  dropped {before - len(d):,} plasmids in lineages smaller than "
        f"{config.MIN_GROUP_SIZE_FOR_TRAINING} (kept {len(d):,})")

    # --- 2. optional taxonomic focus --------------------------------------
    if focus == "enterobacterales":
        before = len(d)
        # Prefer the real NCBI order from the metadata join; fall back to a
        # genus whitelist against the MOB-suite proxy when nuccore.csv is
        # absent. The real taxonomy is both more accurate and more complete.
        if "host_order" in d.columns and d["host_order"].notna().any():
            mask = d["host_order"].astype(str).str.lower().eq("enterobacterales")
            how = "NCBI order"
        else:
            mask = d["host_genus"].isin(ENTEROBACTERALES)
            how = "genus whitelist (proxy)"
        d = d[mask]
        log(f"  --focus enterobacterales: kept {len(d):,} of {before:,} "
            f"via {how}")

    # --- 3. cap members per lineage ---------------------------------------
    before = len(d)
    d = (d.assign(_r=rng.random(len(d)))
           .sort_values("_r")
           .groupby("group", sort=False)
           .head(max_per_group)
           .drop(columns="_r"))
    log(f"  capped at {max_per_group} members per lineage: "
        f"{len(d):,} rows from {before:,} "
        f"({d.group.nunique():,} lineages)")

    if n is None or n >= len(d):
        log(f"  no further subsetting ({len(d):,} plasmids available)")
        return d

    # --- 4. stratified sampling of whole groups ---------------------------
    gl = d.groupby("group").dissemination.first()
    sizes = d.groupby("group").size()
    n_strata = 10
    strata = pd.qcut(gl.rank(method="first"), q=n_strata, labels=False)

    chosen: list[str] = []
    target_per = n / n_strata
    for s in range(n_strata):
        pool = gl.index[strata == s].to_numpy()
        rng.shuffle(pool)
        acc = 0
        for g in pool:
            if acc >= target_per:
                break
            chosen.append(g)
            acc += int(sizes[g])

    out = d[d.group.isin(set(chosen))]
    log(f"  subset               : {len(out):,} plasmids in "
        f"{out.group.nunique():,} lineages across {n_strata} label strata")
    log(f"  label spread in subset: mean {out.dissemination.mean():.1f}, "
        f"sd {out.dissemination.std():.1f}, "
        f"range {out.dissemination.min():.1f}-{out.dissemination.max():.1f}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subset", type=int, default=None,
                    help="approximate number of plasmids to mark for training")
    ap.add_argument("--focus", choices=["enterobacterales"], default=None,
                    help="restrict the training set to a taxonomic slice")
    ap.add_argument("--max-per-group", type=int,
                    default=config.MAX_MEMBERS_PER_GROUP,
                    help="cap on training rows contributed by one lineage")
    ap.add_argument("--label", choices=["full", "taxonomic"], default="full",
                    help="'full' uses country/One-Health/time components when "
                         "nuccore.csv is present; 'taxonomic' uses only "
                         "lineage size, host genera and host range")
    ap.add_argument("--include-stress", action="store_true",
                    help="add metal/biocide resistance gene counts as features")
    ap.add_argument("--out", type=Path, default=config.OUT_DIR)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    typing = data["typing"]

    log("\nGrouping plasmids into lineages")
    groups = grouping.assign_groups(
        typing, config.MASH_SIM,
        config.PLING if Path(config.PLING).exists() else None)

    log("\nJoining metadata")
    md = None
    if Path(config.NUCCORE).exists():
        md = meta_mod.load_metadata(config.NUCCORE, config.TAXONOMY,
                                    config.BIOSAMPLE)
    else:
        log(f"  nuccore.csv not found at {config.NUCCORE}")
        log("  -> label falls back to taxonomic breadth only (no country, no")
        log("     One-Health compartment, no time span). See README section 2.3.")

    log("\nBuilding dissemination labels (on the FULL database)")
    lab = labels.build_labels(typing, groups, md,
                              use_metadata_components=(args.label == "full"))
    log(f"  labelled             : {len(lab):,} plasmids in "
        f"{lab.group.nunique():,} groups")

    log("\nComputing cargo severity from AMRFinderPlus")
    carg = cargo_mod.cargo_table(data["amr"])
    log(f"  plasmids with >=1 AMR/stress hit: {len(carg):,}")

    log("\nBuilding backbone features")
    F = features.build_features(
        typing, data["plasmidfinder"], carg,
        include_stress=args.include_stress)

    # ---- assemble ---------------------------------------------------------
    log("\nAssembling")
    lab_cols = [c for c in ["group", "group_size", "n_host_genera",
                            "n_host_species", "host_rank_max", "n_countries",
                            "n_compartments", "year_span", "first_year",
                            "dissemination", "host_genus", "host_genus_source"]
                if c in lab.columns]
    df = F.join(lab[lab_cols], how="inner")
    df = df.join(groups[["group_source"]])

    # Per-plasmid metadata, carried for filtering, display and plotting.
    # NOT features -- the feature matrix was built and frozen above, and the
    # leakage guard has already run on it.
    if md is not None:
        for c in ["host_species", "host_order", "host_family", "country",
                  "compartments", "year", "ecosystem_tags", "DISEASE_group"]:
            if c in md.columns:
                df[c] = md[c].reindex(df.index)

    carg_cols = ["cargo_severity", "n_amr_genes", "n_amr_classes",
                 "n_stress_genes", "worst_class", "amr_gene_list"]
    cj = carg.reindex(df.index)
    for c in carg_cols:
        df[c] = cj[c]
    df["cargo_severity"] = df["cargo_severity"].fillna(0.0)
    for c in ["n_amr_genes", "n_amr_classes", "n_stress_genes"]:
        df[c] = df[c].fillna(0).astype(int)
    df["worst_class"] = df["worst_class"].fillna("none")
    df["amr_gene_list"] = df["amr_gene_list"].fillna("")

    # Host proxy kept as a plain column for --focus and for plotting. It is
    # NOT a feature -- the feature matrix was already built and frozen above.
    tt = typing.drop_duplicates("NUCCORE_ACC").set_index(
        typing.drop_duplicates("NUCCORE_ACC")["NUCCORE_ACC"].astype(str))
    df["host_genus_proxy"] = (tt["mash_neighbor_identification"]
                              .reindex(df.index).map(labels._genus_of))
    if "host_genus" not in df.columns:
        df["host_genus"] = df["host_genus_proxy"]

    # Raw descriptive columns carried through for display in score.py and for
    # plotting. These are NOT features -- the feature matrix was built and
    # frozen above, and the leakage guard already ran on it.
    for src, dst in [("size", "size"),
                     ("predicted_mobility", "predicted_mobility"),
                     ("rep_type(s)", "rep_type_profile"),
                     ("relaxase_type(s)", "relaxase_profile"),
                     ("mash_nearest_neighbor", "mash_nearest_neighbor")]:
        if src in tt.columns:
            df[dst] = tt[src].reindex(df.index)

    # ---- mark the training subset ----------------------------------------
    log("\nSelecting training subset")
    train = subset_for_training(df, args.subset, args.focus,
                                max_per_group=args.max_per_group)
    df["in_training_set"] = df.index.isin(train.index)

    feature_cols = list(F.columns)
    pd.Series(feature_cols).to_csv(args.out / "feature_columns.txt",
                                   index=False, header=False)

    out_parquet = args.out / "dataset.parquet"
    try:
        df.to_parquet(out_parquet)
        log(f"\nWrote {out_parquet}  ({len(df):,} x {df.shape[1]:,})")
    except Exception as e:                       # pyarrow not installed
        log(f"\n(parquet unavailable: {e}; writing CSV only)")
    df.to_csv(args.out / "dataset.csv.gz", compression="gzip")
    log(f"Wrote {args.out / 'dataset.csv.gz'}")
    log(f"Wrote {args.out / 'feature_columns.txt'}  ({len(feature_cols)} features)")

    # ---- summary ----------------------------------------------------------
    log("\n" + "=" * 62)
    log("SUMMARY")
    log("=" * 62)
    log(f"  plasmids total            : {len(df):,}")
    log(f"  marked for training       : {int(df.in_training_set.sum()):,}")
    log(f"  groups (all)              : {df.group.nunique():,}")
    log(f"  groups (training)         : {df[df.in_training_set].group.nunique():,}")
    log(f"  features                  : {len(feature_cols):,}")
    log(f"  plasmids with AMR cargo   : {int((df.cargo_severity > 0).sum()):,}")
    log(f"  mean cargo severity       : {df.cargo_severity.mean():.3f}")
    log("\n  dissemination label distribution (training rows):")
    tr = df[df.in_training_set]
    if len(tr):
        log("   " + str(tr.dissemination.describe().round(2).to_dict()))
        log("\n  group source breakdown (training rows):")
        for k, v in tr.group_source.value_counts().items():
            log(f"    {k:<14} {v:,}")


if __name__ == "__main__":
    main()
