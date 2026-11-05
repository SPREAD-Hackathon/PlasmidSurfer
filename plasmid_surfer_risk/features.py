"""
Backbone feature matrix.

THE ONE RULE THIS MODULE ENFORCES
----------------------------------
Nothing that went into the LABEL may come out as a FEATURE, and no resistance
cargo may be used to predict dissemination. If either leaks, the model becomes
a tautology and the whole project collapses under the first good question at
the poster.

BANNED, and actively asserted against at the bottom of this file:

  * primary_cluster_id / secondary_cluster_id  -- these ARE the grouping, so
    they encode group size, which is a third of the label.
  * mash_nearest_neighbor / mash_neighbor_identification -- the host proxy used
    to build the label.
  * observed_host_range_* / predicted_host_range_* -- literally the label.
  * every AMR gene from amr.tsv -- this is the cargo term, which is multiplied
    in afterwards. If the model could see it, it would learn "carbapenemase =
    dangerous" and we would have built a very expensive abricate.

WHAT IS ALLOWED: the mobilisation and maintenance machinery, i.e. the vehicle
rather than its payload.

  size, GC, contig count         basic sequence properties
  replicon types (rep_cluster)   from MOB-suite; the backbone's identity
  Inc groups                     from PlasmidFinder; the classical typing
  relaxase (MOB) type            conjugative machinery
  MPF type                       mating pair formation system
  oriT presence                  origin of transfer
  predicted_mobility             MOB-suite's conjugative/mobilisable/non- call
  PMLST scheme present           proxy for being a well-studied epidemic type

ARGUABLE CASE: metal and biocide resistance (copper, mercury, QAC). These are
co-selection drivers and genuinely predict persistence, but they come from the
same file as the antibiotic cargo. They are OFF by default and enabled with
--include-stress. Run it both ways and report the difference; if it barely
moves, say so and leave it off.
"""

from __future__ import annotations

import re
import sys

import numpy as np
import pandas as pd

from . import config

# Columns from typing.csv that must never become features.
BANNED_TYPING_COLS = {
    "primary_cluster_id", "secondary_cluster_id",
    "mash_nearest_neighbor", "mash_neighbor_distance",
    "mash_neighbor_identification",
    "observed_host_range_ncbi_rank", "observed_host_range_ncbi_name",
    "predicted_host_range_overall_rank", "predicted_host_range_overall_name",
    "reported_host_range_lit_rank", "reported_host_range_lit_name",
    "md5", "associated_pmid(s)", "NUCCORE_ACC",
}


def _split_multi(cell) -> list[str]:
    """MOB-suite packs multi-valued fields as comma-separated strings."""
    if not isinstance(cell, str) or not cell.strip() or cell.strip().lower() == "nan":
        return []
    return [p.strip() for p in cell.split(",") if p.strip()]


def _multihot(series: pd.Series, prefix: str, min_count: int) -> pd.DataFrame:
    """One-hot a multi-valued column, keeping only levels seen >= min_count."""
    lists = series.map(_split_multi)
    counts: dict[str, int] = {}
    for vals in lists:
        for v in set(vals):
            counts[v] = counts.get(v, 0) + 1
    keep = sorted(v for v, c in counts.items() if c >= min_count)
    if not keep:
        return pd.DataFrame(index=series.index)
    idx = {v: i for i, v in enumerate(keep)}
    mat = np.zeros((len(series), len(keep)), dtype=np.int8)
    for r, vals in enumerate(lists):
        for v in set(vals):
            j = idx.get(v)
            if j is not None:
                mat[r, j] = 1
    cols = [f"{prefix}__{re.sub(r'[^0-9A-Za-z_]+', '_', v)}" for v in keep]
    return pd.DataFrame(mat, index=series.index, columns=cols)


def inc_family(replicon: str) -> str:
    """
    Collapse a PlasmidFinder allele to its Inc family.

    'IncFII(pRSB107)_1__AJ851089' -> 'IncFII'
    'Col(pHAD28)'                 -> 'Col'
    Keeps the biologically meaningful level and stops 277 allele labels from
    becoming 277 almost-empty columns.
    """
    s = str(replicon)
    s = re.split(r"_\d+__", s)[0]           # drop allele + accession suffix
    s = re.sub(r"\(.*?\)", "", s)           # drop parenthetical variant
    s = re.sub(r"_\d+$", "", s)
    return s.strip() or "unknown"


def build_features(typing: pd.DataFrame,
                   plasmidfinder: pd.DataFrame | None = None,
                   cargo: pd.DataFrame | None = None,
                   include_stress: bool = False,
                   min_count: int = config.MIN_CATEGORY_COUNT,
                   verbose: bool = True) -> pd.DataFrame:
    """Assemble the backbone feature matrix, indexed by NUCCORE_ACC."""
    t = typing.copy()
    t["NUCCORE_ACC"] = t["NUCCORE_ACC"].astype(str)
    t = t.drop_duplicates("NUCCORE_ACC").set_index("NUCCORE_ACC")

    F = pd.DataFrame(index=t.index)

    # --- continuous sequence properties -----------------------------------
    F["log_size"] = np.log10(pd.to_numeric(t["size"], errors="coerce").clip(lower=1))
    F["gc"] = pd.to_numeric(t["gc"], errors="coerce")
    F["num_contigs"] = pd.to_numeric(t["num_contigs"], errors="coerce").fillna(1)
    # A plasmid assembled in many contigs is a data-quality signal, not biology.
    # Kept so the model can discount fragmented records rather than be fooled.
    F["is_fragmented"] = (F["num_contigs"] > 1).astype(int)

    # --- mobility machinery -----------------------------------------------
    F["n_rep_types"] = t["rep_type(s)"].map(lambda c: len(_split_multi(c)))
    F["n_relaxases"] = t["relaxase_type(s)"].map(lambda c: len(_split_multi(c)))
    F["has_relaxase"] = (F["n_relaxases"] > 0).astype(int)
    F["has_mpf"] = t["mpf_type"].notna().astype(int)
    F["has_orit"] = t["orit_type(s)"].notna().astype(int)
    F["has_pmlst"] = t["PMLST_scheme"].notna().astype(int)

    mob = t["predicted_mobility"].astype(str).str.strip().str.lower()
    for lvl in ["conjugative", "mobilizable", "non-mobilizable"]:
        F[f"mobility__{lvl.replace('-', '_')}"] = (mob == lvl).astype(int)

    # --- categorical machinery, multi-hot ----------------------------------
    F = F.join(_multihot(t["relaxase_type(s)"], "mob", min_count))
    F = F.join(_multihot(t["mpf_type"], "mpf", min_count))
    F = F.join(_multihot(t["orit_type(s)"], "orit", min_count))
    F = F.join(_multihot(t["rep_type(s)"], "rep", min_count))

    # --- PlasmidFinder Inc typing ------------------------------------------
    if plasmidfinder is not None and len(plasmidfinder):
        pf = plasmidfinder.copy()
        pf["NUCCORE_ACC"] = pf["NUCCORE_ACC"].astype(str)
        pf["fam"] = pf["typing"].map(inc_family)
        agg = (pf.groupby("NUCCORE_ACC")["fam"]
                 .apply(lambda s: ",".join(sorted(set(s)))))
        agg = agg.reindex(F.index).fillna("")
        F = F.join(_multihot(agg, "inc", min_count))
        F["n_inc_groups"] = agg.map(lambda c: len(_split_multi(c)))

    # --- optional co-selection features ------------------------------------
    if include_stress and cargo is not None:
        F["n_stress_genes"] = (cargo["n_stress_genes"]
                               .reindex(F.index).fillna(0).astype(int))
        F["has_stress"] = (F["n_stress_genes"] > 0).astype(int)

    F = F.fillna(0.0)

    # --- the assertion that protects the whole project ---------------------
    _assert_no_leakage(F, include_stress=include_stress)

    if verbose:
        print(f"  Feature matrix       : {F.shape[0]:,} plasmids "
              f"x {F.shape[1]:,} features", file=sys.stderr)
    return F


def _assert_no_leakage(F: pd.DataFrame, include_stress: bool) -> None:
    """Fail loudly rather than silently training a tautology."""
    bad = []
    forbidden_tokens = [
        "cluster_id", "host_range", "mash_neighbor", "mash_nearest",
        "dissemination", "group_size", "cargo", "risk", "n_amr", "worst_class",
    ]
    for col in F.columns:
        lc = col.lower()
        if any(tok in lc for tok in forbidden_tokens):
            bad.append(col)
    if not include_stress:
        bad += [c for c in F.columns if "stress" in c.lower()]
    if bad:
        raise AssertionError(
            "LEAKAGE GUARD TRIPPED. These columns must not be features:\n  "
            + "\n  ".join(sorted(set(bad)))
            + "\nEither remove them or, if you are certain, edit "
              "features._assert_no_leakage and justify it in the write-up."
        )
