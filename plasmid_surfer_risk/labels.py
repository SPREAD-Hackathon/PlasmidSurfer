"""
The target variable: observed dissemination success of a plasmid lineage.

WHAT WE ARE PREDICTING, AND WHY NOT "RISK" DIRECTLY
----------------------------------------------------
There is no "risk" column anywhere in PLSDB, and inventing one by hand makes
the model circular (train on "has a carbapenemase", predict "has a
carbapenemase"). So we split the problem:

    RISK  =  DISSEMINATION POTENTIAL   x   CARGO SEVERITY
             (learned, from backbone)      (looked up, from AMRFinder/CARD)

This module builds the first term's TARGET. The question it encodes is:

    "How successfully has this plasmid lineage actually spread, judged by the
     breadth of hosts, countries, One-Health compartments and years it has
     been recovered from, and how often it has been independently
     re-observed?"

A plasmid seen once, in one isolate, in one species, in one country scores
low. A lineage recovered from six genera, across a dozen countries, in both
clinical and wastewater samples, over fifteen years, scores high.

The model then learns to predict that from BACKBONE features only, so you can
ask it about a lineage nobody has tracked yet.

LABEL COMPONENTS (all computed per GROUP, never per plasmid)
--------------------------------------------------------------
Always available:
  n_members      distinct PLSDB records in the lineage. Repeated independent
                 recovery is the most direct evidence of success.
  n_host_genera  distinct host genera. Cross-genus spread is the
                 epidemiologically alarming part.
  host_rank      MOB-suite `observed_host_range_ncbi_rank`, ordinal from
                 species (narrow) to multi-phylla (promiscuous).

Added when the nuccore master table is present (see metadata.py):
  n_countries    distinct countries. Geographic spread.
  n_compartments distinct One-Health compartments (human clinical, gut
                 carriage, animal/food, environment, wastewater). THIS is the
                 clinical-community-environment axis.
  year_span      years between first and last record. Persistence over time.

Each component is rank-transformed to [0,1], averaged, then converted to a
percentile. Rank-transforming rather than z-scoring is deliberate: all of
these are heavily right-skewed count variables and only the ordering is
meaningful.

TWO CAVEATS TO PUT ON A SLIDE
-------------------------------
1. SAMPLING BIAS. PLSDB reflects what has been sequenced, not what circulates.
   Clinical Enterobacterales from high-income countries are massively
   over-represented; environmental and African isolates are nearly absent. So
   "dissemination success" partly measures success at being sequenced. This
   does not invalidate the model, it bounds the claim: it predicts spread
   *within the sampled world*.

2. `year` IS NOT A COLLECTION DATE. It derives from NUCCORE_CreateDate, the
   NCBI record creation date, which lags isolation and clusters around
   submission campaigns. Fine as a coarse persistence proxy, not as evidence
   about any particular year.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

# MOB-suite's observed host range rank, as an ordinal scale.
HOST_RANK_ORDER = {
    "species": 0, "genus": 1, "family": 2, "order": 3,
    "class": 4, "phylum": 5, "multi-phylla": 6, "multi-phyla": 6,
    "superkingdom": 7,
}


def _genus_of(name) -> str | float:
    """First token of a binomial is the genus. 'Escherichia coli' -> 'Escherichia'."""
    if not isinstance(name, str) or not name.strip():
        return np.nan
    tok = name.strip().split()[0]
    return tok if tok[:1].isupper() else np.nan


def _count_unique_multi(series: pd.Series) -> int:
    """Distinct values across comma-separated multi-valued cells."""
    seen: set[str] = set()
    for v in series.dropna():
        for part in str(v).split(","):
            p = part.strip()
            if p:
                seen.add(p)
    return len(seen)


def build_labels(typing: pd.DataFrame,
                 groups: pd.DataFrame,
                 metadata: pd.DataFrame | None = None,
                 use_metadata_components: bool = True,
                 min_members_for_metadata: int = 3,
                 verbose: bool = True) -> pd.DataFrame:
    """
    Compute the dissemination label per group, then broadcast to plasmids.

    Parameters
    ----------
    typing   : MOB-suite table (NUCCORE_ACC, mash_neighbor_identification,
               observed_host_range_ncbi_rank)
    groups   : output of grouping.assign_groups, indexed by NUCCORE_ACC
    metadata : optional output of metadata.load_metadata. When given, the
               label gains its geographic, One-Health and temporal components
               and the host genus switches from the MOB-suite proxy to the
               real NCBI taxonomy.

    Returns
    -------
    DataFrame indexed by NUCCORE_ACC with the per-group statistics, the
    0-100 `dissemination` label, `group_size`, and `label_components`
    recording which components were actually used.
    """
    t = typing.copy()
    t["NUCCORE_ACC"] = t["NUCCORE_ACC"].astype(str)
    t = t.drop_duplicates("NUCCORE_ACC").set_index("NUCCORE_ACC")
    t = t.join(groups[["group"]], how="inner")

    # --- host genus: real taxonomy if we have it, MOB-suite proxy if not ---
    proxy = t["mash_neighbor_identification"].map(_genus_of)
    if metadata is not None and "host_genus" in metadata.columns:
        real = metadata["host_genus"].reindex(t.index)
        t["host_genus"] = real.fillna(proxy)
        t["host_genus_source"] = np.where(real.notna(), "ncbi_taxonomy",
                                          "mob_suite_proxy")
        if verbose:
            frac = (real.notna().mean() * 100)
            print(f"  host genus           : {frac:.1f}% real NCBI taxonomy, "
                  f"rest MOB-suite proxy", file=sys.stderr)
    else:
        t["host_genus"] = proxy
        t["host_genus_source"] = "mob_suite_proxy"

    t["host_rank_ord"] = (t["observed_host_range_ncbi_rank"]
                          .astype(str).str.strip().str.lower()
                          .map(HOST_RANK_ORDER))

    # --- attach the metadata columns the label may use --------------------
    for col in ("country", "compartments", "year"):
        if metadata is not None and col in metadata.columns:
            t[col] = metadata[col].reindex(t.index)

    g = t.groupby("group")
    stats = pd.DataFrame({
        "group_size":     g.size(),
        "n_host_genera":  g.host_genus.nunique(),
        "n_host_species": g.mash_neighbor_identification.nunique(),
        "host_rank_max":  g.host_rank_ord.max(),
    })
    # Groups where nobody had a host rank are narrow, not missing; filling
    # with NaN would silently drop them out of the composite.
    stats["host_rank_max"] = stats["host_rank_max"].fillna(0)

    components = ["group_size", "n_host_genera", "host_rank_max"]

    # The metadata-derived statistics are always COMPUTED (they are useful for
    # inspection and plotting) but only enter the composite when asked for.
    meta_components: list[str] = []
    if "country" in t.columns:
        stats["n_countries"] = g.country.nunique()
        meta_components.append("n_countries")
    if "compartments" in t.columns:
        stats["n_compartments"] = g.compartments.apply(_count_unique_multi)
        meta_components.append("n_compartments")
    if "year" in t.columns:
        yr = g.year.agg(["min", "max"])
        stats["year_span"] = (yr["max"] - yr["min"]).fillna(0)
        stats["first_year"] = yr["min"]
        meta_components.append("year_span")

    if use_metadata_components:
        components += meta_components

    # --- composite ---------------------------------------------------------
    #
    # WHY SMALL LINEAGES NEED SPECIAL HANDLING
    # A lineage with two members can only ever show one or two countries, so
    # its n_countries is driven by lineage size rather than by real
    # geographic spread. Feeding that into the composite injects noise exactly
    # where the label is already weakest. For lineages below
    # `min_members_for_metadata`, the composite is therefore built from the
    # always-available components only, and then rescaled onto the same
    # percentile scale as everything else.
    ranked = pd.DataFrame({c: stats[c].rank(pct=True, method="average")
                           for c in components})
    base = [c for c in components if c not in meta_components]

    if use_metadata_components and meta_components:
        small = stats["group_size"] < min_members_for_metadata
        raw = ranked.mean(axis=1)
        raw_small = ranked[base].mean(axis=1)
        stats["dissemination_raw"] = raw.where(~small, raw_small)
        if verbose:
            print(f"  small lineages (<{min_members_for_metadata} members): "
                  f"{int(small.sum()):,} scored on taxonomic components only",
                  file=sys.stderr)
    else:
        stats["dissemination_raw"] = ranked.mean(axis=1)
    # The final 0-100 is the percentile of the composite, so the label is
    # uniform over lineages by construction. That makes a score readable as
    # "this lineage looks like one that has spread more than X% of lineages
    # in PLSDB".
    stats["dissemination"] = stats["dissemination_raw"].rank(pct=True) * 100

    if verbose:
        print(f"  label components     : {', '.join(components)}",
              file=sys.stderr)

    out = t[["group", "host_genus", "host_genus_source"]].join(stats, on="group")
    out.attrs["label_components"] = components
    return out
