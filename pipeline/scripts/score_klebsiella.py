#!/usr/bin/env python3
"""
Three scores for every Klebsiella plasmid: AMR, Mobility, Community.

Rule-based, no machine learning. Produces one table plus a redundancy report
that tells you whether the two community-score variants are actually
measuring different things.

Inputs
------
  --typing       MOB-suite table (required)        -> mobility score
  --pling        PLING community assignments       -> community score
  --kleborate    Kleborate kpsc__amr output        -> AMR score
  --amr          AMRFinder/CARD amr.tsv            -> AMR score (fallback)
  --nuccore      PLSDB master table                -> host genus, compartments
  --mash         Mash similarity edge list         -> community fallback

Usage
-----
    python score_klebsiella.py \
        --typing data/typing.csv \
        --kleborate data/kleborate_amr.tsv \
        --pling data/pling_communities.tsv \
        --nuccore data/nuccore.csv \
        --taxonomy data/taxonomy.csv \
        --biosample data/biosample.csv \
        --out out/klebsiella_scores.tsv

Either --kleborate or --amr is required. If Kleborate is not ready yet, run
with --amr to exercise the whole pipeline on AMRFinder output and swap later;
the AMR score is the only part that changes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError as _e:        # pragma: no cover
    import sys as _s
    _s.exit(
        f"\nERROR: {_e.name} is not installed for this interpreter "
        f"({_s.executable}).\n\n"
        "If you have already run ./setup.sh, activate the environment:\n\n"
        "    source .venv/bin/activate\n\n"
        "If you have not, run ./setup.sh first.\n")

# Make `import plasmidsurfer` work when this script is run directly from
# scripts/, without requiring PYTHONPATH or an install step.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))

from plasmidsurfer import cargo as cargo_mod
from plasmidsurfer import grouping, kleb_scores, metadata as meta_mod
from plasmidsurfer import ranking


# Display label for PLING hub plasmids. They are not a lineage -- see the
# long comment where it is applied.
HUB_LABEL = "hub plasmid"


def log(m: str) -> None:
    print(m, file=sys.stderr, flush=True)


def read_any(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep, low_memory=False)
    if df.shape[1] == 1:
        df = pd.read_csv(path, sep="," if sep == "\t" else "\t", low_memory=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--typing", type=Path, required=True)
    ap.add_argument("--kleborate", type=Path, default=None)
    ap.add_argument("--amr", type=Path, default=None)
    ap.add_argument("--amrfinder-xlsx", type=Path, default=None,
                    help="AMRFinderPlus master workbook. Reads the "
                         "'Raw_AMRFinder' sheet (one row per gene hit) and "
                         "uses it in place of --amr / --kleborate")
    ap.add_argument("--min-length", type=int, default=1000,
                    help="drop plasmids shorter than this many bp "
                         "(default 1000; 0 keeps everything)")
    ap.add_argument("--require-species", default=None,
                    help="abort unless every plasmid has this host species, "
                         'e.g. "Klebsiella pneumoniae"')
    ap.add_argument("--pling", type=Path, default=None,
                    help="PLING typing table (plasmid, type)")
    ap.add_argument("--pling-level", choices=["subcommunity", "community"],
                    default="subcommunity",
                    help="which PLING level is a 'lineage'. subcommunity "
                         "(default) is the right unit; community is far too "
                         "coarse on this dataset -- see grouping.py")
    ap.add_argument("--mash", type=Path, default=None)
    ap.add_argument("--nuccore", type=Path, default=None)
    ap.add_argument("--taxonomy", type=Path, default=None)
    ap.add_argument("--biosample", type=Path, default=None)
    ap.add_argument("--hubs", type=Path, default=None,
                    help="PLING hub_plasmids.csv -- feeds hub connectivity "
                         "into the spread score")
    ap.add_argument("--ids", type=Path, default=None,
                    help="explicit list of plasmids to score (one accession "
                         "per line, or a table whose first column is the "
                         "accession). Use this when your team has a canonical "
                         "set -- it overrides --genus and removes any chance "
                         "of the tool re-deriving a slightly different "
                         "population")
    ap.add_argument("--genus", default="Klebsiella",
                    help="restrict to this host genus (default Klebsiella)")
    ap.add_argument("--out", type=Path, default=Path("out/klebsiella_scores.tsv"))
    args = ap.parse_args()

    if not (args.kleborate or args.amr or args.amrfinder_xlsx):
        raise SystemExit(
            "Provide an AMR source: --amrfinder-xlsx (preferred), "
            "--kleborate, or --amr.")

    # ---- load ------------------------------------------------------------
    log("Loading")
    typing = read_any(args.typing)
    typing["NUCCORE_ACC"] = typing["NUCCORE_ACC"].astype(str)
    log(f"  typing                   : {len(typing):,} plasmids")

    md = None
    if args.nuccore and args.nuccore.exists():
        md = meta_mod.load_metadata(args.nuccore, args.taxonomy, args.biosample)

    # ---- restrict the population ------------------------------------------
    #
    # An explicit id list always wins. Deriving the set from a genus looked
    # convenient but silently disagreed with the team's own list: a genus
    # filter on NCBI taxonomy gives 12,213 Klebsiella plasmids, while their
    # canonical set is 10,116 -- K. pneumoniae only, and only those PLING
    # actually processed. Both are defensible; having two is not.
    if args.ids and args.ids.exists():
        raw = pd.read_csv(args.ids, sep=None, engine="python", dtype=str,
                          header=None, nrows=1)
        first = str(raw.iloc[0, 0]).strip().lower()
        has_header = first in {"plasmid", "plasmid_id", "nuccore_acc", "id",
                               "accession", "hub_plasmids"}
        ids = pd.read_csv(args.ids, sep=None, engine="python", dtype=str,
                          header=0 if has_header else None)
        wanted = set(ids.iloc[:, 0].astype(str).str.strip())
        before = len(typing)
        typing = typing[typing.NUCCORE_ACC.isin(wanted)]
        log(f"  --ids list               : {len(wanted):,} requested, "
            f"{len(typing):,} found in typing ({before:,} available)")
        missing = len(wanted) - len(typing)
        if missing:
            log(f"  !! {missing:,} requested plasmids are not in typing.csv "
                f"and cannot be scored")
    elif md is not None and "host_genus" in md.columns:
        keep = md.index[md["host_genus"].astype(str).str.strip() == args.genus]
        typing = typing[typing.NUCCORE_ACC.isin(set(keep))]
        log(f"  restricted to {args.genus:<12} : {len(typing):,} plasmids")
    else:
        log(f"  !! no metadata and no --ids; scoring everything")

    if not len(typing):
        raise SystemExit(f"No {args.genus} plasmids found.")

    # ---- length filter ----------------------------------------------------
    # Sub-kilobase "plasmids" are mostly assembly fragments and tiny cryptic
    # replicons. They carry no cargo, have no machinery, and their scores are
    # noise that still occupies a dot on the map.
    if args.min_length and "size" in typing.columns:
        sz = pd.to_numeric(typing["size"], errors="coerce")
        before = len(typing)
        typing = typing[sz.fillna(0) >= args.min_length]
        if before - len(typing):
            log(f"  dropped {before - len(typing):,} plasmids shorter than "
                f"{args.min_length:,} bp")

    # ---- species guard ----------------------------------------------------
    if args.require_species and md is not None and "host_species" in md.columns:
        sp = (md["host_species"].reindex(typing.NUCCORE_ACC)
              .astype(str).str.replace("_", " ").str.strip())
        want = args.require_species.strip()
        off = sp[sp.ne(want) & sp.ne("nan")]
        if len(off):
            log(f"  !! {len(off):,} plasmids are not {want}:")
            for k, v in off.value_counts().head(5).items():
                log(f"       {k}  ({v:,})")
            raise SystemExit(
                f"--require-species {want!r} failed. Fix the --ids list, or "
                f"drop the flag to score them anyway.")
        log(f"  species check            : all {len(typing):,} are {want}")

    # ---- AMR score --------------------------------------------------------
    log("\nAMR score (ranked tiers -- see plasmidsurfer/ranking.py)")
    if args.amrfinder_xlsx and args.amrfinder_xlsx.exists():
        # One row per gene hit, which is exactly what the ladder wants. We
        # take the gene SYMBOL rather than the class column: AMRFinder's
        # Class is "BETA-LACTAM" for both blaTEM-1 and blaNDM-5, so scoring
        # off the class alone would put a carbapenemase and a narrow-spectrum
        # penicillinase in the same tier.
        raw = pd.read_excel(args.amrfinder_xlsx, sheet_name="Raw_AMRFinder")
        raw.columns = [str(c).strip() for c in raw.columns]
        sym = next((c for c in raw.columns
                    if c.lower() in ("element symbol", "gene_symbol",
                                     "gene symbol", "symbol")), None)
        if sym is None:
            raise SystemExit(
                f"No gene-symbol column in Raw_AMRFinder. "
                f"Columns: {list(raw.columns)[:12]}")
        typ = next((c for c in raw.columns if c.lower() == "type"), None)
        if typ:
            before = len(raw)
            raw = raw[raw[typ].astype(str).str.upper().eq("AMR")]
            log(f"  kept {len(raw):,} AMR rows of {before:,} "
                f"(dropped STRESS and VIRULENCE)")
        gene_lists = (raw.groupby("plasmid")[sym]
                      .apply(lambda g: ";".join(sorted(set(
                          str(x).strip() for x in g if str(x).strip())))))
        log(f"  source                   : AMRFinderPlus workbook "
            f"({gene_lists.index.nunique():,} plasmids with >=1 hit)")
    elif args.kleborate and args.kleborate.exists():
        kb = read_any(args.kleborate)
        base = kleb_scores.amr_score_from_kleborate(kb)
        gene_lists = base["amr_genes"]
        log(f"  source                   : Kleborate ({len(base):,} rows)")
    else:
        raw = read_any(args.amr)
        ct = cargo_mod.cargo_table(raw)
        gene_lists = ct["amr_gene_list"]
        log(f"  source                   : AMRFinder/CARD fallback "
            f"({len(ct):,} rows)")
        log("  NOTE: swap in --kleborate when ready; only this block changes")

    tiered = [ranking.amr_tier(g) for g in gene_lists]
    amr = pd.DataFrame({
        "amr_tier": [t[0] for t in tiered],
        "amr_severity": [t[1] for t in tiered],
        "n_carbapenemases": [t[2] for t in tiered],
        "amr_genes": gene_lists.values,
    }, index=gene_lists.index)
    amr["worst_class"] = amr["amr_tier"].map(ranking.TIER_LABELS)
    log("  tier distribution:")
    for k, v in amr.amr_tier.value_counts().items():
        log(f"    {ranking.TIER_LABELS[k]:<38} {v:>6,}")

    # ---- mobility score ---------------------------------------------------
    log("\nMobility score (ranked tiers)")
    mob = ranking.mobility_rank_score(typing)
    log("  tier distribution:")
    for k, v in mob.mobility_tier.value_counts().items():
        log(f"    {k:<20} {v:>6,}")
    # Carried through for display in the viewer and the info panel.
    tt = typing.drop_duplicates("NUCCORE_ACC").set_index("NUCCORE_ACC")
    for src, dst in [("size", "length_bp"), ("gc", "gc"),
                     ("num_contigs", "num_contigs")]:
        if src in tt.columns:
            mob[dst] = tt[src].reindex(mob.index)
    log(f"  scored                   : {len(mob):,} plasmids")
    log(f"  mean mobility            : {mob.mobility_score.mean():.3f}")

    # ---- assemble ---------------------------------------------------------
    d = mob.copy()
    for c in ["amr_severity", "amr_tier", "n_carbapenemases",
              "worst_class", "amr_genes"]:
        if c in amr.columns:
            d[c] = amr[c].reindex(d.index)
    d["amr_severity"] = d["amr_severity"].fillna(0.0)
    d["worst_class"] = d["worst_class"].fillna("no AMR detected")
    d["amr_tier"] = d["amr_tier"].fillna("none")
    d["n_carbapenemases"] = d["n_carbapenemases"].fillna(0).astype(int)
    d["amr_genes"] = d["amr_genes"].fillna("")

    if md is not None:
        for c in ["compartments", "country", "year", "host_species"]:
            if c in md.columns:
                d[c] = md[c].reindex(d.index)
    if "compartments" not in d.columns:
        d["compartments"] = ""

    # ---- communities ------------------------------------------------------
    log("\nCommunities")
    if args.pling and args.pling.exists():
        pl = grouping.pling_communities(args.pling, level=args.pling_level)
        d["community"] = pl.reindex(d.index)
        miss = d["community"].isna()
        if miss.any():
            log(f"  !! {int(miss.sum()):,} plasmids absent from the PLING "
                f"table; each becomes its own singleton")
            d.loc[miss, "community"] = ["solo_" + a for a in d.index[miss]]
        log(f"  source                   : PLING ({args.pling_level})")
    elif args.mash and args.mash.exists():
        mc = grouping.mash_components(args.mash)
        d["community"] = mc.reindex(d.index)
        miss = d["community"].isna()
        d.loc[miss, "community"] = ["solo_" + a for a in d.index[miss]]
        log(f"  source                   : Mash components (PLING stand-in)")
    else:
        raise SystemExit("Provide --pling (preferred) or --mash.")

    comm = kleb_scores.community_scores(d)
    d = d.join(comm.drop(columns=["community"]))

    # ---- epidemiological spread ------------------------------------------
    log("\nSpread score (countries + compartments + lineage size + hubs)")
    hub_ids = None
    if args.hubs and args.hubs.exists():
        hb = read_any(args.hubs)
        hub_ids = set(hb.iloc[:, 0].astype(str))
        log(f"  PLING hubs               : {len(hub_ids):,} "
            f"({len(hub_ids & set(d.index)):,} in this genus)")
    else:
        log("  no --hubs given; hub connectivity contributes 0")
    spread = ranking.spread_score(d, hub_ids=hub_ids)

    for c in ["n_countries", "n_compartments", "n_hubs", "spread_score"]:
        if c in d.columns:
            d = d.drop(columns=[c])
    d = d.join(spread.drop(columns=["community", "community_size"]))
    log(f"  spread_score             : mean {d.spread_score.mean():.3f}  "
        f"sd {d.spread_score.std():.3f}")

    # ---- hub plasmids get a display label, not a community ----------------
    #
    # PLING pulls out the plasmids that bridge many communities BEFORE
    # detecting communities, so a hub has no community assignment at all and
    # arrives here as its own singleton. That is correct for SCORING -- a hub
    # has no lineage, so its community score is rightly NaN.
    #
    # But for display they should read as one group, the way PLING presents
    # them. The relabel happens HERE, after community_scores() has already
    # run: doing it earlier would lump all 1,151 hubs into a single
    # pseudo-lineage and compute a "what does my lineage carry" figure over a
    # bag of unrelated bridging plasmids, which would be nonsense.
    if hub_ids:
        is_hub_row = d.index.isin(hub_ids)
        if is_hub_row.any():
            d.loc[is_hub_row, "community"] = HUB_LABEL
            log(f"  relabelled {int(is_hub_row.sum()):,} hub plasmids as "
                f"'{HUB_LABEL}' for display (community scores stay NaN)")

    # ---- percentiles ------------------------------------------------------
    d["AMR_score"] = kleb_scores.to_percentile(d.amr_severity)
    d["Mobility_score"] = kleb_scores.to_percentile(d.mobility_score)
    d["Community_score_A"] = kleb_scores.to_percentile(d.community_A)
    d["Community_score_B"] = kleb_scores.to_percentile(d.community_B)
    d["Ecology_score"] = kleb_scores.to_percentile(d.ecology_score)
    d["Spread_score"] = kleb_scores.to_percentile(d.spread_score)
    # Variant C is a rare-event flag, not a ranking: most plasmids are at 0,
    # so a percentile of it would be meaningless. Report the raw gap and a
    # boolean.
    d["Community_gap"] = d.community_C
    # potential_wave now comes from community_scores(), which requires four
    # conditions rather than a bare gap threshold -- see kleb_scores.py.
    d["potential_wave"] = d["potential_wave"].fillna(False).astype("boolean")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    d.index.name = "plasmid"
    d.to_csv(args.out, sep="\t")
    log(f"\nWrote {args.out}  ({len(d):,} plasmids x {d.shape[1]} columns)")

    # ---- the redundancy report -------------------------------------------
    print("\n" + "=" * 70)
    print("REDUNDANCY REPORT  (Spearman rho)")
    print("=" * 70)
    corr = kleb_scores.compare_variants(d)
    print(corr.round(3).to_string())

    print("\nWhat to look at:")
    a_mob = corr.loc["community_A", "mobility_score"]
    b_amr = corr.loc["community_B", "amr_severity"]
    ab = corr.loc["community_A", "community_B"]
    print(f"  community_A vs mobility_score : {a_mob:+.3f}")
    print("     Variant A includes a neighbour-mobility term. Neighbours share")
    print("     the backbone, so a high value here means A is largely")
    print("     re-reporting the plasmid's own mobility score.")
    print(f"  community_B vs amr_severity   : {b_amr:+.3f}")
    print("     Some correlation is expected and healthy -- neighbours are")
    print("     similar. It should be clearly below 1.")
    print(f"  community_A vs community_B    : {ab:+.3f}")
    print("     If near 1, pick either. If not, they measure different things")
    print("     and the write-up has to say which one you mean.")

    print("\nScore distributions (raw, before percentile):")
    for c in ["amr_severity", "mobility_score", "community_A", "community_B",
              "community_C", "ecology_score"]:
        s = d[c].dropna()
        print(f"  {c:<16} n={len(s):>6,}  mean {s.mean():.3f}  "
              f"sd {s.std():.3f}  min {s.min():.3f}  max {s.max():.3f}")

    # The headline object: empty now, good vehicle, loaded relatives.
    cand = d[(d.potential_wave.fillna(False)) & (d.Mobility_score > 60)]
    print(f"\n'POTENTIAL WAVE' plasmids -- mobile, and their lineage has been seen")
    print(f"carrying something substantially worse than they do: {len(cand):,}")
    if len(cand):
        cols = [c for c in ["AMR_score", "Mobility_score", "Community_gap",
                            "Ecology_score", "community_size", "nb_amr_max",
                            "predicted_mobility", "rep_types"]
                if c in cand.columns]
        print(cand.nlargest(10, "Community_gap")[cols]
              .round(2).to_string(max_colwidth=26))


if __name__ == "__main__":
    main()
