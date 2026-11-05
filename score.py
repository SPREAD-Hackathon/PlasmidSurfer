#!/usr/bin/env python3
"""
STEP 3 -- Score a new plasmid.

INPUT: the same two things PLSDB itself was annotated with, run on your own
assembly:

    mob_typer --infile my_plasmid.fasta --out_file my_typing.txt
    amrfinder -n my_plasmid.fasta --plus -o my_amr.tsv

(the MOB-suite output is tab-separated; this script accepts .csv or .tsv)

OUTPUT: a 0-100 risk score, broken into its two components, plus the nearest
known lineages so a human can sanity-check the number before believing it.

    RISK = 100 x dissemination x (CARGO_FLOOR + (1-CARGO_FLOOR) x cargo)

WHY THE SCORE IS ALWAYS REPORTED AS THREE NUMBERS, NOT ONE
-----------------------------------------------------------
A bare "87/100" invites exactly the kind of false confidence that gets
surveillance tools quietly distrusted by the people who have to act on them.
Reporting dissemination and cargo separately means a reader can immediately
see WHY: a 60 from "very mobile backbone, trivial cargo" is a completely
different object from a 60 from "immobile plasmid carrying NDM-1", and they
warrant different responses.

Usage
-----
    python score.py --typing my_typing.txt --amr my_amr.tsv
    python score.py --typing my_typing.txt --amr my_amr.tsv --json results.json
    python score.py --from-plsdb NZ_CP139252.1     # score an existing record
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from plrisk import cargo as cargo_mod
from plrisk import config, features


def log(m: str) -> None:
    print(m, file=sys.stderr, flush=True)


def read_any(path: Path) -> pd.DataFrame:
    """mob_typer writes TSV; PLSDB dumps are CSV. Sniff it."""
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep, low_memory=False)
    if df.shape[1] == 1:                     # wrong guess, try the other one
        df = pd.read_csv(path, sep="," if sep == "\t" else "\t", low_memory=False)
    return df


def normalise_typing(df: pd.DataFrame) -> pd.DataFrame:
    """
    mob_typer's standalone output uses `sample_id`, the PLSDB dump uses
    `NUCCORE_ACC`. Make them the same, and fill in any column the feature
    builder expects but this run did not produce.
    """
    df = df.copy()
    if "NUCCORE_ACC" not in df.columns:
        for alt in ("sample_id", "id", "filename", "seq_id"):
            if alt in df.columns:
                df["NUCCORE_ACC"] = df[alt].astype(str)
                break
        else:
            df["NUCCORE_ACC"] = [f"query_{i}" for i in range(len(df))]

    needed = ["num_contigs", "size", "gc", "rep_type(s)", "relaxase_type(s)",
              "mpf_type", "orit_type(s)", "predicted_mobility", "PMLST_scheme",
              "primary_cluster_id", "mash_nearest_neighbor",
              "mash_neighbor_identification", "observed_host_range_ncbi_rank"]
    for c in needed:
        if c not in df.columns:
            df[c] = np.nan
    return df


def align_features(F: pd.DataFrame, train_cols: list[str]) -> pd.DataFrame:
    """
    A single query plasmid will not trigger most of the 374 one-hot columns,
    and may trigger a replicon the training set never saw. Reindex to the
    exact training column order: missing -> 0, novel -> dropped.
    """
    novel = [c for c in F.columns if c not in train_cols]
    if novel:
        log(f"  note: {len(novel)} feature(s) in the query were not in the "
            f"training set and are ignored: {', '.join(novel[:6])}"
            f"{' ...' if len(novel) > 6 else ''}")
    return F.reindex(columns=train_cols, fill_value=0.0).astype(float)


def nearest_lineages(typing_row: pd.Series, dataset: pd.DataFrame,
                     k: int = 5) -> pd.DataFrame:
    """
    Context for the score: which known lineages does this plasmid resemble?

    Uses MOB-suite's own mash nearest neighbour when present, otherwise falls
    back to matching the replicon profile. This is deliberately simple -- it
    exists so a microbiologist can say "ah, that is an IncC lineage" and
    decide whether to trust the number, not to be a classifier.
    """
    cols = [c for c in ["group", "group_size", "n_host_genera", "n_countries",
                        "n_compartments", "year_span", "dissemination",
                        "cargo_severity", "host_genus"]
            if c in dataset.columns]

    nn = typing_row.get("mash_nearest_neighbor")
    if isinstance(nn, str) and nn.strip():
        hit = dataset[dataset.index.str.contains(nn.split(".")[0], regex=False,
                                                 na=False)]
        if len(hit):
            return hit.nlargest(k, "group_size")[cols]

    reps = typing_row.get("rep_type(s)")
    if reps is None or (isinstance(reps, float) and pd.isna(reps)):
        reps = typing_row.get("rep_type_profile")
    if isinstance(reps, str) and reps and "rep_type_profile" in dataset.columns:
        same = dataset[dataset.rep_type_profile == reps]
        if len(same):
            return same.nlargest(k, "group_size")[cols]
    return pd.DataFrame()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--typing", type=Path, help="mob_typer output for the query")
    ap.add_argument("--amr", type=Path, help="AMRFinderPlus output for the query")
    ap.add_argument("--from-plsdb", nargs="+", metavar="ACC",
                    help="score accessions already in the built dataset")
    ap.add_argument("--model", type=Path, default=config.OUT_DIR / "model.pkl")
    ap.add_argument("--dataset", type=Path,
                    default=config.OUT_DIR / "dataset.csv.gz")
    ap.add_argument("--json", type=Path, help="also write results as JSON")
    args = ap.parse_args()

    with open(args.model, "rb") as fh:
        bundle = pickle.load(fh)
    model, train_cols = bundle["model"], bundle["features"]
    floor = bundle.get("cargo_floor", config.CARGO_FLOOR)

    dataset = pd.read_csv(args.dataset, index_col=0, low_memory=False)

    # ---------------- path A: score records already in PLSDB --------------
    if args.from_plsdb:
        miss = [a for a in args.from_plsdb if a not in dataset.index]
        if miss:
            log(f"  not found in dataset: {', '.join(miss)}")
        acc = [a for a in args.from_plsdb if a in dataset.index]
        if not acc:
            raise SystemExit("None of the requested accessions are in the dataset.")
        X = dataset.loc[acc, train_cols].astype(float)
        cargo_sev = dataset.loc[acc, "cargo_severity"].astype(float)
        genes = dataset.loc[acc, "amr_gene_list"].fillna("")
        typing_rows = {a: dataset.loc[a] for a in acc}

    # ---------------- path B: score a fresh assembly ----------------------
    else:
        if not args.typing:
            raise SystemExit("Provide --typing (and ideally --amr), or --from-plsdb.")
        typ = normalise_typing(read_any(args.typing))

        amr = read_any(args.amr) if args.amr and args.amr.exists() else None
        if amr is not None and "NUCCORE_ACC" not in amr.columns:
            for alt in ("Name", "Protein identifier", "Contig id", "sample_id"):
                if alt in amr.columns:
                    amr["NUCCORE_ACC"] = amr[alt].astype(str)
                    break
        if amr is not None and "gene_symbol" not in amr.columns:
            for alt in ("Element symbol", "Gene symbol", "symbol"):
                if alt in amr.columns:
                    amr = amr.rename(columns={alt: "gene_symbol"})
                    break
        if amr is not None and "drug_class" not in amr.columns:
            for alt in ("Class", "Subclass", "class"):
                if alt in amr.columns:
                    amr = amr.rename(columns={alt: "drug_class"})
                    break

        # One query plasmid cannot support frequency filtering, so min_count=1
        # and then we align to the training columns.
        F = features.build_features(typ, plasmidfinder=None, cargo=None,
                                    include_stress=False, min_count=1,
                                    verbose=False)
        X = align_features(F, train_cols)

        if amr is not None and len(amr):
            # Every query hit belongs to the query plasmid; if the accession
            # column did not line up, attribute all hits to the single query.
            if X.shape[0] == 1:
                amr = amr.assign(NUCCORE_ACC=X.index[0])
            ct = cargo_mod.cargo_table(amr)
            cargo_sev = ct["cargo_severity"].reindex(X.index).fillna(0.0)
            genes = ct["amr_gene_list"].reindex(X.index).fillna("")
        else:
            log("  no AMR file given -- cargo severity set to 0; the reported "
                "risk is the EMPTY-VEHICLE floor only")
            cargo_sev = pd.Series(0.0, index=X.index)
            genes = pd.Series("", index=X.index)

        typing_rows = {r["NUCCORE_ACC"]: r for _, r in typ.iterrows()}

    # ---------------- combine ---------------------------------------------
    diss = np.clip(model.predict(X), 0, 100)
    cargo_term = floor + (1 - floor) * cargo_sev.to_numpy()
    risk = diss * cargo_term

    results = []
    print()
    for i, acc in enumerate(X.index):
        row = {
            "plasmid": str(acc),
            "risk_score": round(float(risk[i]), 1),
            "dissemination_potential": round(float(diss[i]), 1),
            "cargo_severity": round(float(cargo_sev.iloc[i]), 3),
            "amr_genes": genes.iloc[i],
        }
        results.append(row)

        print("=" * 64)
        print(f"  {acc}")
        print("=" * 64)
        print(f"  RISK SCORE                {row['risk_score']:>6.1f} / 100")
        print(f"    dissemination potential {row['dissemination_potential']:>6.1f} / 100"
              f"   (predicted from backbone)")
        print(f"    cargo severity          {row['cargo_severity']:>6.3f}"
              f"         (looked up from AMR genes)")
        if row["amr_genes"]:
            g = row["amr_genes"].split(";")
            print(f"    AMR genes ({len(g)})          "
                  f"{', '.join(g[:8])}{' ...' if len(g) > 8 else ''}")
        else:
            print("    AMR genes               none detected")

        tr = typing_rows.get(acc)
        if tr is not None:
            print(f"    backbone                size "
                  f"{pd.to_numeric(tr.get('size'), errors='coerce') or float('nan'):,.0f} bp"
                  if pd.notna(tr.get("size")) else "    backbone                size n/a")
            print(f"    mobility                {tr.get('predicted_mobility', 'n/a')}")
            # mob_typer calls it 'rep_type(s)'; the built dataset stores the
            # same values under 'rep_type_profile'. Accept either.
            reps = tr.get("rep_type(s)")
            if reps is None or (isinstance(reps, float) and pd.isna(reps)):
                reps = tr.get("rep_type_profile", "n/a")
            print(f"    replicon(s)             {reps}")
            nb = nearest_lineages(tr, dataset)
            if len(nb):
                print("\n    nearest known lineages in PLSDB:")
                print(nb.to_string().replace("\n", "\n      "))
        print()

    print("  Interpretation: dissemination potential is a percentile against "
          "\n  PLSDB lineages -- 80 means this backbone resembles lineages that "
          "\n  have spread more widely than 80% of those in the database. It is "
          "\n  NOT a probability, and it carries the sampling bias of PLSDB.\n")

    if args.json:
        args.json.write_text(json.dumps(results, indent=2))
        log(f"Wrote {args.json}")


if __name__ == "__main__":
    main()
