"""
Grouping plasmids into near-identical families.

WHY THIS MODULE IS THE MOST IMPORTANT ONE IN THE REPO
-----------------------------------------------------
PLSDB is full of near-duplicates: the same epidemic plasmid sequenced in 40
different isolates. If you do a random train/test split, the same plasmid ends
up on both sides and your R^2 goes to 0.95 while the model has learned
nothing. Every published plasmid-ML paper that reports suspiciously good
numbers has made this mistake.

So we need groups, and we use them for TWO things:

  1. Cross-validation: whole groups go to train or to test, never split.
  2. The label itself: "dissemination success" is a property of a plasmid
     FAMILY (how widely has this lineage been observed?), not of one record.

Three possible sources of grouping, in order of preference:

  PLING communities      -- best. Uses DCJ-indel distance, which respects
                            rearrangements. Use as soon as the run finishes.
  Mash similarity graph  -- plsdb_mashdb_sim.tsv gives ~525k similar pairs over
                            ~32k plasmids. Connected components of that graph
                            are a good stand-in and are available RIGHT NOW.
  MOB-suite clusters     -- primary_cluster_id in typing.csv, 100% coverage.
                            Used as a fallback for the ~40k plasmids that have
                            no edge in the Mash graph.

We combine Mash components with MOB fallback so every plasmid gets a group.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


def _union_find(pairs) -> dict:
    """Connected components without pulling in networkx. Path-compressed."""
    parent: dict = {}

    def find(x):
        parent.setdefault(x, x)
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:          # path compression
            parent[x], x = root, parent[x]
        return root

    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return {node: find(node) for node in parent}


def mash_components(mash_path: Path) -> pd.Series:
    """
    Connected components of the Mash similarity graph.

    plsdb_mashdb_sim.tsv is a headerless two-column edge list of plasmid
    accessions that Mash considers similar. Returns Series: ACC -> 'mash_<id>'.
    """
    edges = pd.read_csv(mash_path, sep="\t", header=None,
                        names=["a", "b"], dtype=str)
    # Guard against the file having a header after all.
    if edges.iloc[0].a in {"NUCCORE_ACC", "a", "source"}:
        edges = edges.iloc[1:]
    comp = _union_find(edges.itertuples(index=False, name=None))
    s = pd.Series(comp, name="group")
    # Rename roots to stable short ids
    codes = {root: f"mash_{i}" for i, root in enumerate(sorted(set(s.values)))}
    return s.map(codes)


def pling_communities(pling_path: Path) -> pd.Series | None:
    """
    Read a PLING community table if it exists.

    Accepts any 2+ column delimited file, with or without a header. Column 1 is
    the plasmid id, column 2 is the community id. If your PLING output names
    plasmids by filename (e.g. 'NZ_CP139252.1.fna'), the suffix is stripped.
    """
    if not Path(pling_path).exists():
        return None
    df = pd.read_csv(pling_path, sep=None, engine="python", dtype=str, header=None)
    # Drop a header row if the first row does not look like an accession.
    if not str(df.iloc[0, 0]).replace("_", "").replace(".", "").isalnum() or \
            str(df.iloc[0, 0]).lower() in {"plasmid", "plasmid_id", "nuccore_acc", "id"}:
        df = df.iloc[1:]
    ids = (df.iloc[:, 0].astype(str)
           .str.replace(r"\.(fna|fasta|fa|gbk|gff)(\.gz)?$", "", regex=True))
    comm = "pling_" + df.iloc[:, 1].astype(str)
    s = pd.Series(comm.values, index=ids.values, name="group")
    return s[~s.index.duplicated(keep="first")]


def assign_groups(typing: pd.DataFrame,
                  mash_path: Path,
                  pling_path: Path | None = None,
                  verbose: bool = True) -> pd.DataFrame:
    """
    Give every plasmid in `typing` a group id and record where it came from.

    Returns DataFrame indexed by NUCCORE_ACC with columns [group, group_source].
    """
    accs = typing["NUCCORE_ACC"].astype(str)
    out = pd.DataFrame(index=pd.Index(accs, name="NUCCORE_ACC"))
    out["group"] = pd.NA
    out["group_source"] = pd.NA

    # 1. PLING, if available -- highest priority.
    if pling_path is not None:
        pl = pling_communities(pling_path)
        if pl is not None:
            hit = out.index.intersection(pl.index)
            out.loc[hit, "group"] = pl.reindex(hit).values
            out.loc[hit, "group_source"] = "pling"
            if verbose:
                print(f"  PLING communities    : {len(hit):>7,} plasmids "
                      f"({pl.loc[hit].nunique():,} communities)", file=sys.stderr)

    # 2. Mash components for whatever PLING did not cover.
    if Path(mash_path).exists():
        mc = mash_components(mash_path)
        todo = out.index[out.group.isna()]
        hit = todo.intersection(mc.index)
        out.loc[hit, "group"] = mc.reindex(hit).values
        out.loc[hit, "group_source"] = "mash"
        if verbose:
            print(f"  Mash components      : {len(hit):>7,} plasmids "
                  f"({pd.Series(mc.reindex(hit)).nunique():,} components)",
                  file=sys.stderr)

    # 3. MOB-suite primary cluster as the catch-all.
    mob = ("mob_" + typing["primary_cluster_id"].astype(str)).values
    mob = pd.Series(mob, index=accs.values)
    mob = mob[~mob.index.duplicated(keep="first")]
    todo = out.index[out.group.isna()]
    out.loc[todo, "group"] = mob.reindex(todo).values
    out.loc[todo, "group_source"] = "mob_cluster"
    if verbose:
        print(f"  MOB-suite fallback   : {len(todo):>7,} plasmids "
              f"({pd.Series(mob.reindex(todo)).nunique():,} clusters)",
              file=sys.stderr)

    # Any plasmid still without a group becomes its own singleton.
    orphan = out.group.isna()
    out.loc[orphan, "group"] = ["solo_" + a for a in out.index[orphan]]
    out.loc[orphan, "group_source"] = "singleton"

    if verbose:
        sizes = out.group.value_counts()
        print(f"  TOTAL                : {len(out):>7,} plasmids in "
              f"{out.group.nunique():,} groups "
              f"(median size {sizes.median():.0f}, max {sizes.max():,})",
              file=sys.stderr)
    return out
