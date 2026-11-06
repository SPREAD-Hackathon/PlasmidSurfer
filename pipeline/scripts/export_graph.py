#!/usr/bin/env python3
"""
Export PlasmidSurfer's `data.js` -- the bridge between the scoring pipeline
and the web viewer.

Reads the three-score table produced by `score_klebsiella.py` plus a community
edge list (PLING, or Mash as a stand-in), computes a layout, and writes the
exact `window.PLASMID_GRAPH = {...}` structure `index.html` expects:

    metadata_columns : [str]            # drives the "Colour by" selector
    nodes            : [{id, accession, pling_type, community, hub, x, y, ...}]
    edges            : [[source, target, distance], ...]

Every column listed in `metadata_columns` appears in the viewer's "Colour by"
dropdown, so adding the scores there is what puts the risk map on screen.

Usage
-----
    python export_graph.py \
        --scores out/klebsiella_scores.tsv \
        --pling data/pling_communities.tsv \
        --out ../plasmidsurfer/data.js

    # Mash stand-in while PLING is pending
    python export_graph.py --scores out/klebsiella_scores.tsv \
        --mash data/plsdb_mashdb_sim.tsv --out data.js

Note on size: `data.js` is loaded synchronously by the browser, so keep it
under a few MB. `--max-nodes` drops the least connected plasmids first; isolated
singletons carry no community information and are the right thing to lose.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

# Running this with the system python instead of the project environment is
# the single most common way to get stuck, and the bare ModuleNotFoundError
# does not say what to do about it. Catch it and say.
try:
    import networkx as nx
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as _e:        # pragma: no cover
    import sys as _s
    _s.exit(
        f"\nERROR: {_e.name} is not installed for this interpreter "
        f"({_s.executable}).\n\n"
        "If you have already run ./setup.sh, you just need to activate the\n"
        "environment in this terminal:\n\n"
        "    source .venv/bin/activate\n\n"
        "If you have not, run ./setup.sh first.\n")

# Make `import plasmidsurfer` work when this script is run directly from
# scripts/, without requiring PYTHONPATH or an install step.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def log(m: str) -> None:
    print(m, file=sys.stderr, flush=True)


# Columns offered in the viewer's "Colour by" selector, in order. The scores
# come first because they are the point of the tool.
METADATA_COLUMNS = [
    "verdict", "overall_score",
    "amr_score", "mobility_score", "spread_score",
    "amr_tier", "is_hub", "ecology_score",
    "community_gap", "potential_wave",
    # worst_class is deliberately absent: it holds the same string as
    # amr_tier, and listing both printed every plasmid's AMR class twice.
    "mobility", "replicon", "host",
    "country", "collection_date", "compartments", "amr_genes", "length_bp",
]

# Verdict bands. These MUST match the BANDS table in web/index.html -- the
# viewer recomputes the verdict client-side so it stays correct if someone
# edits a score by hand, and the two must not disagree.
BANDS = [(30, "Low concern"), (55, "Monitor"),
         (75, "Concerning"), (101, "High concern")]


def overall_and_verdict(amr: float, mob: float, eco: float,
                        potential_wave: bool) -> tuple[float, str]:
    """
    Weighted combination of the three axes, escalated by the potential-wave flag.

    WEIGHTS ARE NOT EQUAL, DELIBERATELY. Epidemiological spread is the
    weakest of the three as evidence about a single plasmid: it is a property
    of the whole lineage, every plasmid in a lineage gets the same value, and
    it is the axis most contaminated by sampling bias -- a lineage looks
    "widespread" partly because somebody sequenced it a lot. Letting it carry
    a third of the score made common lineage members look dangerous on the
    strength of their neighbours' popularity.

        AMR cargo   0.45    what this plasmid actually carries
        Mobility    0.35    whether it can actually move
        Spread      0.20    context, not evidence about this plasmid

    The community gap is a rare-event flag rather than an axis, so it ADDS
    rather than averaging: a clean plasmid whose lineage carries a
    carbapenemase must not be diluted back to "normal" by its own low AMR.

    Reporting a BAND rather than only a number is deliberate. "High concern"
    is defensible given these scores rest on a hand-written weight table and
    PLSDB's sampling; "73.4 / 100" implies a precision that is not there.
    """
    pairs = [(v, w) for v, w in ((amr, 0.45), (mob, 0.35), (eco, 0.20))
             if v == v]
    if not pairs:
        return float("nan"), "unknown"
    # Re-normalise so a missing axis does not silently drag the score down.
    score = sum(v * w for v, w in pairs) / sum(w for _, w in pairs)
    if potential_wave:
        score = min(100.0, score + 15.0)
    for lim, label in BANDS:
        if score < lim:
            return score, label
    return score, BANDS[-1][1]


def _val(row, col, default=""):
    v = row.get(col, default)
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return default
    return v


def build_node(acc: str, r, xy: tuple[float, float], is_graph_hub: bool) -> dict:
    """One plasmid as the viewer wants it.

    Pulled out of main()'s loop so that scripts/add_isolates.py can build the
    unconnected plasmids with EXACTLY this code. Two copies of this mapping
    would drift, and the drift would be silent -- a plasmid would score
    differently depending on whether it happened to have an edge.
    """
    val = _val
    x, y = xy
    _amr = float(val(r, "AMR_score", float("nan")) or 0)
    _mob = float(val(r, "Mobility_score", float("nan")) or 0)
    _eco = float(val(r, "Spread_score", val(r, "Ecology_score", float("nan"))) or 0)
    _gun = bool(val(r, "potential_wave", False) is True)
    _overall, _verdict = overall_and_verdict(_amr, _mob, _eco, _gun)
    return {
        "overall_score": round(_overall, 1),
        "verdict": _verdict,
        "id": acc,
        "accession": acc,
        "community": str(r.get("community", "")),
        # Display form: "pling_community_0_subcommunity_7131" is 36 characters
        # of mostly boilerplate and overflows the panel.
        "pling_type": (str(r.get("community", ""))
                       .replace("pling_community_", "c")
                       .replace("_subcommunity_", " · sub ")),
        "hub": is_graph_hub,
        # --- the three scores, plus the community flag ---------------------
        "amr_score": round(float(val(r, "AMR_score", 0)), 1),
        "mobility_score": round(float(val(r, "Mobility_score", 0)), 1),
        "ecology_score": round(float(val(r, "Ecology_score", 0)), 1),
        "spread_score": round(float(val(r, "Spread_score", 0)), 1),
        "amr_tier": str(val(r, "worst_class", "no AMR detected")),
        "is_hub": bool(val(r, "is_hub", 0)),
        "n_carbapenemases": int(val(r, "n_carbapenemases", 0) or 0),
        "community_gap": round(float(val(r, "Community_gap", 0)), 3),
        "potential_wave": bool(val(r, "potential_wave", False) is True),
        # --- supporting detail for the info panel --------------------------
        "worst_class": str(val(r, "worst_class", "none")),
        "amr_genes": str(val(r, "amr_genes", "")),
        "mobility": str(val(r, "predicted_mobility", "")),
        "replicon": str(val(r, "rep_types", "")),
        "host": str(val(r, "host_species", "")).replace("_", " "),
        "country": str(val(r, "country", "")),
        "collection_date": str(val(r, "year", ""))[:4],
        "compartments": str(val(r, "compartments", "")),
        "community_size": int(val(r, "community_size", 1) or 1),
        "length_bp": str(val(r, "length_bp", "")),
        "x": round(float(x), 2),
        "y": round(float(y), 2),
    }


def load_edges(pling: Path | None, mash: Path | None,
               keep: set[str], chunk_rows: int = 2_000_000) -> pd.DataFrame:
    """
    Edge list restricted to the plasmids being exported.

    READ IN CHUNKS. A full PLING distance file is ~420 MB and tens of millions
    of rows; loading it whole costs several GB of RAM and will take down a
    laptop mid-demo. We stream it and keep only rows where BOTH endpoints are
    in `keep`, which for one genus discards the overwhelming majority.

    PLING output has three columns (plasmid_1, plasmid_2, distance) with a
    header. Mash gives only pairs, so every edge is assigned distance 1 --
    which means the viewer's "Max DCJ distance" slider does nothing until real
    PLING distances arrive. Say so in the demo rather than letting someone
    discover a dead control.
    """
    if pling and pling.exists():
        kept, total = [], 0
        for chunk in pd.read_csv(pling, sep="\t", dtype=str,
                                 chunksize=chunk_rows):
            total += len(chunk)
            chunk.columns = [c.strip().lower() for c in chunk.columns]
            cols = list(chunk.columns)
            a_col = next((c for c in cols if c in
                          ("plasmid_1", "plasmid1", "source", "a")), cols[0])
            b_col = next((c for c in cols if c in
                          ("plasmid_2", "plasmid2", "target", "b")), cols[1])
            d_col = next((c for c in cols if c in
                          ("distance", "dist", "dcj")), cols[2] if len(cols) > 2
                         else None)
            sub = chunk[chunk[a_col].isin(keep) & chunk[b_col].isin(keep)]
            if len(sub):
                out = sub[[a_col, b_col]].copy()
                out.columns = ["a", "b"]
                out["dist"] = (pd.to_numeric(sub[d_col], errors="coerce")
                               if d_col else 1.0)
                kept.append(out)
            log(f"    ...{total:,} rows scanned, "
                f"{sum(len(k) for k in kept):,} kept")
        df = (pd.concat(kept, ignore_index=True) if kept
              else pd.DataFrame(columns=["a", "b", "dist"]))
        df["dist"] = df["dist"].fillna(1.0)
        log(f"  edges from PLING     : {len(df):,} of {total:,} scanned "
            f"(real DCJ distances)")
    elif mash and mash.exists():
        df = pd.read_csv(mash, sep="\t", header=None, names=["a", "b"],
                         dtype=str)
        df["dist"] = 1.0
        df = df[df.a.isin(keep) & df.b.isin(keep)]
        log(f"  edges from Mash      : {len(df):,} "
            f"(distance fixed at 1 -- the DCJ slider will be inert)")
    else:
        raise SystemExit("Provide --pling or --mash.")

    df = df[df.a != df.b]
    log(f"  edges within export  : {len(df):,}")
    return df


def prune_edges(G: nx.Graph, max_degree: int) -> nx.Graph:
    """
    Thin out the edge list without breaking any community apart.

    Dense communities are the whole problem: a 300-member Mash component can
    contribute ~45,000 near-identical edges, which bloats `data.js` far beyond
    what a browser will load comfortably and tells the viewer nothing a few
    hundred edges would not.

    Two passes, so connectivity is never sacrificed to the degree cap:
      1. Keep a minimum spanning forest (by distance). Every community stays
         in one piece, which is what the viewer's community colouring and the
         neighbour panel depend on.
      2. Add back the remaining edges shortest-first, skipping any whose
         endpoints are already at `max_degree`.
    """
    before = G.number_of_edges()
    keep = nx.Graph()
    keep.add_nodes_from(G.nodes)

    forest = nx.minimum_spanning_tree(G, weight="dist")
    for a, b, d in forest.edges(data=True):
        keep.add_edge(a, b, **d)

    rest = sorted(((d.get("dist", 1.0), a, b)
                   for a, b, d in G.edges(data=True)
                   if not keep.has_edge(a, b)), key=lambda t: t[0])
    for dist, a, b in rest:
        if keep.degree(a) >= max_degree or keep.degree(b) >= max_degree:
            continue
        keep.add_edge(a, b, dist=dist)

    log(f"  pruned edges         : {before:,} -> {keep.number_of_edges():,} "
        f"(max degree {max_degree}, spanning forest preserved)")
    return keep


def layout(G: nx.Graph, communities: pd.Series, seed: int = 42,
           iterations: int = 60, gap: float = 1.35) -> dict:
    """
    One compact blob per community, blobs packed side by side on a shelf.

    WHY NOT ONE BIG LAYOUT. Running ForceAtlas2 over the whole graph at once
    pulls every community into a single hairball: the layout optimises global
    edge length, so weakly-linked clusters get dragged together and the
    picture stops showing you the thing the tool is about. It is also slow --
    one N-squared-ish layout instead of many tiny ones.

    So each community is laid out ON ITS OWN, scaled to a radius that grows
    with its membership, and then the blobs are shelf-packed left to right,
    largest first. Communities end up visibly separate, their internal shape
    is preserved, and the whole thing runs in a fraction of the time.

    Small communities get a circle rather than a force layout: with 2-6
    members a force layout is noise, and a ring reads as "these are related,
    and there are five of them" at a glance.
    """
    pos: dict = {}
    rng = np.random.default_rng(seed)

    groups: dict[str, list] = {}
    for n in G.nodes:
        groups.setdefault(str(communities.get(n, "unassigned")), []).append(n)
    order = sorted(groups.items(), key=lambda kv: -len(kv[1]))
    log(f"  laying out {len(order):,} communities separately")

    blobs: list[tuple[str, dict, float]] = []
    for cid, members in order:
        sub = G.subgraph(members)
        k = len(members)
        # Radius grows with sqrt(n) so dot density stays roughly constant
        # instead of big communities becoming solid discs.
        r = 3.0 + 7.0 * np.sqrt(k)

        if k == 1:
            local = {members[0]: (0.0, 0.0)}
        elif k <= 6 or sub.number_of_edges() == 0:
            local = {}
            for i, n in enumerate(sorted(members)):
                a = 2 * np.pi * i / k
                local[n] = (r * 0.6 * np.cos(a), r * 0.6 * np.sin(a))
        else:
            try:
                local = nx.forceatlas2_layout(sub, max_iter=iterations,
                                              seed=seed)
            except Exception:
                local = nx.spring_layout(sub, iterations=min(iterations, 50),
                                         seed=seed)
            xy = np.array(list(local.values()), dtype=float)
            span = max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1])) or 1.0
            cx, cy = xy[:, 0].mean(), xy[:, 1].mean()
            local = {n: ((p[0] - cx) / span * 2 * r,
                         (p[1] - cy) / span * 2 * r)
                     for n, p in local.items()}
            # Re-measure: the force layout's extent is not the nominal radius.
            xy = np.array(list(local.values()), dtype=float)
            r = float(max(np.abs(xy).max(), 1.0))
        blobs.append((cid, local, r))

    return pack_blobs(blobs, gap)


def pack_blobs(blobs: list, gap: float = 1.35) -> dict:
    """Shelf-pack per-community blobs into one tidy, roughly square figure.

    `blobs` is [(community_id, {node: (x, y) local}, radius), ...]. Width is
    chosen so the result is roughly square; blobs go left to right, largest
    first, and a shelf is as tall as its tallest blob.

    Separated from layout() so scripts/add_isolates.py can re-pack with the
    unconnected plasmids included. Adding them to an ALREADY packed figure
    makes satellites collide with neighbouring communities and destroys the
    grid -- the packing has to be redone, not patched.
    """
    pos: dict = {}
    total_area = sum((2 * r * gap) ** 2 for _, _, r in blobs)
    width = np.sqrt(total_area) * 1.15

    x_cursor = 0.0
    y_cursor = 0.0
    shelf_h = 0.0
    for cid, local, r in blobs:
        w = 2 * r * gap
        if x_cursor > 0 and x_cursor + w > width:
            x_cursor = 0.0
            y_cursor -= shelf_h
            shelf_h = 0.0
        cx = x_cursor + w / 2
        cy = y_cursor - r * gap
        for n, (px, py) in local.items():
            pos[n] = (cx + px, cy + py)
        x_cursor += w
        shelf_h = max(shelf_h, w)

    # Centre the whole figure on the origin.
    if pos:
        xs = np.array([p[0] for p in pos.values()])
        ys = np.array([p[1] for p in pos.values()])
        mx, my = (xs.min() + xs.max()) / 2, (ys.min() + ys.max()) / 2
        pos = {n: (p[0] - mx, p[1] - my) for n, p in pos.items()}
        # np.ptp(arr), not arr.ptp() -- the method was removed in NumPy 2.0.
        log(f"  canvas extent        : {np.ptp(xs):,.0f} x {np.ptp(ys):,.0f}")
    return pos


def full_community_table(pling: Path, nuccore: Path, taxonomy: Path | None,
                         biosample: Path | None) -> pd.DataFrame | None:
    """
    Metadata for EVERY member of each PLING community, not just the genus
    being exported.

    This matters more than it sounds. A PLING community is a lineage that
    exists across the whole database; restricting its statistics to the
    Klebsiella members hides the single most interesting thing a community
    can tell you -- that the same backbone has also been recovered from
    Escherichia and Citrobacter. Computed on the Klebsiella subset alone,
    "host genera" is always 1 and "crosses host genera" never fires, because
    we filtered the genus before asking the question.

    So community figures are computed over the full community and the panel
    says how many of those members are reachable in the exported graph.
    """
    if not (pling and Path(pling).exists()
            and nuccore and Path(nuccore).exists()):
        return None
    from plasmidsurfer import metadata as meta_mod

    pl = pd.read_csv(pling, sep=None, engine="python", dtype=str)
    cols = [c.strip().lower() for c in pl.columns]
    pl.columns = cols
    pid = next((c for c in cols if c in ("plasmid", "plasmid_id", "id")), cols[0])
    ptype = next((c for c in cols if c in ("type", "community", "subcommunity")),
                 cols[1])
    pl = pl[[pid, ptype]].dropna()
    pl.columns = ["plasmid", "community"]
    pl["plasmid"] = (pl["plasmid"].astype(str)
                     .str.replace(r"\.(fna|fasta|fa|gbk|gff)(\.gz)?$", "",
                                  regex=True))
    pl["community"] = "pling_" + pl["community"].astype(str)

    md = meta_mod.load_metadata(nuccore, taxonomy, biosample, verbose=False)
    out = pl.set_index("plasmid").join(md, how="left")
    log(f"  full community table : {len(out):,} plasmids across "
        f"{out.community.nunique():,} communities "
        f"({out['host_genus'].notna().mean() * 100:.0f}% with a host genus)")
    return out


def community_epidemiology(s: pd.DataFrame, G: nx.Graph,
                           full: pd.DataFrame | None = None) -> dict:
    """
    The per-community figures the info panel shows when a plasmid is clicked.

    Emitted once per community rather than copied onto every node -- with
    ~1,100 communities and ~12,000 plasmids, repeating the country list on
    each node would roughly double the size of data.js for no benefit.

    Fields:
      members        plasmids in this community that are actually IN the
                     exported graph, i.e. how many relatives you can reach
                     from here -- not the raw PLING membership, which may
                     include plasmids dropped by --connected-only.
      countries      how many distinct countries, plus the names
      host_genera    how many distinct host genera, plus the names
      years          first and last year of isolation, and the span
      metadata       percentage of the metadata fields that are populated
                     across the community -- the honesty bar: a community
                     at 30% is one you should not draw conclusions from
    """
    META_FIELDS = ["country", "year", "host_species", "compartments"]
    out = {}
    present = set(G.nodes)

    # Figures come from the full community when we have it, falling back to
    # the exported subset when we do not.
    source = full if full is not None else s
    src_groups = {cid: g for cid, g in
                  source.groupby(source["community"].astype(str))}

    for cid, sub_exported in s.groupby(s["community"].astype(str)):
        # The hub pseudo-community is a bag of plasmids that bridge many
        # lineages, not a lineage. "What does my community carry" is not a
        # meaningful question about it, so it gets no epidemiology panel.
        if cid == "hub plasmid":
            continue
        reachable = [i for i in sub_exported.index if i in present]
        if not reachable:
            continue
        sub = src_groups.get(cid, sub_exported.loc[reachable])

        countries = sorted({str(c).strip() for c in sub.get(
            "country", pd.Series(dtype=str)).dropna() if str(c).strip()})
        # host_species arrives as "Klebsiella_pneumoniae" -- underscores, not
        # spaces. Splitting on whitespace alone returned the whole binomial,
        # so every species counted as its own "genus" and a pure-Klebsiella
        # community reported 8 genera and "crosses host genera". Split on
        # both separators.
        genera = sorted({re.split(r"[ _]", str(h).strip())[0]
                         for h in sub.get("host_species",
                                          pd.Series(dtype=str)).dropna()
                         if str(h).strip() and str(h)[:1].isupper()})
        yrs = pd.to_numeric(sub.get("year", pd.Series(dtype=float)),
                            errors="coerce").dropna()
        yrs = yrs[(yrs > 1950) & (yrs < 2100)]

        filled = total = 0
        for f in META_FIELDS:
            if f in sub.columns:
                col = sub[f].astype(str).str.strip()
                filled += int((col.notna() & (col != "") & (col != "nan")).sum())
                total += len(sub)
        meta_pct = round(100.0 * filled / total, 0) if total else 0.0

        species = sorted({str(h).strip().replace("_", " ")
                          for h in sub.get("host_species",
                                           pd.Series(dtype=str)).dropna()
                          if str(h).strip()})

        out[cid] = {
            "members": int(len(sub)),
            "reachable": len(reachable),
            "n_host_species": len(species),
            "host_species": species[:6],
            "n_countries": len(countries),
            # Cap the displayed list; the count is the honest number.
            "countries": countries[:6],
            "n_host_genera": len(genera),
            "host_genera": genera[:6],
            "year_min": int(yrs.min()) if len(yrs) else None,
            "year_max": int(yrs.max()) if len(yrs) else None,
            "n_years": int(yrs.nunique()) if len(yrs) else 0,
            "metadata_pct": meta_pct,
            "crosses_genera": len(genera) > 1,
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--pling", type=Path, default=None)
    ap.add_argument("--mash", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=Path("data.js"))
    ap.add_argument("--max-nodes", type=int, default=0,
                    help="cap the export, dropping least-connected first "
                         "(0 = no cap)")
    ap.add_argument("--max-degree", type=int, default=12,
                    help="cap each plasmid's edges, keeping the closest "
                         "neighbours and a spanning forest so no community "
                         "fragments (0 = keep every edge)")
    ap.add_argument("--layout-iter", type=int, default=60,
                    help="ForceAtlas2 iterations. Pure-Python FA2 is slow on "
                         "big graphs -- 60 is a reasonable hackathon setting, "
                         "20 is fine for a first look, 200+ only overnight")
    ap.add_argument("--nuccore", type=Path, default=None,
                    help="PLSDB master table. With --pling-typing, community "
                         "statistics are computed over EVERY member of the "
                         "community, not just the exported genus -- which is "
                         "what makes cross-genus spread visible")
    ap.add_argument("--pling-typing", type=Path, default=None,
                    help="PLING typing table (plasmid, type) covering all "
                         "genera")
    ap.add_argument("--taxonomy", type=Path, default=None)
    ap.add_argument("--biosample", type=Path, default=None)
    ap.add_argument("--connected-only", action="store_true",
                    help="drop plasmids with no edges at all")
    args = ap.parse_args()

    log(f"Reading {args.scores}")
    s = pd.read_csv(args.scores, sep="\t", index_col=0, low_memory=False)
    s.index = s.index.astype(str)
    log(f"  plasmids             : {len(s):,}")

    edges = load_edges(args.pling, args.mash, set(s.index))

    G = nx.Graph()
    G.add_nodes_from(s.index)
    for a, b, d in edges.itertuples(index=False):
        G.add_edge(a, b, dist=float(d))

    # ---- optional pruning -------------------------------------------------
    if args.connected_only:
        iso = [n for n in G if G.degree(n) == 0]
        G.remove_nodes_from(iso)
        s = s.loc[list(G.nodes)]
        log(f"  dropped {len(iso):,} isolated plasmids")
    if args.max_nodes and G.number_of_nodes() > args.max_nodes:
        order = sorted(G.nodes, key=lambda n: G.degree(n), reverse=True)
        keep = set(order[:args.max_nodes])
        G = G.subgraph(keep).copy()
        s = s.loc[list(G.nodes)]
        log(f"  capped at {args.max_nodes:,} nodes (kept the best connected)")

    if args.max_degree:
        G = prune_edges(G, args.max_degree)

    comm = s["community"].astype(str)
    pos = layout(G, comm, iterations=args.layout_iter)

    # ---- hubs: the most connected plasmid in each community ---------------
    deg = pd.Series({n: G.degree(n) for n in G.nodes})
    hubs = set(deg.groupby(comm.reindex(deg.index)).idxmax().dropna())
    log(f"  hubs                 : {len(hubs):,} (one per community)")

    # ---- node records -----------------------------------------------------
    nodes = []
    for acc, r in s.iterrows():
        nodes.append(build_node(acc, r, pos.get(acc, (0.0, 0.0)), acc in hubs))

    edge_list = [[a, b, float(G.edges[a, b]["dist"])] for a, b in G.edges]

    full = full_community_table(args.pling_typing, args.nuccore,
                                args.taxonomy, args.biosample)
    if full is None:
        log("  (no --nuccore/--pling-typing: community figures cover only "
            "the exported genus, so host range will read as 1 genus)")
    epi = community_epidemiology(s, G, full)
    log(f"  community epidemiology: {len(epi):,} communities summarised")

    payload = {"metadata_columns": METADATA_COLUMNS,
               "nodes": nodes, "edges": edge_list, "communities": epi}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("window.PLASMID_GRAPH = ")
        json.dump(payload, fh, separators=(",", ":"))
        fh.write(";\n")

    mb = args.out.stat().st_size / 1e6
    log(f"\nWrote {args.out}  ({len(nodes):,} nodes, {len(edge_list):,} edges, "
        f"{mb:.1f} MB)")
    if mb > 8:
        log("  !! over 8 MB -- data.js loads synchronously and will stall the")
        log("     page on a conference network. Use --max-nodes or "
            "--connected-only.")

    log("\nScores in the export:")
    for c in ["amr_score", "mobility_score", "ecology_score"]:
        v = [n[c] for n in nodes]
        log(f"  {c:<16} mean {np.mean(v):5.1f}  max {np.max(v):5.1f}")
    log(f"  potential_wave flags : {sum(1 for n in nodes if n['potential_wave']):,}")


if __name__ == "__main__":
    main()
