#!/usr/bin/env python3
"""
Put the unconnected plasmids back into an existing web/data.js.

WHY THIS EXISTS
---------------
`export_graph.py --connected-only` drops every plasmid with no PLING edge.
That makes a tidier picture and a smaller file, and it quietly hides a third
of the dataset: 10,114 plasmids scored, 5,970 drawn. If someone searches for
an accession and the viewer says it does not exist, the tool has lied to them.

THE RIGHT FIX is to re-run export_graph.py WITHOUT --connected-only, which
needs the PLING distance file. Use this script only when you have the scores
table and an already-built data.js but not the 421 MB distance file.

It adds one node per missing plasmid, built with export_graph.build_node() --
the same function the connected nodes go through, so the scores, verdict and
panel fields are computed identically. Isolated plasmids have no edges by
definition, so no edge data is needed or invented.

LAYOUT
------
The whole figure is RE-PACKED, with export_graph.pack_blobs(), the same
function that laid it out in the first place. An earlier version of this
script bolted the new plasmids onto the finished picture as satellite discs;
that collided with neighbouring communities and destroyed the tidy grid of
separated clusters. Blobs have to be re-measured and re-packed, not patched.

Inside a community, the connected members keep the exact relative positions
ForceAtlas2 gave them -- that shape is real information and is preserved.
Isolates fill a ring just outside the members. A community whose members are
ALL isolates gets an evenly spaced circle, the same treatment layout() gives
any edgeless community.

USAGE
-----
    python3 scripts/add_isolates.py \\
        --scores data/kp_scores.tsv \\
        --data   web/data.js
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np                                             # noqa: E402
import pandas as pd                                            # noqa: E402

from export_graph import build_node, pack_blobs                # noqa: E402


def log(m: str) -> None:
    print(m, file=sys.stderr, flush=True)


def read_data_js(path: Path) -> tuple[str, dict, str]:
    """Return (prefix, parsed object, suffix) so we can rewrite it in place."""
    raw = path.read_text()
    i = raw.index("{")
    j = raw.rindex("}") + 1
    return raw[:i], json.loads(raw[i:j]), raw[j:]


def ring(k: int, n: int, r: float) -> tuple[float, float]:
    """k-th of n points on a circle of radius r."""
    a = 2 * math.pi * k / max(n, 1)
    return r * math.cos(a), r * math.sin(a)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--data", type=Path, default=Path("web/data.js"))
    ap.add_argument("--gap", type=float, default=1.35,
                    help="spacing between packed communities (as in layout())")
    args = ap.parse_args()

    prefix, obj, suffix = read_data_js(args.data)
    nodes = obj["nodes"]
    by_id = {n["id"]: n for n in nodes}
    log(f"  already in data.js   : {len(by_id):,}")

    s = pd.read_csv(args.scores, sep="\t", index_col=0, low_memory=False)
    s.index = s.index.astype(str)
    missing = [a for a in s.index if a not in by_id]
    log(f"  in scores table      : {len(s):,}")
    log(f"  missing from graph   : {len(missing):,}")
    if not missing:
        log("Nothing to add.")
        return

    # --- group every plasmid, drawn or not, by community -------------------
    comm_of = {a: str(s.at[a, "community"]) for a in s.index} \
        if "community" in s.columns else {a: "" for a in s.index}
    groups: dict[str, list[str]] = {}
    for a in s.index:
        groups.setdefault(comm_of[a], []).append(a)

    # --- rebuild each blob in LOCAL coordinates ----------------------------
    blobs = []
    kept_shape = added_ring = ringed_whole = 0
    for cid, members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        drawn = [a for a in members if a in by_id]
        iso = [a for a in members if a not in by_id]

        if not drawn:
            # No member was ever laid out: evenly spaced circle, exactly what
            # layout() does for a community with no edges.
            k = len(iso)
            r = 3.0 + 7.0 * math.sqrt(k)
            local = {a: ring(i, k, r * 0.6) for i, a in enumerate(sorted(iso))}
            if k == 1:
                local = {iso[0]: (0.0, 0.0)}
            ringed_whole += k
            blobs.append((cid, local, max(r * 0.6, 1.0)))
            continue

        # Keep the force layout's shape: recentre the drawn members on their
        # own centroid and reuse those offsets unchanged.
        cx = sum(by_id[a]["x"] for a in drawn) / len(drawn)
        cy = sum(by_id[a]["y"] for a in drawn) / len(drawn)
        local = {a: (by_id[a]["x"] - cx, by_id[a]["y"] - cy) for a in drawn}
        kept_shape += len(drawn)
        r_in = max((math.hypot(*p) for p in local.values()), default=1.0)

        if iso:
            # Isolates fill an annulus around the members. Concentric rings of
            # a FIXED count look like a sunburst -- every ring lines its dots
            # up on the same angles and the spokes read as structure. A golden
            # -angle spiral with radius ~ sqrt(k) gives constant density and no
            # spokes, and the annulus is sized so that density matches the
            # blob's interior rather than being chosen by eye.
            n = len(iso)
            r0 = r_in * 1.10 + 4.0
            per_node = (math.pi * r_in ** 2) / max(len(drawn), 1)
            r1 = math.sqrt(r0 ** 2 + n * per_node / math.pi)
            ga = math.pi * (3.0 - math.sqrt(5.0))
            for i, a in enumerate(sorted(iso)):
                rr = math.sqrt(r0 ** 2 + (i + 0.5) / n * (r1 ** 2 - r0 ** 2))
                local[a] = (rr * math.cos(i * ga), rr * math.sin(i * ga))
            added_ring += n

        r = max((math.hypot(*p) for p in local.values()), default=1.0)
        blobs.append((cid, local, float(r)))

    log(f"  kept force-layout    : {kept_shape:,} plasmids")
    log(f"  ringed around members: {added_ring:,}")
    log(f"  edgeless communities : {ringed_whole:,} plasmids")
    log(f"  re-packing {len(blobs):,} communities")

    pos = pack_blobs(blobs, args.gap)

    # --- write every node at its new position ------------------------------
    out = []
    for acc in s.index:
        x, y = pos.get(acc, (0.0, 0.0))
        if acc in by_id:
            n = dict(by_id[acc])
            n["x"], n["y"] = round(float(x), 2), round(float(y), 2)
        else:
            n = build_node(acc, s.loc[acc], (x, y), False)
        out.append(n)

    obj["nodes"] = out
    args.data.write_text(prefix + json.dumps(obj, separators=(",", ":")) + suffix)
    log(f"\nwrote {args.data}  ({len(out):,} nodes, "
        f"{len(obj.get('edges', [])):,} edges, "
        f"{args.data.stat().st_size / 1e6:.1f} MB)")
    if len(out) != len(s):
        log(f"  !! {len(s) - len(out):,} unaccounted for")


if __name__ == "__main__":
    main()
