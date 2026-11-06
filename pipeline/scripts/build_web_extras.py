#!/usr/bin/env python3
"""
Build the two community files the website's maps and timelines read:

    community_countries.js   plasmids per country, for EVERY member of each
                             PLING community in PLSDB (not only K. pneumoniae)
    community_timeline.js    plasmids per year and host species, same members

Both are computed over the WHOLE community in PLSDB, so the "Show on map" and
"Show timeline" buttons can show when and where a lineage reached other genera.
Only communities that appear in the website's data.js are written.

Country and year come from plasmidsurfer.metadata.load_metadata(), the same
cleaning the scores and the community-epidemiology panel use, so the numbers
agree with the panel. "Year" is the year of the NCBI nucleotide record
(NUCCORE_CreateDate): PLSDB has no true collection date.

Usage
-----
    python scripts/build_web_extras.py \
        --pling-typing data/typing.tsv \
        --nuccore      data/nuccore.csv \
        --taxonomy     data/taxonomy.csv \
        --biosample    data/biosample.csv \
        --data-js      wesbite/data.js \
        --out-dir      wesbite/
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from plasmidsurfer import metadata as meta_mod  # noqa: E402

# Same merges as the community panel / world map.
COUNTRY_NORM = {"Hong Kong": "China", "US": "USA",
                "Czechia": "Czech Republic", "The Gambia": "Gambia"}
TOP_SPECIES = 7          # the rest fold into "Other species" (one grey bar)


def read_communities(data_js: Path) -> set[str]:
    s = data_js.read_text()
    d = json.loads(s[s.index("{"):s.rindex("}") + 1])
    return set(d.get("communities", {}))


def read_typing(path: Path) -> dict[str, str]:
    """PLING typing.tsv (plasmid <tab> type) -> {plasmid: 'pling_' + type}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            p, _, t = line.rstrip("\n").partition("\t")
            if p and p != "plasmid" and t:
                out[p] = "pling_" + t
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pling-typing", type=Path, required=True)
    ap.add_argument("--nuccore", type=Path, required=True)
    ap.add_argument("--taxonomy", type=Path, required=True)
    ap.add_argument("--biosample", type=Path, required=True)
    ap.add_argument("--data-js", type=Path, required=True,
                    help="the website's data.js; only its communities are written")
    ap.add_argument("--out-dir", type=Path, default=Path("."))
    a = ap.parse_args()

    want = read_communities(a.data_js)
    typ = read_typing(a.pling_typing)
    md = meta_mod.load_metadata(a.nuccore, a.taxonomy, a.biosample, verbose=False)

    countries: dict[str, Counter] = defaultdict(Counter)
    years: dict[str, dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    no_year: Counter = Counter()
    missing = 0
    for p, comm in typ.items():
        if comm not in want:
            continue
        if p not in md.index:
            missing += 1
            continue
        r = md.loc[p]
        c = r.get("country")
        c = "" if c is None or c != c else str(c).strip()
        countries[comm][COUNTRY_NORM.get(c, c)] += 1

        sp = r.get("host_species")
        sp = "Unknown species" if sp is None or sp != sp or not str(sp).strip() else str(sp)
        y = r.get("year")
        if y is None or y != y:
            no_year[comm] += 1
        else:
            years[comm][sp][int(y)] += 1

    # ---- community_countries.js -------------------------------------------
    cc = {k: dict(v.most_common()) for k, v in countries.items()}
    (a.out_dir / "community_countries.js").write_text(
        "// Plasmids per country for EVERY member of each PLING community (all of PLSDB, not only\n"
        "// the K. pneumoniae plasmids drawn in the graph). The empty key holds plasmids with no country.\n"
        "// Built by scripts/build_web_extras.py; Hong Kong is counted under China, as in data.js.\n"
        "window.COMMUNITY_COUNTRIES = " + json.dumps(cc, separators=(",", ":")) + ";\n")

    # ---- community_timeline.js --------------------------------------------
    tl = {}
    for comm, by_sp in years.items():
        ys = sorted({y for cnt in by_sp.values() for y in cnt})
        y0, y1 = ys[0], ys[-1]
        totals = Counter({s: sum(c.values()) for s, c in by_sp.items()})
        top = [s for s, _ in totals.most_common(TOP_SPECIES)]
        series = {s: [by_sp[s].get(y, 0) for y in range(y0, y1 + 1)] for s in top}
        rest = [sum(by_sp[s].get(y, 0) for s in by_sp if s not in top)
                for y in range(y0, y1 + 1)]
        if any(rest):
            series["Other species"] = rest
        tl[comm] = {"y0": y0, "n_species": len(by_sp), "series": series,
                    "no_year": no_year[comm]}
    (a.out_dir / "community_timeline.js").write_text(
        "// Plasmids per year and host species for EVERY member of each PLING community (all of PLSDB).\n"
        "// Year = year of the NCBI nucleotide record (NUCCORE_CreateDate), the same year data.js shows.\n"
        f"// series: top {TOP_SPECIES} species (by total) + 'Other species'; arrays run from y0, one value per year.\n"
        "window.COMMUNITY_TIMELINE = " + json.dumps(tl, separators=(",", ":")) + ";\n")

    print(f"communities in data.js : {len(want):,}", file=sys.stderr)
    print(f"written                : {len(cc):,} (countries), {len(tl):,} (timeline)", file=sys.stderr)
    if missing:
        print(f"!! {missing:,} typed plasmids absent from nuccore.csv", file=sys.stderr)


if __name__ == "__main__":
    main()
