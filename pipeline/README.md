# Rebuilding the website data

These scripts regenerate every data file the website in `wesbite/` loads:

| Website file | Made by | Step |
|---|---|---|
| `data.js`: graph, metadata, scores, community stats | `score_klebsiella.py` → `export_graph.py` | 2–3 |
| `community_countries.js`: per-country counts for each community | `build_web_extras.py` | 4 |
| `community_timeline.js`: per-year × species counts for each community | `build_web_extras.py` | 4 |
| `world.js`: map outline and country centres | `make_world_js.js` | 5 (only if you want a different map) |

How the scores are calculated is described in [`../docs/METHODS.md`](../docs/METHODS.md).

```
pipeline/
├── plasmidsurfer/            Python package: scoring rules, metadata cleaning, PLING grouping
├── scripts/
│   ├── score_klebsiella.py   step 2: AMR, Mobility, Spread + community-gap flag -> kp_scores.tsv
│   ├── export_graph.py       step 3: scores + PLING distances -> data.js (layout + panel stats)
│   ├── add_isolates.py       optional: add unlinked plasmids to an existing data.js
│   ├── build_web_extras.py   step 4: community_countries.js + community_timeline.js
│   └── make_world_js.js      step 5: world.js (Node)
└── requirements.txt          pandas, numpy, networkx, openpyxl
```

---

## 0. Install

```bash
cd pipeline
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## 1. Put the inputs in `pipeline/data/`

| File | What it is | Where it comes from |
|---|---|---|
| `nuccore.csv`, `taxonomy.csv`, `biosample.csv` | PLSDB metadata tables | [PLSDB download](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download) |
| `typing.tsv` | PLING community of every plasmid (`plasmid`, `type`) | PLING run on all of PLSDB |
| `hub_plasmids.csv` | PLING hub plasmids | PLING run |
| `all_plasmids_distances.tsv` | PLING pairwise DCJ distances (`plasmid_1`, `plasmid_2`, `distance`; ~420 MB) | PLING run |
| `typing.csv` | MOB-suite typing: mobility, relaxase, MPF, oriT, replicons | MOB-suite on the plasmids |
| `Pling_AMR_Kleborate_Master_10116.xlsx` | AMRFinderPlus hits, sheet `Raw_AMRFinder` | AMRFinderPlus on the plasmids |
| `kp_plasmid_ids.tsv` | the 10,116 *K. pneumoniae* plasmid accessions (first column) | team list |

> The large files (PLSDB tables, PLING output) are not in this repository. Download or regenerate them, or ask the team. [Add a Zenodo / Drive link here.]

## 2. Score every plasmid

```bash
python scripts/score_klebsiella.py \
    --typing          data/typing.csv \
    --amrfinder-xlsx  data/Pling_AMR_Kleborate_Master_10116.xlsx \
    --pling           data/typing.tsv \
    --hubs            data/hub_plasmids.csv \
    --nuccore         data/nuccore.csv \
    --taxonomy        data/taxonomy.csv \
    --biosample       data/biosample.csv \
    --ids             data/kp_plasmid_ids.tsv \
    --require-species "Klebsiella pneumoniae" \
    --min-length      1000 \
    --out             data/kp_scores.tsv
```

- **Output:** `data/kp_scores.tsv`, with 10,114 plasmids.
- **`--require-species`** stops the run if any non-*K. pneumoniae* plasmid slips in.
- **`--min-length`** drops the 2 plasmids under 1 kb.
- **Redundancy report:** the script ends with a Spearman matrix showing how far the scores overlap.

## 3. Build `data.js`

```bash
python scripts/export_graph.py \
    --scores        data/kp_scores.tsv \
    --pling         data/all_plasmids_distances.tsv \
    --pling-typing  data/typing.tsv \
    --nuccore       data/nuccore.csv \
    --taxonomy      data/taxonomy.csv \
    --biosample     data/biosample.csv \
    --out           ../wesbite/data.js
```

- **What it does:** lays out the graph (one blob per PLING community, ForceAtlas2), keeps each plasmid's closest links, and computes the verdict. It also writes the community-epidemiology figures over **all** PLSDB members of each community. These need `--pling-typing` together with `--nuccore`.
- **`--layout-iter`** controls the layout. Use 25 for a quick look and 60 (the default) for the final version.
- **Do not add `--connected-only`.** It drops every plasmid with no PLING link, about 4,000 of them.

## 4. Build the community map and timeline files

```bash
python scripts/build_web_extras.py \
    --pling-typing data/typing.tsv \
    --nuccore      data/nuccore.csv \
    --taxonomy     data/taxonomy.csv \
    --biosample    data/biosample.csv \
    --data-js      ../wesbite/data.js \
    --out-dir      ../wesbite/
```

This writes `community_countries.js` and `community_timeline.js`, one entry per community in `data.js`. It uses the same country and year cleaning as the community panel, so the numbers agree with it.

## 5. (Optional) Rebuild the world map

```bash
npm install world-atlas@2.0.2 topojson-client@3.1.0
node scripts/make_world_js.js > ../wesbite/world.js
```

## 6. Check it

```bash
cd ../wesbite && python3 -m http.server 8000     # open http://localhost:8000
```

The graph, panels, maps and timeline work locally. The PLSDB search needs the Vercel deployment, because it goes through `api/plsdb.js`.
