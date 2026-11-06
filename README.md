<p align="center">
  <img src="wesbite/logo.png" alt="PlasmidSurfer" width="420">
</p>

<h3 align="center">Find a plasmid's relatives, see where they have been, and judge how worried to be.</h3>

<p align="center">
  <a href="https://plasmidsurfer.vercel.app/"><b>🌐 Open the website</b></a> ·
  <a href="#what-you-can-do-on-the-website">Features</a> ·
  <a href="#how-the-risk-verdict-works">Risk verdict</a> ·
  <a href="#run-it-yourself">Run it yourself</a> ·
  <a href="#team">Team</a>
</p>

---

## What PlasmidSurfer is

PlasmidSurfer is an interactive map of **10,114 *Klebsiella pneumoniae* plasmids from PLSDB**, built for infection-control and surveillance teams.

Each dot is a plasmid. Plasmids are linked when **PLING** finds them evolutionarily related: a short DCJ distance means few structural rearrangements separate them. Related plasmids group into PLING communities, which work as plasmid lineages.

Paste a plasmid sequence and PlasmidSurfer searches **PLSDB**. It lights up the closest matches in the map and shows, for each one:

- **What it carries:** AMR genes, replicon and mobility.
- **Where its family has been:** countries, host species and years, for every member of its PLING community in PLSDB, not only *K. pneumoniae*.
- **How worried to be:** a transparent risk verdict built from resistance, mobility and spread.

### The question it answers

> A carbapenem-resistant *K. pneumoniae* turns up on a ward. You sequence it and paste its plasmid into PlasmidSurfer. In under a minute you know:
>
> - whether it belongs to a known epidemic plasmid family;
> - which countries and host species that family has reached;
> - whether its relatives have been seen carrying something worse.

---

## What you can do on the website

| Feature | What it shows |
|---|---|
| **Relatedness graph** | Every plasmid in the set, laid out by PLING community. Colour by community, AMR tier, risk verdict, country or year. Click a legend entry to isolate that group. |
| **Find** | Jump to any plasmid by accession. |
| **PLSDB search** | Paste or upload a FASTA file and search PLSDB with **Mash dist** (whole plasmids, fast), **Mash screen** (containment) or **BLASTn** (genes or fragments of 100–5,000 bp, up to 10 records). Filter by minimum identity and coverage. |
| **Hits in the graph** | Matching plasmids light up and everything else greys out. Click a hit to see its details and PLING relatives. |
| **Hide hub plasmids** | Drops the 1,149 PLING hub plasmids from the results. Hubs link many communities and can crowd out the family your query actually belongs to. |
| **World map of hits** | One dot per country, sized by the number of hits. Click a country to list its plasmids. |
| **Plasmid panel** | Risk verdict with the resistance, mobility and spread scores. Also AMR genes, replicon, mobility, host, country, year, length and One-Health compartment. |
| **Community epidemiology** | For the plasmid's whole PLING community in PLSDB: number of members, countries, host range and years. |
| **Community map** | *Show on map* (next to Countries): where the whole community has been found, with one dot per country sized by its number of plasmids. |
| **Community timeline** | *Show timeline* (next to Years): a stacked bar chart of the community's plasmids per year, coloured by host species, to show when it crossed into other species. |

<!-- Optional: add screenshots, e.g. ![Search results](docs/img/search.png) -->

---

## How the risk verdict works

The verdict is **rule-based, not a black box**. Every number traces back to a weight table you can argue with. The full write-up is in [`docs/METHODS.md`](docs/METHODS.md).

| Axis | Question | Weight |
|---|---|---|
| **AMR cargo** | What does this plasmid carry right now? Carbapenemases and *mcr* score highest, *sul*/*dfr* lowest. | 0.45 |
| **Mobility** | Can it move? Conjugative > mobilisable > non-mobilisable. | 0.35 |
| **Spread** | How far has its PLING community spread across countries, hosts and One-Health compartments? | 0.20 |

```
overall = 0.45 × AMR + 0.35 × Mobility + 0.20 × Spread     (0–100)
        + 15 if the "community gap" flag fires
```

**Community gap** is a rare-event flag. It fires when the plasmid's own lineage has been seen carrying worse resistance than this plasmid does. It tells you what the backbone tends to acquire, not just what it carries today.

The overall score is reported as a band rather than a precise number:

| Overall | Verdict |
|---|---|
| < 30 | Low concern |
| 30–55 | Monitor |
| 55–75 | Concerning |
| ≥ 75 | High concern |

Spread gets the lowest weight on purpose. It describes the lineage rather than this plasmid, and it is the axis most affected by sampling bias.

---

## Data

| Source | Used for |
|---|---|
| [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) | Plasmid sequences and metadata: host taxonomy, country, dates, One-Health compartment. Also the live sequence search, through the PLSDB API. |
| [PLING](https://github.com/iqbal-lab-org/pling) | Plasmid relatedness: DCJ-indel distances, communities and hub plasmids. |
| MOB-suite | Replicon, relaxase, MPF, oriT and predicted mobility. |
| AMRFinderPlus / Kleborate | AMR genes and resistance tiers. |
| [Natural Earth](https://www.naturalearthdata.com/) via `world-atlas` | World map outline (public domain). |

**Caveats to keep in mind**

- PLSDB reflects what has been sequenced, not everything that circulates. High-income countries and clinical isolates are over-represented.
- PLSDB has no true collection date. "Year" is the year the sequence record was created in NCBI, which lags isolation.
- The verdict is a genotypic indicator for triage. It is not a clinical prediction.

---

## How it is built

```
Browser (static page)                          Vercel
┌────────────────────────────────────┐   ┌────────────────────────┐
│ index.html  sigma.js + graphology  │   │ api/plsdb.js           │
│ data.js     graph, metadata,       │──▶│ same-origin proxy  ───▶│──▶ PLSDB API
│             scores                 │   │ (submit + poll job)    │   (Mash / BLASTn)
│ world.js, community_*.js           │◀──│                        │
└────────────────────────────────────┘   └────────────────────────┘
```

The graph layout and all scores are computed offline in Python, so the browser only draws. Sequence searches go to the PLSDB API through a small serverless proxy, so no BLAST server needs hosting.

### Repository layout

```
wesbite/                      the website (deployed on Vercel)
├── index.html                the viewer
├── data.js                   10,114 plasmids: layout, metadata, scores, community stats
├── world.js                  world outline + country centres for the maps
├── community_countries.js    plasmids per country for every PLING community (all of PLSDB)
├── community_timeline.js     plasmids per year and host species for every community
├── api/plsdb.js              serverless proxy to the PLSDB sequence-search API
├── vercel.json               Vercel routing
├── logo.png, mark.png        logos
└── lib/                      sigma.js, graphology (vendored)

```

---

## Run it yourself

**View locally.** The graph and panels work. The PLSDB search needs the Vercel proxy.

```bash
cd wesbite
python3 -m http.server 8000
# open http://localhost:8000
```

**Deploy.**

1. Import this repository into [Vercel](https://vercel.com).
2. Set the root directory to `wesbite/`.
3. Deploy. `api/plsdb.js` becomes the `/api/plsdb` endpoint automatically.

**Rebuild the data.** The scoring and graph-export pipeline (`score_klebsiella.py`, `export_graph.py`) takes the PLSDB tables, PLING output, MOB-suite and AMRFinder results. See `docs/METHODS.md`. [Add the scripts here, or a link to where they live.]

---

## Team

<!-- Add names, affiliations, roles and links -->

| Name | GitHub user | ORCID ID |
|---|---|---|
| Manuel Ares-Arroyo | [ares-arroyo](https://github.com/ares-arroyo) |  |
| Riccardo Bollini | [riccabolla](https://github.com/riccabolla) | |
| Riccardo Polani | [riccardopolani-max](https://github.com/riccardopolani-max) | |
| Anurag Kumar Bari | [AnuragKumarBari](https://github.com/AnuragKumarBari)||
| Laura Monge-Olivares | [lauramonge9](https://github.com/lauramonge9) | |
| Angel Molina Guillen | [anmogu](https://github.com/anmogu)| |
| Miriam Russell | [miriam-russell](https://github.com/miriam-russell) | |
| Christian Owusu-Nyantaki | [reiychris](https://github.com/Reiychris) | |
| George Kitundu | [eorgekit](https://github.com/EorgeKit) | |

Developed at the **[SPREAD Hackathon, place, date]**.

## Acknowledgments

<!-- Funders, organisers, mentors, data providers -->

- [Hackathon organisers / funders]
- The PLSDB team, for the database and the public sequence-search API
- The PLING developers
- [Mentors, institutions, compute providers]

## How to cite

[Citation or preprint, when available]

## Contact

[Email or GitHub issues link]

## License

See [`LICENSE`](LICENSE).
