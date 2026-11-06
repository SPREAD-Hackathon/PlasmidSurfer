# How PlasmidSurfer scores a plasmid

This page describes exactly what the website shows: the code in `plasmidsurfer/ranking.py`, `plasmidsurfer/kleb_scores.py`, `scripts/score_klebsiella.py` and `scripts/export_graph.py`, and the viewer in `wesbite/index.html`.

The scoring is **rule-based**. There is no machine learning in it. Every number traces back to one line in a table, so anyone can disagree with a specific weight. The separate machine-learning dissemination model is described in [`RISK_MODEL.md`](RISK_MODEL.md); it does **not** feed the website.

---

## 1. Population

- **10,114 *Klebsiella pneumoniae* plasmids** from PLSDB.
- They come from the team's list of 10,116 plasmids that PLING processed: 8,965 typed plasmids plus 1,151 hub plasmids. Two plasmids under 1 kb were dropped.
- **Lineage** means the PLING **subcommunity** (for example `community_0_subcommunity_7133`).
- **Hub plasmids** are the plasmids PLING sets aside because they bridge many communities. They have no lineage of their own and are shown as one "hub plasmid" group.

## 2. The three scores

Each score is first computed on a 0–1 scale. It is then converted to a **percentile within the 10,114 plasmids** (0–100). So "AMR 21" means this plasmid's AMR score is higher than about 21 % of *K. pneumoniae* plasmids in the set. It is not "21 % dangerous".

### 2.1 AMR cargo: what does it carry now?

- **Input:** AMRFinderPlus gene symbols, keeping AMR hits only. Stress genes (metals, biocides) and virulence genes are ignored.
- **Rule:** each gene is placed on a ladder, and the **worst tier present decides** the score.
- **Exception:** two or more *different* carbapenemases promote the plasmid to the top tier.

| Tier | Examples | Score |
|---|---|---|
| 2+ carbapenemases | KPC + NDM, NDM + OXA-48 | 1.00 |
| Metallo-carbapenemase | NDM, IMP, VIM | 0.95 |
| Serine carbapenemase | KPC, OXA-48-like, other OXA carbapenemases, GES/IMI | 0.90 |
| Colistin | *mcr* | 0.85 |
| ESBL | CTX-M, ESBL variants of SHV/TEM, PER, VEB | 0.70 |
| Plasmid-borne AmpC | CMY, DHA, FOX, ACC | 0.60 |
| Quinolone | *qnr*, *aac(6')-Ib-cr*, *oqxAB* | 0.45 |
| Other beta-lactam | TEM-1, SHV-1/-11, narrow OXA (OXA-1, OXA-10) | 0.30 |
| Any other AMR gene | *sul*, *dfr*, *tet*, aminoglycoside genes | 0.15 |
| None detected | | 0.00 |

SHV, TEM and OXA variants are resolved against explicit lists, because their numbering is not systematic. OXA-48 is a carbapenemase but OXA-1 is not, and SHV-12 is an ESBL but SHV-1 is not. An unknown variant falls into the lower tier.

### 2.2 Mobility: can it move?

- **Input:** MOB-suite's `predicted_mobility` call sets the base tier.
- **Rule:** small bonuses are added for transfer machinery that was actually detected.

| MOB-suite call | Base |
|---|---|
| Conjugative | 0.90 |
| Mobilisable | 0.60 |
| Indeterminate / missing | 0.35 |
| Non-mobilisable | 0.10 |

**Bonuses:** +0.04 if a relaxase was detected, +0.03 for an MPF system, +0.03 for an oriT. The total is capped at 1.0.

A missing call scores above "non-mobilisable" on purpose. A plasmid MOB-suite could not call is not proof that it cannot move.

### 2.3 Spread: how far has its lineage gone?

Spread is a property of the **lineage** (PLING subcommunity). It is shared by every plasmid in the lineage and counted over the *K. pneumoniae* plasmids in the scored set.

| Component | Weight |
|---|---|
| Number of countries | 0.35 |
| Number of One-Health compartments (human clinical, gut carriage, animal/food, environment, wastewater) | 0.35 |
| Lineage size (members in the subcommunity) | 0.30 |

Each component is rank-transformed to 0–1 first, so lineage size cannot dominate. The weighted lineage value then makes up 80 % of the score. The other 20 % is whether the plasmid **is a PLING hub**, because a hub, by definition, bridges lineages:

```
spread = 0.80 × lineage value + 0.20 × (1 if hub plasmid else 0)
```

## 3. The verdict

```
overall = 0.45 × AMR + 0.35 × Mobility + 0.20 × Spread       (percentiles, 0–100)
          + 15 if the potential wave flag fires (capped at 100)
```

If one score is missing, the remaining weights are re-normalised.

| Overall | Verdict | Plasmids in the set |
|---|---|---|
| < 30 | Low concern | 2,009 |
| 30–55 | Monitor | 3,831 |
| 55–75 | Concerning | 2,683 |
| ≥ 75 | High concern | 1,591 |

**Why the weights are unequal.** AMR is what this plasmid carries, and mobility is whether it can move. Spread describes its relatives and is the score most affected by who sequenced what, so it gets the smallest weight.

**Why a band and not a number.** The scores rest on a hand-written weight table and on PLSDB's sampling. "High concern" is defensible; "73.4 / 100" claims a precision that isn't there.

## 4. Potential waves 

This flag asks: has this plasmid's lineage been seen carrying something much worse than it carries now? It is computed **leave-one-out**, meaning over the plasmid's lineage-mates and excluding the plasmid itself. It fires only when **all four** conditions hold:

1. A lineage-mate carries a carbapenemase-grade gene (worst neighbour AMR ≥ 0.95).
2. This plasmid does not (own AMR < 0.50).
3. This plasmid is mobile (mobility > 0.60, so conjugative or well-equipped mobilisable).
4. At least 20 % of the lineage carries high-grade resistance (AMR ≥ 0.70), so it is not one odd neighbour.

It fires for **164 plasmids** (about 1.6 %). Plasmids alone in their lineage and hub plasmids never get it.

## 5. Community epidemiology panel, maps and timeline

These describe the plasmid's **whole PLING community in PLSDB**, including members from other genera, not only the *K. pneumoniae* plasmids drawn in the graph:

- members;
- countries (*Show on map*);
- host genera and species;
- years (*Show timeline*);
- the share of members with metadata.

> **Note:** the **Spread score** (section 2.3) is computed over the *K. pneumoniae* set only. The **community panel** counts all of PLSDB. They answer different questions: how widespread the lineage is among *K. pneumoniae*, and where the whole family has been seen.

## 6. Limitations

- **Sampling bias.** PLSDB reflects what has been sequenced. Clinical isolates from high-income countries dominate, so "spread" partly measures how often a lineage was sequenced.
- **Dates.** "Year" is the year the NCBI sequence record was created (`NUCCORE_CreateDate`). PLSDB has no true collection date.
- **Genotype only.** Gene presence is not expression, and the verdict is not a clinical prediction. It is a triage aid for deciding what to look at first.
- **The weights are judgements.** They are written down so they can be challenged and changed.
