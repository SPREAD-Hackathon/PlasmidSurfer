# PlasmidSurfer — a plasmid epidemiological risk score

A tool that takes an assembled plasmid and returns a 0–100 risk score for
antimicrobial-resistance epidemiology, trained on PLSDB.

---

## 1. The idea in one paragraph

A naive version of this project ends badly. If you define "risk" as *carries a
carbapenemase* and then train a model on gene content, the model learns
`blaKPC present → risk = 100` and you have built a very expensive `abricate`.
The circularity is invisible until someone asks about it at the poster.

So we split risk into two terms that are estimated in completely different ways:

```
RISK  =  DISSEMINATION POTENTIAL   ×   CARGO SEVERITY
         learned, from backbone        looked up, from AMR genes
         (the vehicle)                 (the payload)
```

* **Cargo severity** needs no model. AMRFinderPlus and CARD already tell you
  which drug classes a determinant defeats; a transparent weight table turns
  that into 0–1. Carbapenemases and *mcr* at the top, *sul*/*dfr* at the bottom.
* **Dissemination potential** is the interesting, genuinely predictive part:
  *given only the plasmid's backbone — replicon, relaxase, MPF, oriT,
  maintenance machinery, with all resistance genes masked — how good a vehicle
  is this?*

The question that makes this epidemiologically useful is the one the split lets
you ask: **if this backbone captures a carbapenemase tomorrow, how far will it
go?** No existing tool answers that.

---

## 2. What we have, what it's good for, and what is missing

### What you gave me

| File | Rows | What it actually is | Verdict |
|---|---|---|---|
| `typing.csv` | 72,556 | **MOB-suite output for all of PLSDB**, one row per plasmid: size, GC, replicon types, relaxase, MPF, oriT, predicted mobility, MOB clusters, observed host range, PMLST | **The backbone of the project.** Nearly every feature and the label come from here |
| `amr.tsv` | 251,138 | AMR/stress gene hits, long format | **Cargo term.** See the warning below |
| `plasmidfinder.csv` | 60,619 | PlasmidFinder Inc typing, 40,208 plasmids, 277 replicon alleles | Good complementary backbone features |
| `plsdb_mashdb_sim.tsv` | 525,604 | Headerless edge list of Mash-similar plasmid pairs, 32,336 plasmids | **The de-risking find — see §3** |
| `typing_markers.csv` | 333,649 | oriT/relaxase BLAST hits, long format | Mostly redundant with `typing.csv`; skip for now |
| `nuccore.csv` | 72,556 | **The master table.** The only bridge between accession-keyed and UID-keyed files | **The keystone — see §2.3** |
| `biosample.csv` | 29,222 | Country, lat/lon, ecosystem tags, disease tags | Country + One-Health compartment, via `nuccore.csv` |
| `taxonomy.csv` | 7,440 | NCBI lineage per TAXONOMY_UID | Real host taxonomy, via `nuccore.csv` |
| `assembly.csv` | 22,065 | Assembly status, coverage, dates | Nice-to-have QC |
| `nucc_identical.csv` | 3,906 | Only the *duplicated* nuccore records | Superseded by `nuccore.csv` |
| `disease_terms.csv` | 330 | SYMP ontology definitions | Not per-plasmid. Ignore |
| `protein.csv` (PGAP) | — | Per-gene annotation of every plasmid | **The route to an actual language model — see §9** |

### 2.1 Two warnings about `amr.tsv`

**It contains two annotators merged.** 201,920 rows from `amrfinderplus`
(27,777 plasmids) and 49,218 from `rgi`/CARD (17,450 plasmids), overlapping on
15,469 plasmids. They disagree on everything cosmetic: `blaKPC-2` vs `KPC-2`,
`mph(A)` vs `mphA`, uppercase `BETA-LACTAM` vs lowercase
`monobactam; carbapenem; cephalosporin; penam`. If you don't deduplicate,
dual-annotated plasmids get double gene counts and an inflated breadth of
co-resistance. `plrisk/cargo.py` normalises and merges them.

**Half the rows are not antibiotics.** ~50,000 hits are copper, mercury,
arsenic, tellurium and quaternary ammonium. They do not belong in a clinical
severity score. They are tracked separately as co-selection features
(`--include-stress`).

**Do not classify carbapenemases with a regex on gene names.** We tried; it is
a trap. OXA numbering is not systematic — OXA-48 is a carbapenemase, OXA-10 is
not, OXA-72 is (an OXA-24/40 variant), OXA-437 is (OXA-48-like). Severity is
therefore driven by the drug classes NCBI and CARD assign, which are curated by
people who track this full time.

### 2.2 `typing.csv` field coverage

```
predicted_mobility                100.0%     primary/secondary_cluster_id  100.0%
size, gc, num_contigs             100.0%     mash_neighbor_identification  100.0%
observed_host_range_ncbi_*         86.8%     rep_type(s)                    73.8%
relaxase_type(s)                   54.2%     mpf_type                       34.4%
orit_type(s)                       24.2%     PMLST_scheme                    6.2%
```

`predicted_mobility` splits 28,816 non-mobilizable / 22,588 conjugative /
21,152 mobilizable — a well-balanced, fully populated feature.

### 2.3 The master table, and what it unlocks

`nuccore.csv` is the keystone. Everything annotation-side (`typing.csv`,
`amr.tsv`, `plasmidfinder.csv`) is keyed on `NUCCORE_ACC`; everything
metadata-side (`biosample.csv`, `taxonomy.csv`) is keyed on integer UIDs.
`nuccore.csv` is the only bridge.

```
typing.csv / amr.tsv / plasmidfinder.csv      keyed on NUCCORE_ACC
                  |
             nuccore.csv        ACC <-> NUCCORE_UID, TAXONOMY_UID, BIOSAMPLE_UID
                  |
      +-----------+-----------+
 taxonomy.csv            biosample.csv
 (host lineage)          (country, ecosystem, disease)
```

**The join closes at 72,556 / 72,556 — 100%.** Resulting coverage:

| Field | Coverage |
|---|---|
| `TAXONOMY_species`, `NUCCORE_CreateDate` | 100.0% |
| `TAXONOMY_order` | 99.4% |
| `TAXONOMY_genus` | 99.3% |
| `LOCATION_name` → country | 79.6% |
| `ECOSYSTEM_tags` → One-Health compartment | 76.2% |
| `DISEASE_group` | 14.1% |

Top hosts: *Escherichia* 14,804, *Klebsiella* 12,213, *Salmonella* 3,524,
*Enterococcus* 3,231, *Staphylococcus* 3,180. Top countries: USA 2,492, UK
2,219, **Spain 1,818**, China 1,601. Records span 1982–2024.

Two fields need cleaning before use, both handled in `plrisk/metadata.py`:

**`LOCATION_name` is free text at wildly varying granularity.** Real values
include `USA`, `USA,VA`, `China,Chengdu,Sichuan`, the full geocoded
`Leipzig Markt, Markt, center, Mitte, Leipzig, Saxony, 04109, Germany`, and
raw coordinates `37.09N;121E`. Taken as-is there are 4,663 distinct values, so
counting them as countries would score a lineage seen in Seville and Madrid as
spanning two countries. `country_of()` resolves the country from any of these
shapes (country-first, country-last, embedded) and returns missing rather than
a guess for coordinates — better a gap than a wrong value in a diversity count.

**`ECOSYSTEM_tags` is 162 tags at mixed levels of abstraction.**
`compartments_of()` collapses them into One-Health compartments:
`human_clinical`, `gut_carriage`, `animal_food`, `environment`, `wastewater`,
`anthropogenic_other`. Note that the two most frequent tags — `host_associated`
(18,102) and `cell_culture` (6,849) — deliberately map to nothing: they say a
sample came from a host or a culture, not which compartment it belongs to, and
guessing would manufacture signal.

**One caveat to state on a slide:** `NUCCORE_CreateDate` is the NCBI *record
creation* date, not a collection date. It lags isolation and clusters around
submission campaigns. Fine as a coarse "how long has this lineage been in the
record" proxy; not evidence about any particular year.

### 2.4 "We wanted to use the sequences"

You don't need them for v1, and for a 27-hour deadline that is good news.
Everything above is a *derived representation* of those sequences produced by
tools that are better at it than anything you'd train tomorrow. MOB-suite has
already found the replicons and relaxases; PGAP has already called the genes.
Raw FASTA becomes necessary only for §9.

---

## 3. The single most important methodological point

**PLSDB is full of near-duplicates.** The same epidemic plasmid appears in
dozens of isolates. The largest Mash component has 467 members; the largest
MOB-suite cluster has 1,703.

If you split train/test at random, a plasmid and its 40 near-identical siblings
end up on both sides, and the model gets graded on sequences it has effectively
memorised. We measured it on this exact dataset:

```
random 5-fold CV     R² = 0.485      <-- wrong, and meaningless
grouped 5-fold CV    R² = 0.078      <-- real
                     inflation = +0.407
```

**Put that comparison on a slide.** It is the most convincing methods point
available to you, it takes ten seconds to explain, and most hackathon projects
will have the inflated number without knowing it. Reproduce it with
`python train.py --random-cv-demo`.

Groups come from three sources, in priority order:

1. **PLING communities** — best (DCJ-indel distance respects rearrangements).
   Drop the file in and the pipeline picks it up automatically.
2. **Mash components** — connected components of `plsdb_mashdb_sim.tsv`.
   Covers 28,141 plasmids in 6,176 components. **Available right now**, which
   is why you are not blocked on PLING finishing.
3. **MOB-suite `primary_cluster_id`** — 100% coverage, used as the fallback for
   the remaining 44,415 plasmids.

Total: 72,556 plasmids in 12,900 lineages.

---

## 4. The label

There is no `risk` column in PLSDB, so the target is **observed dissemination
success of a lineage**, computed per group.

**Always available (`--label taxonomic`):**

| Component | Meaning |
|---|---|
| `group_size` | distinct PLSDB records in the lineage — direct evidence of repeated independent recovery |
| `n_host_genera` | distinct host genera — cross-genus spread is the epidemiologically alarming part |
| `host_rank_max` | MOB-suite `observed_host_range_ncbi_rank`, ordinal from species → multi-phylla |

**Added with `nuccore.csv` (`--label full`, the default):**

| Component | Meaning |
|---|---|
| `n_countries` | geographic spread |
| `n_compartments` | distinct One-Health compartments — **the clinical–community–environment axis** |
| `year_span` | persistence in the record over time |

With the metadata join, the host genus also switches from the MOB-suite
nearest-neighbour *proxy* to real NCBI taxonomy for 99.3% of plasmids.

Each component is rank-transformed to [0,1], averaged, and converted to a
percentile. A score of 80 means *"this lineage has spread more widely than 80%
of lineages in PLSDB"*. Rank-transforming rather than z-scoring is deliberate:
these are all heavily right-skewed count variables and only the ordering is
meaningful.

### Small lineages are handled separately, and here is why

A lineage with two members can only ever show one or two countries, so its
`n_countries` is driven by lineage *size* rather than by real geographic
spread. Feeding that into the composite injects noise exactly where the label
is already weakest — and 7,072 of 12,900 lineages have fewer than three
members.

So for lineages below `min_members_for_metadata` (default 3), the composite is
built from the always-available components only, then rescaled onto the same
percentile scale. This was not a theoretical worry: the first version without
it produced a label the model could not learn at all (see §7).

### The bias you must state out loud

PLSDB reflects **what has been sequenced**, not what circulates. Clinical
Enterobacterales from high-income countries are massively over-represented;
environmental and African isolates are nearly absent. So "dissemination success"
partly measures *success at being sequenced*.

This does not invalidate the model. It bounds the claim: it predicts spread
*within the sampled world*. Say so before a reviewer says it for you — it reads
as rigour, not weakness.

---

## 5. Features — and the leakage guard

**Allowed** (the vehicle): `log_size`, `gc`, `num_contigs`, `is_fragmented`,
replicon types (MOB-suite `rep_cluster` + PlasmidFinder Inc families), relaxase
(MOB) type, MPF type, oriT presence, `predicted_mobility` one-hot, PMLST scheme
presence, and counts of each. **374 features** after dropping categorical
levels seen fewer than 30 times.

**Banned, and actively asserted against** in `features._assert_no_leakage()`,
which raises rather than letting a tautology train silently:

* `primary_cluster_id` / `secondary_cluster_id` — these *are* the grouping, so
  they encode group size, which is a third of the label
* `mash_nearest_neighbor` / `mash_neighbor_identification` — the host proxy the
  label is built from
* `observed_host_range_*` / `predicted_host_range_*` — literally the label
* **every AMR gene** — this is the cargo term, multiplied in afterwards. If the
  model could see it, it would learn "carbapenemase = dangerous".

**Arguable case:** metal and biocide resistance. Genuine co-selection drivers,
but they come from the same file as the cargo. Off by default, enabled with
`--include-stress`. Run it both ways and report the difference.

---

## 6. Answering your actual question: how to subset for 1,000 plasmids

### Rule 1 — compute labels on all 72,556, subset only for training

The label is a *lineage-level* property. If you subset first, a family that is
"6 genera across 40 assemblies" becomes "2 records, 1 genus" and you have
destroyed the signal you're trying to learn. **Subset last, never first.** The
pipeline enforces this.

### Rule 2 — 1,000 plasmids is the wrong unit. Count lineages.

This is the most important practical thing in this document. Your effective
sample size is the number of independent **lineages**, not plasmids. 1,000
plasmids drawn naively gives you maybe 40 lineages, because epidemic plasmids
dominate. With the per-lineage cap we get **1,007 plasmids spanning 326
lineages** — and 326 is still small.

### Rule 3 — cap members per lineage, and weight by 1/n

Without a cap, three epidemic plasmids supply most training rows and the loss
is dominated by memorising them. We keep at most 5 members per lineage and
weight each row `1/n_members`, so every lineage contributes total weight 1.

### Rule 4 — don't filter hard on group size

Tempting, looks like QC, but group size is **one third of the label**, so
filtering on it is range restriction on the target in disguise:

| min group size | Enterobacterales n | label mean | label SD |
|---|---|---|---|
| ≥1 | 38,563 | 87.9 | 16.2 |
| ≥2 | 38,182 | 88.5 | 15.1 |
| ≥3 | 34,666 | 91.8 | **10.5** |
| ≥5 | 31,337 | 93.6 | **8.1** |

At ≥3 the label collapses into the top quartile. We use **≥2** — enough to drop
pure singletons, not enough to flatten the target.

### Rule 5 — focus on Enterobacterales

Concentrates the signal where the clinically relevant conjugative plasmids are,
removes cyanobacterial and plant-symbiont plasmids whose "dissemination" means
something different, and makes host range comparable across rows. 34,780 of
61,253 eligible plasmids.

### Rule 6 — stratify the sample across label deciles

A uniform draw returns mostly mid-range lineages. Sampling whole groups
stratified by label decile gives a usable spread: mean 71.5, SD 19.5, range
17.2–100.0.

### My actual recommendation: run 1,000 first, then scale up immediately

I ran both. The comparison is the most useful result in this repo:

| training set | lineages | plasmids | R² (grouped CV) | Spearman ρ |
|---|---|---|---|---|
| `--subset 1000` | 326 | 1,007 | 0.078 | 0.373 |
| all Enterobacterales | **4,175** | 13,919 | **0.241** | **0.402** |

**Tripling R² by adding lineages means you are data-limited, not
signal-limited.** The backbone genuinely carries information about
dissemination; 326 lineages is just not enough to extract it.

So: build the 1,000-plasmid version first to have an end-to-end result early,
then run the full set — it takes about ten minutes and costs you nothing. Show
both. "We checked whether we were data-limited, and we were" is a better slide
than either number alone.

---

## 7. Results

**1,000-plasmid subset — 326 lineages**

```
Grouped 5-fold CV:   R² = 0.078    MAE = 14.6    Spearman ρ = 0.373
  mean only          R² = -0.001                 ρ = -0.043
  size only          R² = -0.044                 ρ = +0.066
  size + GC + mob    R² = -0.152                 ρ = +0.171
  FULL MODEL         R² = +0.078                 ρ = +0.373
```

**All Enterobacterales — 4,175 lineages, 13,919 plasmids**

```
Grouped 5-fold CV:   R² = 0.241    MAE = 14.2    Spearman ρ = 0.402
  mean only          R² = -0.000                 ρ = -0.010
  size only          R² = +0.020                 ρ = +0.157
  FULL MODEL         R² = +0.241                 ρ = +0.402
  model adds R² +0.221 over size alone
```

**Read this honestly.** ρ ≈ 0.40 on held-out lineages is a real but modest
signal. The model is not a finished predictor. What it *is*:

* genuinely better than every baseline, including plasmid size alone (+0.22 R²)
  — so it is a plasmid model, not a ruler
* evaluated in a way most of the room will not have bothered with
* clearly still data-limited: R² tripled from 326 → 4,175 lineages, so the
  curve has not flattened and more lineages (or PLING's better grouping) should
  help further

Top features by permutation importance: `n_rep_types`, `log_size`, `gc`,
`n_inc_groups`, `rep__IncC`, `mob__MOBQ`, `n_relaxases`, `has_orit`,
`rep__IncI_gamma_K1`, `inc__IncI1_I`.

That IncC, IncI1 and MOBQ surface near the top is a reassuring sanity check —
those are exactly the backbones epidemiologists worry about.

**Do not quote a number without saying it is cross-validated by lineage.** A
bare "ρ = 0.37" invites someone to compare it against a leaky 0.9 from another
team.

---

## 8. Running it

```bash
pip install pandas numpy scikit-learn scipy        # that's all; no GPU, no torch

# point config at your files (or edit plrisk/config.py)
export PLRISK_TYPING=data/typing.csv
export PLRISK_AMR=data/amr.tsv
export PLRISK_PLASMIDFIND=data/plasmidfinder.csv
export PLRISK_MASH=data/plsdb_mashdb_sim.tsv

# 1. build the dataset (labels on all 72k, training flag on the subset)
python build_dataset.py --subset 1000 --focus enterobacterales

# 2. train + evaluate, including the leakage demonstration
python train.py --random-cv-demo

# 3. score a new plasmid
mob_typer --infile my_plasmid.fasta --out_file my_typing.txt
amrfinder -n my_plasmid.fasta --plus -o my_amr.tsv
python score.py --typing my_typing.txt --amr my_amr.tsv

# or score something already in PLSDB
python score.py --from-plsdb NZ_CP139252.1
```

Full run, all Enterobacterales lineages:

```bash
python build_dataset.py --focus enterobacterales --out out_full
python train.py --dataset out_full/dataset.csv.gz \
                --features out_full/feature_columns.txt --out out_full
```

### Swapping in PLING

Write a two-column file (plasmid id, community id) to
`data/pling_communities.tsv` and re-run `build_dataset.py`. It is picked up
automatically and takes priority over Mash. Filename suffixes like `.fna` are
stripped. Nothing else changes — **this is why you are not blocked on PLING.**

### Output

```
out/dataset.csv.gz           all 72,556 plasmids, features + labels + cargo
out/feature_columns.txt      the 374 feature names, in order
out/model.pkl                trained model + feature order + cargo floor
out/metrics.json             grouped CV, per fold, and all baselines
out/feature_importance.csv   permutation importance
out/oof_predictions.csv      out-of-fold predictions for plotting
out/figures/*.png            slide-ready figures (python figures.py)
```

### Figures

```bash
python figures.py                      # reads out/
python figures.py --out sweep_full_all # any run directory
```

| File | The question it answers |
|---|---|
| `fig1_leakage.png` | Why grouped CV. **The methods slide.** |
| `fig2_label.png` | What "dissemination" is actually made of |
| `fig3_performance.png` | Predicted vs observed, and the baseline comparison |
| `fig4_importance.png` | Which backbone features the model reads |
| `fig5_risk_space.png` | The two-component risk plane. **The product slide.** |

---

## 9. If you want an actual language model (the stretch goal)

Your original instinct was right, and `protein.csv` from PGAP is the way in —
but **not** by tokenising nucleotides. A 150 kb plasmid does not fit in any
reasonable context window, and pretraining a DNA LM overnight will not happen.

The published approach treats a plasmid as a **sentence whose words are genes**.
PLASMe does exactly this: it encodes plasmid sequences as a language over a
protein-cluster token set, so a Transformer learns the importance of proteins
and their correlations through positional token embeddings and attention. That
is precisely the DeepLoc analogy you started from, with the "amino acid"
promoted to "gene".

Concretely, if there is time on day 2:

1. `mmseqs2 easy-cluster` the PGAP proteins at ~50% identity → a vocabulary of
   protein families
2. each plasmid becomes an ordered token sequence over that vocabulary
3. either a small attention-pooling model trained from scratch, or
   PlasmidGPT's pretrained embeddings as a drop-in (it is on HuggingFace with a
   ready `embeddings.py`; trained on 153,208 Addgene plasmids, and its
   representations are reported to generalise to natural plasmids — mind the
   engineered-vs-natural domain shift)
4. same grouped CV, same baselines, compared head-to-head against the
   gradient-boosting model here

**Do not throw away the tabular model to do this.** With ~4,000 lineages the
boosted trees will likely still win, and a credible negative result
("embeddings did not beat engineered features at this sample size") is a far
better slide than a transformer with no baseline.

---

## 10. Limitations slide — write these down before someone else does

1. **Sampling bias.** PLSDB measures what gets sequenced. The label partly
   measures success at being sequenced. Bounds the claim to the sampled world.
2. **No geography, no time, no ecosystem yet.** Blocked on the missing master
   table (§2.3). The current label is taxonomic breadth, not full
   epidemiological spread.
3. **Host taxonomy is a proxy.** Derived from MOB-suite's nearest-neighbour
   identification, not the real NCBI host of the record.
4. **The severity table is a judgement call.** Transparent and arguable by
   design — that is a feature, but it is not objective.
5. **Modest performance.** ρ ≈ 0.37 on 326 lineages. Useful for triage and
   ranking; not a diagnostic.
6. **A score is not a probability.** It is a percentile against PLSDB lineages.
   `score.py` always prints the two components separately, because a 60 from
   "mobile backbone, trivial cargo" and a 60 from "immobile plasmid carrying
   NDM-1" demand completely different responses.
7. **Correlation, not causation.** Backbone features that predict observed
   spread are not proven *causes* of spread. Conjugation assays would be the
   experiment.

---

## Module map

```
plrisk/config.py     paths + tunable constants — the only file to edit
plrisk/metadata.py   nuccore join, country normalisation, One-Health mapping
plrisk/grouping.py   PLING / Mash / MOB-suite lineage assignment
plrisk/labels.py     dissemination success, per lineage
plrisk/features.py   backbone feature matrix + the leakage guard
plrisk/cargo.py      AMR severity, two-annotator deduplication
build_dataset.py     step 1 — assemble everything, mark the training subset
train.py             step 2 — grouped CV, baselines, importance, final fit
score.py             step 3 — score a new plasmid
```
