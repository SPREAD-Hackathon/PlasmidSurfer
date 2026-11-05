"""
Cargo severity: how clinically bad is what this plasmid is carrying?

TWO THINGS ABOUT amr.tsv YOU MUST KNOW BEFORE USING IT
-------------------------------------------------------
1. IT CONTAINS TWO ANNOTATORS, MERGED.
       amrfinderplus : ~202k hits over ~27.8k plasmids.
                       `drug_class` is UPPERCASE and coarse ("BETA-LACTAM");
                       `antimicrobial_agent` carries the useful detail
                       ("CARBAPENEM"). Gene symbols keep the bla prefix
                       (blaKPC-2).
       rgi (CARD)    : ~49k hits over ~17.5k plasmids.
                       `drug_class` is lowercase and semicolon-separated
                       ("monobactam; carbapenem; cephalosporin; penam");
                       `antimicrobial_agent` is empty. Gene symbols drop the
                       bla prefix (KPC-2).
   ~15.5k plasmids are annotated by BOTH. If you do not deduplicate, those
   plasmids get double the apparent gene count and an inflated "breadth of
   co-resistance". This module normalises gene names (strips bla, uppercases)
   and takes the union, so each determinant is counted once.

2. HALF THE ROWS ARE NOT ANTIBIOTICS. Copper, mercury, arsenic, tellurium and
   quaternary ammonium account for ~50k hits. They do not belong in a clinical
   severity score. They are tracked separately as co-selection features.

DESIGN DECISION: SEVERITY COMES FROM THE PHENOTYPE, NOT FROM A GENE REGEX
--------------------------------------------------------------------------
An earlier version of this file matched gene names with regexes. That was a
mistake and it is worth saying why, because it is a trap anyone writing this
would fall into: OXA numbering is not systematic. blaOXA-48 is a carbapenemase,
blaOXA-10 is not, blaOXA-72 is (it is an OXA-24/40 variant), blaOXA-437 is
(OXA-48-like). No regex over OXA numbers is maintainable, and getting it wrong
silently scores a carbapenemase as an ordinary beta-lactamase.

So severity is driven by the drug classes NCBI and CARD themselves assign.
Those curations are maintained by people who track the literature full time.
A small gene-level override table remains, but only for determinants whose
class strings genuinely understate them.

The severity numbers below are a scientific judgement, not a constant of
nature. Put the table on a slide and let people argue with it.
"""

from __future__ import annotations

import re

import pandas as pd

# ---------------------------------------------------------------------------
# Severity by implicated drug class / agent, 0-1.
# Matched case-insensitively as substrings against both the AMRFinderPlus
# `antimicrobial_agent` + `drug_class` and the CARD `drug_class` strings.
# A determinant's severity is the MAXIMUM over every class it matches.
# ---------------------------------------------------------------------------
CLASS_SEVERITY: list[tuple[str, float]] = [
    ("carbapenem",        1.00),   # loss of the last broadly usable beta-lactam
    ("colistin",          0.95),
    ("peptide antibiotic", 0.80),  # CARD's label for mcr and friends
    ("glycopeptide",      0.90),   # vanA/B
    ("vancomycin",        0.90),
    ("oxazolidinone",     0.85),   # linezolid: cfr, optrA, poxtA
    ("linezolid",         0.85),
    ("glycylcycline",     0.80),   # tigecycline: tet(X)
    ("tigecycline",       0.80),
    ("cephalosporin",     0.65),   # ESBL / AmpC territory
    ("cephamycin",        0.60),
    ("monobactam",        0.55),
    ("penem",             0.50),
    ("fluoroquinolone",   0.45),
    ("quinolone",         0.40),
    ("fosfomycin",        0.40),
    ("rifamycin",         0.35),
    ("aminoglycoside",    0.30),
    ("macrolide",         0.25),
    ("penam",             0.25),
    ("beta-lactam",       0.25),   # unqualified: narrow-spectrum by default
    ("lincosamide",       0.20),
    ("streptogramin",     0.20),
    ("phenicol",          0.20),
    ("tetracycline",      0.15),
    ("diaminopyrimidine", 0.10),
    ("trimethoprim",      0.10),
    ("sulfonamide",       0.10),
    ("nitroimidazole",    0.15),
    ("bleomycin",         0.05),
    ("streptothricin",    0.05),
    ("kasugamycin",       0.02),
]

# Gene-level overrides, applied AFTER the class lookup and taking the maximum.
# Only for determinants whose curated class strings understate them.
GENE_OVERRIDE: list[tuple[str, float]] = [
    (r"^MCR-?\d",                 0.95),   # mobile colistin resistance
    (r"^VAN[ABDEGMN]\b",          0.90),
    (r"^(OPTRA|CFR|POXTA)\b",     0.85),
    (r"^TET\(?X",                 0.80),
    (r"^(ARMA|RMT[A-H]|NPMA)\b",  0.70),   # 16S methylases: pan-aminoglycoside
]

# Non-antibiotic determinants: excluded from severity, counted separately.
STRESS_TOKENS = (
    "copper", "silver", "mercury", "arsenic", "nickel", "tellurium", "zinc",
    "cadmium", "quaternary ammonium", "biocide", "gold", "antimony",
    "chromium", "disinfecting agents", "efflux", "organomercury",
    "phenylmercury", "acid",
)

_OVERRIDE = [(re.compile(p, re.I), w) for p, w in GENE_OVERRIDE]


def normalise_gene(symbol) -> str:
    """
    Display name: strip the bla prefix and uppercase.

    'blaKPC-2' -> 'KPC-2'. This is what gets shown to the user.
    """
    s = str(symbol).strip()
    if not s or s.lower() in {"nan", "none", "na", "-"}:
        return ""
    s = re.sub(r"^bla", "", s, flags=re.I)
    return s.upper()


def dedup_key(gene_norm: str) -> str:
    """
    Matching key used to merge the two annotators' calls for the same gene.

    AMRFinderPlus and CARD punctuate differently -- mph(A) vs mphA,
    msr(E) vs msrE, aac(6')-Ib vs AAC(6')-IB -- so all parentheses, primes,
    hyphens and spaces are removed before comparison.

    KNOWN LIMITATION: this does not merge allele-numbering differences such as
    rmtB vs rmtB1, or aadA vs aadA1. Those still count as two determinants, so
    `n_amr_genes` is a slight over-count for dual-annotated plasmids. It does
    not affect `cargo_severity` meaningfully, because severity is driven by the
    single worst determinant plus a saturating class-breadth term, and the
    duplicated pair shares a class. Worth fixing if you ever report gene counts
    as a headline number.
    """
    return re.sub(r"[^A-Z0-9/]", "", gene_norm.upper())


def _class_text(drug_class, agent) -> str:
    bits = [str(x) for x in (drug_class, agent)
            if x is not None and str(x).lower() not in {"nan", "none", ""}]
    return " ; ".join(bits).lower()


def is_stress(class_text: str) -> bool:
    return any(tok in class_text for tok in STRESS_TOKENS)


def severity_of(gene_norm: str, class_text: str) -> float:
    """Severity of one determinant: worst implicated class, then overrides."""
    if is_stress(class_text):
        return 0.0
    sev = 0.0
    for token, w in CLASS_SEVERITY:
        if token in class_text:
            sev = max(sev, w)
    for rx, w in _OVERRIDE:
        if rx.search(gene_norm):
            sev = max(sev, w)
    return sev


def worst_class_of(class_text: str) -> str:
    """Human-readable label for the most serious class implicated."""
    best, best_w = "unclassified", -1.0
    for token, w in CLASS_SEVERITY:
        if token in class_text and w > best_w:
            best, best_w = token, w
    return best


def combine(severities) -> float:
    """
    Collapse per-determinant severities into one 0-1 cargo score.

    Dominated by the single worst determinant -- one NDM already makes a
    plasmid maximally concerning -- plus a saturating bonus for breadth of
    co-resistance. This avoids the classic failure where ten tet/sul genes
    outrank one carbapenemase.
    """
    sev = [s for s in severities if s > 0]
    if not sev:
        return 0.0
    worst = max(sev)
    breadth = 1.0 - 0.85 ** (len(set(round(s, 3) for s in sev)) - 1)
    return float(min(1.0, worst + 0.15 * breadth * (1.0 - worst)))


def cargo_table(amr: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse long-format AMR annotation to one row per plasmid.

    Handles the amrfinderplus/rgi mixture and deduplicates determinants.

    Input  : amr.tsv as a DataFrame. Needs NUCCORE_ACC, gene_symbol and
             drug_class; antimicrobial_agent is used when present.
    Output : indexed by NUCCORE_ACC, columns
             cargo_severity, n_amr_genes, n_amr_classes, n_stress_genes,
             worst_class, amr_gene_list, n_annotators
    """
    df = amr.copy()
    if "antimicrobial_agent" not in df.columns:
        df["antimicrobial_agent"] = None
    if "analysis_software_name" not in df.columns:
        df["analysis_software_name"] = "unknown"

    df["gene_norm"] = df["gene_symbol"].map(normalise_gene)
    df = df[df.gene_norm != ""]
    df["gene_key"] = df["gene_norm"].map(dedup_key)
    df["class_text"] = [
        _class_text(dc, ag)
        for dc, ag in zip(df["drug_class"], df["antimicrobial_agent"])
    ]

    # --- DEDUPLICATION ----------------------------------------------------
    # One row per (plasmid, determinant), matched on the punctuation-free key.
    # Where both annotators called the same gene, their class strings are
    # concatenated so neither annotator's detail is lost, and the determinant
    # is scored once.
    merged = (df.groupby(["NUCCORE_ACC", "gene_key"], sort=False)
                .agg(gene_norm=("gene_norm", "first"),
                     class_text=("class_text", lambda s: " ; ".join(sorted(set(s)))),
                     n_annotators=("analysis_software_name", "nunique"))
                .reset_index())

    merged["severity"] = [severity_of(g, c)
                          for g, c in zip(merged.gene_norm, merged.class_text)]
    merged["stress"] = merged.class_text.map(is_stress)
    merged["worst_class"] = merged.class_text.map(worst_class_of)

    abx = merged[~merged.stress]
    rows = []
    for acc, t in abx.groupby("NUCCORE_ACC", sort=False):
        rows.append({
            "NUCCORE_ACC": acc,
            "cargo_severity": combine(t.severity.tolist()),
            "n_amr_genes": int(t.gene_norm.nunique()),
            "n_amr_classes": int(t.worst_class.nunique()),
            "worst_class": t.loc[t.severity.idxmax(), "worst_class"],
            "amr_gene_list": ";".join(sorted(set(t.gene_norm))),
            "n_annotators": int(t.n_annotators.max()),
        })
    out = pd.DataFrame(rows).set_index("NUCCORE_ACC") if rows else pd.DataFrame(
        columns=["cargo_severity", "n_amr_genes", "n_amr_classes",
                 "worst_class", "amr_gene_list", "n_annotators"])

    stress_counts = (merged[merged.stress]
                     .groupby("NUCCORE_ACC").gene_norm.nunique())

    # Plasmids that appear only via stress genes still get a row (severity 0).
    only_stress = set(stress_counts.index) - set(out.index)
    if only_stress:
        extra = pd.DataFrame(index=sorted(only_stress))
        extra["cargo_severity"] = 0.0
        extra["n_amr_genes"] = 0
        extra["n_amr_classes"] = 0
        extra["worst_class"] = "none"
        extra["amr_gene_list"] = ""
        extra["n_annotators"] = 1
        out = pd.concat([out, extra])

    out["n_stress_genes"] = (stress_counts.reindex(out.index)
                             .fillna(0).astype(int))
    out.index.name = "NUCCORE_ACC"
    return out
