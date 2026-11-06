"""
Joining PLSDB's metadata tables, and cleaning the two messy fields.

THE JOIN
--------
`nuccore.csv` is the master table and the only bridge between the
accession-keyed annotation files and the UID-keyed metadata files:

    typing.csv / amr.tsv / plasmidfinder.csv      keyed on NUCCORE_ACC
                      |
                 nuccore.csv        NUCCORE_ACC <-> NUCCORE_UID,
                      |                            TAXONOMY_UID, BIOSAMPLE_UID
          +-----------+-----------+
     taxonomy.csv            biosample.csv
     (host lineage)          (country, ecosystem, disease)

With it, coverage over the 72,556 plasmids in typing.csv is:

    TAXONOMY_species     100.0%      NUCCORE_CreateDate   100.0%
    TAXONOMY_order        99.4%      LOCATION_name         79.6%
    TAXONOMY_genus        99.3%      ECOSYSTEM_tags        76.2%
                                     DISEASE_group         14.1%

TWO FIELDS NEED CLEANING BEFORE THEY ARE USABLE

1. LOCATION_name is free text at wildly varying granularity. Real examples:
       'USA'
       'USA,VA'
       'China,Chengdu,Sichuan'
       'Leipzig Markt, Markt, center, Mitte, Leipzig, Saxony, 04109, Germany'
       '37.09N;121E'
   Taken raw there are 4,663 distinct values, so counting them as "countries"
   would score a lineage seen in Seville and Madrid as spanning two countries.
   `country_of()` resolves a country name from any of these shapes.

2. ECOSYSTEM_tags is a comma-separated set drawn from 162 tags, at mixed
   levels of abstraction ('host_associated', 'blood', 'wastewater', 'soil').
   `compartments_of()` collapses them into One-Health compartments.

A CAVEAT YOU MUST STATE: NUCCORE_CreateDate IS NOT A COLLECTION DATE
---------------------------------------------------------------------
It is the date the record was created at NCBI, which can lag isolation by
years and clusters around submission campaigns. Used here as a coarse proxy
for "how long has this lineage been circulating in the record", which is
defensible for a time *span* but not for any claim about a specific year.
Label it as such on any figure.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Country resolution
# ---------------------------------------------------------------------------
# Canonical names plus the aliases that actually appear in PLSDB.
_COUNTRY_ALIASES = {
    "usa": "USA", "united states": "USA", "united states of america": "USA",
    "u.s.a.": "USA", "us": "USA",
    "uk": "United Kingdom", "united kingdom": "United Kingdom",
    "great britain": "United Kingdom", "england": "United Kingdom",
    "scotland": "United Kingdom", "wales": "United Kingdom",
    "northern ireland": "United Kingdom",
    "south korea": "South Korea", "republic of korea": "South Korea",
    "korea": "South Korea", "korea, south": "South Korea",
    "north korea": "North Korea",
    "russia": "Russia", "russian federation": "Russia",
    "viet nam": "Vietnam", "vietnam": "Vietnam",
    "czech republic": "Czech Republic", "czechia": "Czech Republic",
    "ivory coast": "Cote d'Ivoire", "cote d'ivoire": "Cote d'Ivoire",
    "côte d'ivoire": "Cote d'Ivoire",
    "burma": "Myanmar", "myanmar": "Myanmar",
    "holland": "Netherlands", "netherlands": "Netherlands",
    "the netherlands": "Netherlands",
    "uae": "United Arab Emirates",
    "united arab emirates": "United Arab Emirates",
    "hong kong": "China", "macau": "China", "macao": "China",
    "taiwan": "Taiwan",
    "republic of the congo": "Congo",
    "democratic republic of the congo": "DR Congo",
    "drc": "DR Congo",
    "bosnia and herzegovina": "Bosnia and Herzegovina",
    "trinidad and tobago": "Trinidad and Tobago",
    "papua new guinea": "Papua New Guinea",
    "new zealand": "New Zealand",
    "saudi arabia": "Saudi Arabia",
    "south africa": "South Africa",
    "sri lanka": "Sri Lanka",
    "costa rica": "Costa Rica",
    "puerto rico": "USA",
    "el salvador": "El Salvador",
    "dominican republic": "Dominican Republic",
    "burkina faso": "Burkina Faso",
    "sierra leone": "Sierra Leone",
    "south sudan": "South Sudan",
    "cape verde": "Cabo Verde",
    "swaziland": "Eswatini",
    "macedonia": "North Macedonia",
    "north macedonia": "North Macedonia",
    "timor-leste": "Timor-Leste", "east timor": "Timor-Leste",
}

_COUNTRIES = """
Afghanistan Albania Algeria Andorra Angola Argentina Armenia Australia Austria
Azerbaijan Bahamas Bahrain Bangladesh Barbados Belarus Belgium Belize Benin
Bhutan Bolivia Botswana Brazil Brunei Bulgaria Burundi Cambodia Cameroon Canada
Chad Chile China Colombia Comoros Congo Croatia Cuba Cyprus Denmark Djibouti
Dominica Ecuador Egypt Eritrea Estonia Eswatini Ethiopia Fiji Finland France
Gabon Gambia Georgia Germany Ghana Greece Greenland Grenada Guatemala Guinea
Guyana Haiti Honduras Hungary Iceland India Indonesia Iran Iraq Ireland Israel
Italy Jamaica Japan Jordan Kazakhstan Kenya Kosovo Kuwait Kyrgyzstan Laos
Latvia Lebanon Lesotho Liberia Libya Liechtenstein Lithuania Luxembourg
Madagascar Malawi Malaysia Maldives Mali Malta Mauritania Mauritius Mexico
Moldova Monaco Mongolia Montenegro Morocco Mozambique Myanmar Namibia Nepal
Netherlands Nicaragua Niger Nigeria Norway Oman Pakistan Palau Palestine Panama
Paraguay Peru Philippines Poland Portugal Qatar Romania Rwanda Samoa Senegal
Serbia Seychelles Singapore Slovakia Slovenia Somalia Spain Sudan Suriname
Sweden Switzerland Syria Taiwan Tajikistan Tanzania Thailand Togo Tonga Tunisia
Turkey Turkmenistan Uganda Ukraine Uruguay Uzbekistan Vanuatu Venezuela Vietnam
Yemen Zambia Zimbabwe
""".split()

_COUNTRY_LOOKUP: dict[str, str] = {c.lower(): c for c in _COUNTRIES}
_COUNTRY_LOOKUP.update(_COUNTRY_ALIASES)

_COORD_RE = re.compile(r"^\s*-?\d+(\.\d+)?\s*[NS]?\s*[;,]\s*-?\d+(\.\d+)?\s*[EW]?\s*$",
                       re.I)


def country_of(location) -> str | float:
    """
    Pull a country out of a PLSDB LOCATION_name string.

    Handles 'USA', 'USA,VA', 'China,Chengdu,Sichuan' (country first) and
    long geocoded addresses ending in the country (country last). Returns NaN
    for raw coordinates and for anything unrecognised -- better a missing
    value than a wrong one, since these feed a diversity count.
    """
    if not isinstance(location, str) or not location.strip():
        return np.nan
    s = location.strip()
    if _COORD_RE.match(s):
        return np.nan

    parts = [p.strip().lower() for p in re.split(r"[,;]", s) if p.strip()]
    if not parts:
        return np.nan
    # Country-first is the most common PLSDB shape, country-last covers
    # fully geocoded addresses. Check both ends before scanning the middle.
    for cand in (parts[0], parts[-1]):
        hit = _COUNTRY_LOOKUP.get(cand)
        if hit:
            return hit
    for p in parts:
        hit = _COUNTRY_LOOKUP.get(p)
        if hit:
            return hit
    # Last resort: a country name embedded in a longer free-text token.
    low = s.lower()
    for name, canon in _COUNTRY_LOOKUP.items():
        if len(name) > 4 and re.search(rf"\b{re.escape(name)}\b", low):
            return canon
    return np.nan


# ---------------------------------------------------------------------------
# One-Health compartments
# ---------------------------------------------------------------------------
# PLSDB's 162 ecosystem tags collapsed into the compartments that matter for
# AMR surveillance across the clinical-community-environment interface.
# A sample can belong to several; that is the point.
COMPARTMENT_TAGS: dict[str, set[str]] = {
    "human_clinical": {
        "blood", "circulatory_system", "urinary_system", "respiratory_system",
        "wound", "abscess", "nasal", "skin", "oral", "trachea", "bronchi",
        "lung", "liver", "kidney", "nervous_system", "vaginal", "mamary",
        "hospital", "disease", "cerebrospinal", "sputum", "catheter",
        "peritoneal", "pus", "biopsy", "clinical",
    },
    "gut_carriage": {
        "fecal", "gastrointestinal_system", "rectal", "cecal", "stool",
    },
    "animal_food": {
        "food", "meat", "farm", "fermented", "drink", "aquaculture",
        "poultry", "livestock", "dairy", "milk", "egg", "feed",
    },
    "environment": {
        "soil", "aquatic", "freshwater", "marine", "sediment", "river",
        "lake", "sea", "wetland", "estuary", "pond", "spring", "terrestrial",
        "rhizosphere", "saline", "groundwater", "plant", "air", "dust",
    },
    "wastewater": {
        "wastewater", "wwtp", "sewage", "sludge", "effluent", "manure",
    },
    "anthropogenic_other": {
        "anthropogenic", "contaminated", "construction", "artificial",
        "laboratory", "industrial", "building",
    },
}

_TAG_TO_COMPARTMENT = {t: c for c, ts in COMPARTMENT_TAGS.items() for t in ts}


def compartments_of(tags) -> str:
    """
    Map a comma-separated ECOSYSTEM_tags string to One-Health compartments.

    Returns a comma-separated set, or '' if nothing maps. Note that
    'host_associated' and 'cell_culture', the two most frequent tags, map to
    nothing on purpose: they say a sample came from a host or a culture, not
    which compartment it belongs to, and guessing would manufacture signal.
    """
    if not isinstance(tags, str) or not tags.strip():
        return ""
    found = set()
    for t in tags.split(","):
        c = _TAG_TO_COMPARTMENT.get(t.strip().lower())
        if c:
            found.add(c)
    return ",".join(sorted(found))


# ---------------------------------------------------------------------------
# The join
# ---------------------------------------------------------------------------
def load_metadata(nuccore_path: Path,
                  taxonomy_path: Path | None = None,
                  biosample_path: Path | None = None,
                  verbose: bool = True) -> pd.DataFrame:
    """
    Join nuccore -> taxonomy + biosample and derive the cleaned fields.

    Returns a DataFrame indexed by NUCCORE_ACC with (where available):
        host_genus, host_species, host_family, host_order, host_phylum,
        country, year, ecosystem_tags, compartments, disease_group
    Missing optional tables simply mean missing columns -- nothing raises.
    """
    nc = pd.read_csv(nuccore_path, low_memory=False)
    nc["NUCCORE_ACC"] = nc["NUCCORE_ACC"].astype(str)
    keep = [c for c in ["NUCCORE_ACC", "NUCCORE_UID", "TAXONOMY_UID",
                        "BIOSAMPLE_UID", "NUCCORE_CreateDate", "NUCCORE_Length",
                        "NUCCORE_Topology", "NUCCORE_Completeness"]
            if c in nc.columns]
    m = nc[keep].copy()
    if verbose:
        print(f"  nuccore.csv          : {len(m):>7,} records", file=sys.stderr)

    if taxonomy_path and Path(taxonomy_path).exists():
        tx = pd.read_csv(taxonomy_path, low_memory=False)
        tcols = {c: c.replace("TAXONOMY_", "host_")
                 for c in ["TAXONOMY_genus", "TAXONOMY_species",
                           "TAXONOMY_family", "TAXONOMY_order",
                           "TAXONOMY_phylum"] if c in tx.columns}
        m = m.merge(tx[["TAXONOMY_UID"] + list(tcols)].rename(columns=tcols),
                    on="TAXONOMY_UID", how="left")
        if verbose:
            print(f"  taxonomy.csv         : joined "
                  f"({m['host_genus'].notna().mean() * 100:.1f}% genus)",
                  file=sys.stderr)

    if biosample_path and Path(biosample_path).exists():
        bcols = ["BIOSAMPLE_UID", "BIOSAMPLE_ACC", "ECOSYSTEM_tags",
                 "DISEASE_group", "LOCATION_name"]
        bs = pd.read_csv(biosample_path, low_memory=False)
        bs = bs[[c for c in bcols if c in bs.columns]]
        bs = bs.drop_duplicates("BIOSAMPLE_UID")
        m = m.merge(bs, on="BIOSAMPLE_UID", how="left")
        m["country"] = m["LOCATION_name"].map(country_of) \
            if "LOCATION_name" in m.columns else np.nan
        m["compartments"] = m["ECOSYSTEM_tags"].map(compartments_of) \
            if "ECOSYSTEM_tags" in m.columns else ""
        if verbose:
            print(f"  biosample.csv        : joined "
                  f"({m['country'].notna().mean() * 100:.1f}% country, "
                  f"{(m['compartments'] != '').mean() * 100:.1f}% compartment)",
                  file=sys.stderr)

    if "NUCCORE_CreateDate" in m.columns:
        m["year"] = pd.to_datetime(m["NUCCORE_CreateDate"],
                                   errors="coerce").dt.year

    m = m.drop_duplicates("NUCCORE_ACC").set_index("NUCCORE_ACC")
    return m
