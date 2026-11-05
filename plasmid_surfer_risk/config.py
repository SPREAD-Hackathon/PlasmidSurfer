"""
Central configuration: file paths and tunable constants.

EDIT THE PATHS AT THE TOP to point at your copies of the PLSDB tables.
Everything else in the pipeline reads from here, so this is the only file
you need to touch when the data moves.
"""

from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------------
# INPUT FILES
# --------------------------------------------------------------------------
# Override any of these with an environment variable of the same name, e.g.
#   export PLRISK_TYPING=/data/plsdb/typing.csv
# so nobody has to edit the file on a shared machine.

DATA_DIR = Path(os.environ.get("PLRISK_DATA", "data"))

TYPING      = Path(os.environ.get("PLRISK_TYPING",      DATA_DIR / "typing.csv"))
AMR         = Path(os.environ.get("PLRISK_AMR",         DATA_DIR / "amr.tsv"))
PLASMIDFIND = Path(os.environ.get("PLRISK_PLASMIDFIND", DATA_DIR / "plasmidfinder.csv"))
MASH_SIM    = Path(os.environ.get("PLRISK_MASH",        DATA_DIR / "plsdb_mashdb_sim.tsv"))

# OPTIONAL — the pipeline runs without these, but they make the label better.
# NUCCORE is the PLSDB master table (NUCCORE_ACC, NUCCORE_UID, TAXONOMY_UID,
# BIOSAMPLE_UID, NUCCORE_CreateDate...). Without it you cannot join plasmid
# accessions to biosample metadata -- see README, "The missing table".
NUCCORE   = Path(os.environ.get("PLRISK_NUCCORE",   DATA_DIR / "nuccore.csv"))
BIOSAMPLE = Path(os.environ.get("PLRISK_BIOSAMPLE", DATA_DIR / "biosample.csv"))
TAXONOMY  = Path(os.environ.get("PLRISK_TAXONOMY",  DATA_DIR / "taxonomy.csv"))

# A PLING community assignment, if/when your teammate's run finishes.
# Expected format: two columns, plasmid id and community id, any delimiter,
# with or without a header. If present it REPLACES the Mash components as the
# grouping. See README, "Swapping in PLING".
PLING = Path(os.environ.get("PLRISK_PLING", DATA_DIR / "pling_communities.tsv"))

OUT_DIR = Path(os.environ.get("PLRISK_OUT", "out"))

# --------------------------------------------------------------------------
# MODELLING CONSTANTS
# --------------------------------------------------------------------------

# A group (Mash component / MOB cluster / PLING community) needs at least this
# many members before its "dissemination success" is treated as a real
# observation rather than noise. Groups below this are kept for scoring but
# excluded from TRAINING, because a label computed from 1 genome is mostly
# sampling accident.
#
# Keep this LOW. Raising it looks like quality control but is actually
# range restriction: group size is one third of the label, so filtering on
# group size filters on the target. At >=3 the training label distribution
# collapses into the top quartile and the model has almost nothing to learn.
MIN_GROUP_SIZE_FOR_TRAINING = 2

# At most this many members of any one lineage enter the training set.
# PLSDB contains lineages with 1,700 near-identical records; without a cap,
# three epidemic plasmids would supply most of the training rows and the model
# would simply memorise them. Combined with the 1/group_size sample weights
# used in train.py, this makes every lineage count roughly equally.
MAX_MEMBERS_PER_GROUP = 5

# Keep only categorical levels seen at least this many times, to stop the
# feature matrix exploding into thousands of near-empty one-hot columns.
MIN_CATEGORY_COUNT = 30

# Cross-validation folds (grouped -- never split a group across folds).
N_FOLDS = 5

RANDOM_STATE = 42

# How the two components combine into the final 0-100 risk score.
# risk = 100 * dissemination * (CARGO_FLOOR + (1 - CARGO_FLOOR) * cargo)
# CARGO_FLOOR > 0 means a highly mobile plasmid with no resistance cargo still
# scores above zero -- it is an empty vehicle, not a harmless one. Defensible
# either way; say which you chose and why in the talk.
CARGO_FLOOR = 0.25
