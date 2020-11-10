"""The main configuration file."""

import os

WORK_DIR = os.environ['WORK_DIR']
DATA_DIR = "data"
TOOLS_DIR = "tools"
OUTPUT_DIR = "outputs"
INDEXES_DIR = "indexes"
RAW_DATA_DIR = "dataset"

TARGET_LEN = 23
PAM_LEN = 3
# Used to separate targets in the artificial genome.
NUM_N_BETWEEN_GUIDES = 50

# The format of the normalised output files.
TARGET_COLUMN = 1
SCORE_COLUMN = 5

# Standard file names.
REFGENE_FNAME = "chr%d_refGene.txt"
NORMALISED_OUT = "chr%d_targets.normalised"
RAW_OUT_CSV = "chr%d_out.csv"
RAW_DATA_CSV = "chr%d_%s.csv"
AGG_OUT_NAME = "chr%d_agg_results.csv"
