"""Produces a features file for models to be trained on.

The first value in each line is the target sequence. Next is the GenomeCRISPR
effect label. Then follow all the scores and signals.
"""

from os.path import join as join_path
from itertools import chain

import consts
import common
from tool_utils import run_mm10db as mm10db
from aggregate_outputs import DEFAULT_SCORE


# Marks missing scores.
UNHANDLED_MARKER = DEFAULT_SCORE
# A score to use for a guide that was not considered by the tool (except for
# mm10db). We operate under the assumption that the scores are (roughly) between
# -1 and 1.
UNHANDLED_SCORE = -3

OUT_FILE = "%s%d_features.csv"

# The list of possible mm10db decisions.
MM10DB_DECISIONS = [
    mm10db.REJECTIONS[mm10db.AT],               # 0
    mm10db.REJECTIONS[mm10db.MULTI_MATCH],      # 1
    mm10db.REJECTIONS[mm10db.STRUCT_ENERGY],    # 2
    mm10db.REJECTIONS[mm10db.REVERSE_PRIMER],   # 3
    mm10db.REJECTIONS[mm10db.TTTT],             # 4
    mm10db.REJECTIONS[mm10db.OFF_TARGET],       # 5
    1,                                          # 6 (accepted)
]


def decisions_to_encoding(decisions):
    """Converts a list of possible decisions to a one-hot encoding mapping.

    Args:
        decisions: A list of possible decisions.

    Returns:
        A mapping from decision to its one-hot encoding.
    """
    encodings = {}
    for idx, val in enumerate(decisions):
        vector = [0] * len(decisions)
        vector[idx] = 1
        encodings[val] = vector
    return encodings


def encode_score(score):
    """Encodes a standard score by filling in the missing scores."""
    if score == UNHANDLED_MARKER:
        return [UNHANDLED_SCORE]
    return [score]

def encode_mm10db(score):
    """Encodes an mm10db decision using one-hot encoding."""
    encodings = decisions_to_encoding(MM10DB_DECISIONS)
    # We attempted to minimise off-target considerations, but some effect
    # remains, so some targets are rejected for this reason. This is equivalent
    # to not having any efficiency prediction. Thus, we treat unhandled targets
    # in the same way we treat those rejected due to off-target considerations.
    if score == UNHANDLED_MARKER:
        return encodings[mm10db.REJECTIONS[mm10db.OFF_TARGET]]
    return encodings[float(score)]


def get_column_handlers():
    """Returns a list of encoding functions matching the columns of scores."""
    handlers = \
        [encode_score for _ in range(8)] +\
        [encode_mm10db] +\
        [encode_score for _ in range(3)]
    return handlers


def parse_features(path, handlers, fd_out, label_getter):
    """Parses the input to produce the feature representations.

    Args:
        path: Path to the input.
        handlers: A list of encoding functions for the columns of the input.
        fd_out: A file descriptor for the output file.
        label_getter: A function which gets the GenomeCRISPR effect label for a
            given target.
    """
    with open(path, 'r') as fd:
        lines = fd.readlines()
        for line in lines[1:]:
            values = common.from_csv_line(line)
            target = values[0]
            encoded = list(chain(*[handlers[i](values[i+1])\
                                   for i in range(len(handlers))]))
            entry = [target, label_getter(target)] + encoded
            fd_out.write(common.to_csv_line(*entry))


def get_args():
    parser = common.default_parser(
        "Produces a features file for models to be trained on.")
    return parser.parse_args()


def main():
    """Produces a features file for models to be trained on.

    For command line help, run with the '-h' flag.

    Writes:
        An output CSV file with the targets, labels and feature representations.
    """
    args = get_args()
    dataset, chrom, _ = common.process_args(args)
    if not dataset: return

    out_file = join_path(
        args.path, consts.DATA_DIR, OUT_FILE % (dataset.name, chrom.num))
    fd_out = open(out_file, 'w')
    handlers = get_column_handlers()

    out_dir = dataset.get_out_path()
    in_file = join_path(out_dir, consts.AGG_OUT_NAME % chrom.num)
    def label_getter(target):
        return dataset.get_value(args.chr, target, dataset.genome_label_idx)
    parse_features(in_file, handlers, fd_out, label_getter)

    fd_out.close()


if __name__ == "__main__":
    main()
