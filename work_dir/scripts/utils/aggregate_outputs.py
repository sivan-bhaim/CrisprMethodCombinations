"""
Aggregates the outputs of all the tools for a specific dataset.
"""

from os.path import join as join_path
from collections import defaultdict

import consts
import common
from tools import TOOLS

DEFAULT_SCORE = "x"  # Used as a filler for missing scores.


def get_args():
    parser = common.default_parser(
        "Aggregates the outputs of all the tools for a specific dataset.")
    return parser.parse_args()


def get_spacer(target_sequence):
    spacer_len = consts.TARGET_LEN - consts.PAM_LEN
    return target_sequence[:spacer_len]


def add_scores(args, tool, scored_targets):
    """Adds the scores from the given tool.

    Args:
        args: The arguments given to the script.
        tool: A Tool instance for the tool to read the scores from.
        scored_targets: A dictionary mapping target sequences to dictionaries of
            scores.
    """
    # Reads the output of the tool.
    data_path = tool.get_normalised_out_path(args.chr, args.dataset)
    with open(data_path, 'r') as fd:
        data = fd.readlines()

    for line in data:
        values = common.from_csv_line(line)
        target = get_spacer(values[consts.TARGET_COLUMN])
        if target not in scored_targets: continue
        target_scores = values[consts.SCORE_COLUMN:]
        # The length and order of target scores and the number of scores for the
        # tool should match
        for i, score_name in enumerate(tool.scores.values()):
            scored_targets[target][tool.name][score_name] = target_scores[i]


def get_headers():
    """Returns the headers for the output CSV file."""
    headers = ["target"]
    for tool in TOOLS.values():
        for score in tool.scores.values():
            headers.append("%s_%s" % (tool.name, score))
    return headers


def aggregate(targets, scored_targets):
    """Aggregates the scores.

    Args:
        targets: A list of target sequences.
        scored_targets: A dictionary mapping target sequences to dictionaries of
            scores.

    Returns:
        A list where each item is a record: the record always starts with the
        target sequence, followed by all the scores for the target.
    """
    records = []
    for target in targets:
        seq = get_spacer(target)
        record = [target]
        for tool in TOOLS.values():
            for score_name in tool.scores.values():
                score = scored_targets[seq][tool.name].get(score_name,
                                                           DEFAULT_SCORE)
                record.append(score)
        records.append(record)

    return records


def main():
    """Aggregates the outputs of all the tools for a specific dataset.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset, chrom, _ = common.process_args(args)
    if not dataset: return

    # Collects the scores to the scored_targets dictionary, where the first ket
    # is the target spacer, the second is the tool and the third is the score's
    # name.
    targets = dataset.get_targets(args.chr)
    scored_targets = {}
    for target in targets:
        scored_targets[get_spacer(target)] = defaultdict(defaultdict)
    for tool in TOOLS.values():
        add_scores(args, tool, scored_targets)

    # Converts the mapping to a list of records, each one containing all of the
    # scores of its corresponding target.
    records = aggregate(targets, scored_targets)
    # Writes the output.
    aggregate_path = join_path(dataset.get_out_path(), consts.AGG_OUT_NAME % args.chr)
    with open(aggregate_path, 'w') as out:
        out.write(common.to_csv_line(*get_headers()))
        for record in records:
            out.write(common.to_csv_line(*record))


if __name__ == "__main__":
    main()
