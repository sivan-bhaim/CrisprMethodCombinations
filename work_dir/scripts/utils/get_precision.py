"""Get precision stats for a tool against experimental data."""

import consts
import common


EPSILON = 1e-10


def get_args():
    parser = common.default_parser(
        "Get precision stats for a tool according to experimental data.")
    parser.add_argument(
        "-t", "--tool", type=str,
        help="Prefix of one of the tools' names.")
    parser.add_argument(
        "-s", "--score", type=str, default=None,
        help="Name of the score to use (e.g. g20, doench_2014, decision, etc.)")
    parser.add_argument(
        "--no_pam", action="store_true",
        help="Compare without PAM (should be used for mm10db)")
    parser.add_argument(
        "--relative", action="store_true",
        help="Get precision based on relative scores")
    parser.add_argument(
        "--threshold", type=float, default=0,
        help="A threshold score above which a target is considered accepted")
    parser.add_argument(
        "--non_negative", action="store_true",
        help="Targets are accepted if and only if they have a non-negative "
             "score. Equivalent to setting the threshold at some negative "
             "epsilon.")

    return parser.parse_args()


def print_stat(description, stat):
    print "%s: %s" % (description.ljust(40), stat)


def get_overlap(found, experimental):
    """Returns the intersection between found and experimental targets.

    Args:
        found: A dictionary mapping targets to their scores.
        experimental: A set of experimental targets.

    Returns:
        A set with the targets in the intersection of the found and experimental
        targets.
    """
    found_set = set(found.keys())
    return found_set.intersection(experimental)


def abs_precision(found, efficient, inefficient, threshold):
    """Computes absolute precision.

    The absolute precision is the percentage of experimentally efficient guides
    out of all the guides accepted by the tool.

    Args:
        found: A dictionary mapping targets to their scores according to the
            tool.
        efficient: A set of efficient targets.
        inefficient: A set of inefficient targets.
        threshold: A threshold score above which the guide is accepted by the
            tool.

    Returns:
        The absolute precision.
    """
    rejected, accepted = 0, 0
    accepted_exp, accepted_eff = 0, 0

    for target, score in found.iteritems():
        if score <= threshold:
            rejected += 1
            continue
        accepted += 1
        if target in efficient:
            accepted_eff += 1
            accepted_exp += 1
        if target in inefficient:
            accepted_exp += 1

    print_stat("Number of targets accepted", accepted)
    print_stat("Number of targets rejected", rejected)
    print_stat("Number of efficient targets accepted", accepted_eff)
    print_stat("Number of experimental targets accepted", accepted_exp)

    return 100 * float(accepted_eff) / accepted_exp


def relative_precision(found, efficient, inefficient):
    """Computes relative precision (as needed for sgRNA scorer 2.0).

    Consider all the pairs of guides where one is efficient and the other
    inefficient. The tool is correct about a pair if it gives the efficient
    guide a higher score than the inefficient one. The relative precision is
    the percentage of these pairs about which the tool is correct.

    Args:
        found: A dictionary mapping targets to their scores according to the
            tool.
        efficient: A set of efficient targets.
        inefficient: A set of inefficient targets.

    Returns:
        The relative precision of the tool.
    """
    found_efficient = get_overlap(found, efficient)
    print_stat("Number of efficient targets found", len(found_efficient))

    total, correct = 0, 0
    for eff_target in efficient:
        if eff_target not in found: continue
        for ineff_target in inefficient:
            if ineff_target not in found: continue
            total += 1
            if found[eff_target] >= found[ineff_target]:
                correct += 1
    return 100 * float(correct) / total


def main():
    """Gets precision stats for a tool against experimental data.

    For command line help, run with the '-h' flag.

    Prints:
        Precision statistics for the tool.
    """
    args = get_args()
    dataset, chrom, tool = common.process_args(args, tool_name=args.tool)
    if not dataset or not tool: return

    # If a specific scoring signal is wanted rather than just the first one,
    # adjust the column.
    score_start_idx = consts.SCORE_COLUMN
    score_shift_idx = 0
    if args.score:
        score_shift_idx = tool.get_score_idx(args.score)
    if score_shift_idx == -1:
        print "unknown score %s for tool %s" % (args.score, args.tool)
        return
    score_idx = score_start_idx + score_shift_idx

    # Reads the normalised output of the tool.
    found = {}
    output_path = tool.get_normalised_out_path(chrom.num, dataset.name)
    with open(output_path, 'r') as fd:
        for line in fd:
            values = line.split(',')
            target = values[consts.TARGET_COLUMN]
            score = values[score_idx] if len(values) > score_idx else 0
            if args.no_pam: target = target[:-consts.PAM_LEN]
            found[target] = float(score)
    print_stat("Number of targets found", len(found))

    efficient = set(dataset.get_efficient_targets(args.chr, args.no_pam))
    inefficient = set(dataset.get_inefficient_targets(args.chr, args.no_pam))

    if args.relative:
        precision = relative_precision(found, efficient, inefficient)
    else:
        threshold = -EPSILON if args.non_negative else args.threshold
        precision = abs_precision(found, efficient, inefficient, threshold)
    print_stat("\nPrecision", "%.02f%%" % precision)


if __name__ == '__main__':
    main()