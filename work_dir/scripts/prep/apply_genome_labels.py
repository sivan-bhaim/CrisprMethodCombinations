"""Adds the GenomeCRISPR effect labels to the dataset."""

import numpy as np
from collections import OrderedDict

import common


EFFICIENCY_THRESHOLD = -1


def get_args():
    parser = common.default_parser(
        "Adds the GenomeCRISPR effect labels to the dataset.")
    parser.add_argument(
        "--positive", action="store_true",
        help=("Use this option if the activity metric is non-negative (like "
              "Chari)."))
    parser.add_argument(
        "--shift", type=float, default=0,
        help=("Shifts all the values before setting percentiles. If used along "
              "with --invert, the shift occurs before the inversion."))
    parser.add_argument(
        "--invert", action="store_true",
        help=("Inverts all the labels. Use this when the activity score has "
              "both positive and negative values, and positive scores "
              "represent high efficiency (e.g. the gene rank for adjusted "
              "Doench)."))
    parser.add_argument(
        "--labels", type=int, nargs='+', default=list(range(10)),
        help=("A list of positive labels to use, from low to high. The "
              "negative labels will be the same, only inverted. Given as a "
              "sequence of space separated values, for example: "
              "--labels 2 3 4 5."))
    return parser.parse_args()


def get_label(activity, pos_thresholds, neg_thresholds, pos_only, labels):
    """Returns the label which the given activity belongs to.

    Args:
        activity: An activity score.
        pos_thresholds: The minimum thresholds for positive labels.
        neg_thresholds: The maximum thresholds for negative labels.
        pos_only: A bool indicating if the activity is always non-negative.
        labels: The positive labels available, from low to high. The negative
            labels are assumed to be the inversions of the positive labels.

    Returns:
        The label of the activity score.
    """
    if activity >= 0:
        # Looks from the highest threshold and down.
        for i in range(len(pos_thresholds)-1, 0, -1):
            if activity >= pos_thresholds[i]:
                label = labels[i]
                # When the activity is always non-negative, only negative labels
                # are used. For example, the Chari dataset has mutation rates
                # that are always non-negative. But any mutation rate is a
                # measure of "good" activity, and the "good" labels are the
                # negative labels, since we deal with negative screens. Hence,
                # in this case, the label is inverted.
                if pos_only: label = -labels[i]
                return label
        return labels[0]

    # Look from the lowest threshold and up.
    for i in range(len(neg_thresholds) - 1, 0, -1):
        if activity <= neg_thresholds[i]:
            return -labels[i]
    return -labels[0]


def get_pos_thresholds(vals, num_labels):
    """Produces a list of minimum thresholds for positive activity scores.

    Args:
        vals: The activity scores.
        num_labels: The number of labels (which is also the number of
            thresholds).

    Returns:
        The list of positive thresholds.
    """
    pos_vals = [val for val in vals if val >= 0]
    quantile = 100 // num_labels
    return [np.percentile(pos_vals, quantile*i) for i in range(num_labels)]


def get_neg_thresholds(vals, num_labels):
    """Produces a list of maximum thresholds for negative activity scores.

    Args:
        vals: The activity scores.
        num_labels: The number of labels (which is also the number of
            thresholds).

    Returns:
        The list of negative thresholds.
    """
    neg_vals = [-val for val in vals if val < 0]
    quantile = 100 // num_labels
    return [-np.percentile(neg_vals, quantile*i) for i in range(num_labels)]


def shift_values(activity, shift):
    """Shifts all the activity scores by a fixed number.

    Args:
        activity: A map from targets to their activity scores.
        shift: The amount to shift by.

    Returns:
        The shifted activity mapping.
    """
    if not shift: return activity
    for target, score in activity.iteritems():
        activity[target] = score + shift
    return activity


def main():
    """Adds the GenomeCRISPR effect labels to the dataset.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset, chrom, _ = common.process_args(args)
    if not dataset: return

    targets = dataset.get_targets(chrom.num)
    activity = OrderedDict()
    for target in targets:
        activity[target] = float(
            dataset.get_value(chrom.num, target, dataset.raw_activity_idx))

    num_labels = len(args.labels)
    shifted_activity = shift_values(activity, args.shift)
    pos_thresholds = get_pos_thresholds(shifted_activity.values(), num_labels)
    neg_thresholds = None if args.positive\
        else get_neg_thresholds(shifted_activity.values(), num_labels)

    labels = OrderedDict()
    for target, freq in shifted_activity.iteritems():
        label = get_label(
            freq, pos_thresholds, neg_thresholds, args.positive, args.labels)
        if args.invert: label = -label
        efficient = 1 if label < EFFICIENCY_THRESHOLD else 0
        labels[target] = [label, efficient]

    chrom.add_values(["GenomeCRISPRlabel", "genome_efficient"], labels)


if __name__ == "__main__":
    main()
