"""Samples targets from the GenomeCRISPR DB into CSV dataset file."""

import os
import random

import consts
import common
from utils.datasets import DATASETS


OUT_FILE = "chr%d_%s.csv"
FULL_DATASET = "GenomeCRISPR_full05112017.csv"
DATASET_NAME = "genome"
EFFICIENCY_THRESHOLD = -1

# The index of the column which specifies the screen type.
SCREEN_TYPE_COL = 15
# The index of the column which specifies the effect label.
SCORE_COL = 13
# Number of progress prints to provide.
PRINT_NUM = 10

# Reasons for rejecting targets.
WRONG_END_MINUS_START = "end-start is not the expected target length"
BAD_PAM = "PAM is not NGG"
TARGET_BAD_LEN = "The target is of the wrong length"
NOT_NEG_SCREEN = "The screening is not negative"
EXISTING_SPACER = "Spacer is already included"
# Collects the rejection statistics.
STATS = {
    WRONG_END_MINUS_START : 0,
    BAD_PAM : 0,
    TARGET_BAD_LEN: 0,
    NOT_NEG_SCREEN: 0,
    EXISTING_SPACER: 0,
}


def get_args(raw_args):
    parser = common.default_parser(
        "Samples targets from the GenomeCRISPR DB into CSV dataset file.")
    parser.set_defaults(dataset=DATASET_NAME)
    parser.add_argument(
        "-s", "--size", type=int, default=1000,
        help="Number of targets to sample")
    parser.add_argument(
        "-i", "--input", type=str, default=FULL_DATASET,
        help="Name of the full DB file")
    parser.add_argument(
        "-o", "--output", type=str,
        help=("Name of the output file. By default, the name is "
              "chr<chromosome number>_<sample size>.csv."))
    return parser.parse_args(raw_args)


def format_number(number):
    """Simplifies the number by converting 1,000 to k.

    Args:
        number: A number to simplify, if possible.

    Returns:
        If the number is at least 1,000 and is divisible by 1,000, it is
        simplified (e.g. 3,000 => 3k and 150,000 => 150k).
        Otherwise, the number is returned (unchanged) as a string.
    """
    if number >= 1000 and number % 1000 == 0:
        return "%dk" % (number/1000)
    else:
        return str(number)


def add_headers(headers, out):
    """Writes the DB headers along with an efficiency header to the output file.

    Args:
        headers: The existing headers in the GenomeCRISPR DB.
        out: A file descriptor for the output file.
    """
    out.write(common.to_csv_line(headers, "efficient"))


def is_valid(values, dataset):
    """Filters targets which are invalid.

    Args:
        values: The values in the GenomeCRISPR entry.
        dataset: A Dataset instance, used to parse the values.

    Returns:
        A bool indicating if the entry is valid.
    """
    # Only includes negative screens.
    if values[SCREEN_TYPE_COL] != "negative selection":
        STATS[NOT_NEG_SCREEN] += 1
        return False
    # Targets must have the correct length.
    if int(values[dataset.end_idx]) - int(values[dataset.start_idx]) !=\
            consts.TARGET_LEN:
        STATS[WRONG_END_MINUS_START] += 1
        return False

    target = dataset.get_target(values)
    # Targets must have an NGG PAM sequence.
    if not target.endswith("GG"):
        STATS[BAD_PAM] += 1
        return False
    # Another safety measure against targets with the wrong length.
    if len(target) != consts.TARGET_LEN:
        STATS[TARGET_BAD_LEN] += 1
        return False
    return True


def sample(entries, spacers, num, dataset):
    """Sample targets from the GenomeCRISPR DB.

    Args:
        entries: Entries to sample from.
        spacers: Spacers already considered (where a spacer is a target without
            its PAM).
        num: Number to sampling attempts to make.
        dataset: A Dataset instance used to parse entries.

    Returns:
        The chosen entries. Since entries are filtered out after sampling, it is
        possible to sample fewer entries than the number of sampling attempts.
    """
    sampled = random.sample(entries, num)
    valid_entries = []
    for entry in sampled:
        entries.remove(entry)
        values = common.from_csv_line(entry)
        # Filters bad entries.
        if not is_valid(values, dataset): continue

        spacer = dataset.get_target(values)[:-consts.PAM_LEN]
        # Avoids duplicate spacers.
        if spacer in spacers:
            STATS[EXISTING_SPACER] += 1
            continue
        spacers.add(spacer)

        # Adds efficiency label and saves entry.
        efficient = 1 if int(values[SCORE_COL]) <= EFFICIENCY_THRESHOLD else 0
        valid_entries.append(values + [efficient])

    return valid_entries


def main(raw_args=None):
    """Samples targets from the GenomeCRISPR DB into a CSV dataset file.

    For command line help, run with the '-h' flag.

    Writes:
        A dataset CSV file with the sampled targets.
    """
    args = get_args(raw_args)
    dataset = DATASETS[args.dataset]
    dataset.set_work_dir(args.path)

    dataset_dir = os.path.join(dataset.get_data_path(), consts.RAW_DATA_DIR)
    input_file = os.path.join(dataset_dir, args.input)
    if args.output:
        out_name = args.output
    else:
        out_name = OUT_FILE % (args.chr, format_number(args.size))
    output_file = os.path.join(dataset_dir, out_name)

    print("Output will be in: %s" % output_file)

    entries = set()
    out = open(output_file, 'w')
    input = open(input_file, 'r')

    is_first = True
    for line in input:
        if is_first:
            add_headers(line.strip(), out)
            is_first = False
        else:
            entries.add(line)
    input.close()
    print("Done reading")

    # Keeps track of spacers that were already included.
    spacers = set()
    # Keeps track of the number of samples already obtained.
    sampled = 0
    # Keeps track of progress for printing.
    chunk = 1

    # While we have not sampled enough...
    while sampled < args.size:
        remaining = args.size - sampled
        # Sample more entries.
        new_entries = sample(entries, spacers, remaining, dataset)
        sampled += len(new_entries)
        for entry in new_entries:
            out.write(common.to_csv_line(*entry))

        # If we have passed a progress checkpoint - prints progress.
        if (args.size / PRINT_NUM) * chunk < sampled:
            print("Sampled: %d" % sampled)
            chunk += 1

    print("Done")
    print(STATS)


if __name__ == "__main__":
    main()
