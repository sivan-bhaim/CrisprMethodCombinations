"""Returns the length of a chromosome sequence."""

import common


def get_args():
    parser = common.default_parser(
        description="Returns the length of a chromosome sequence.")
    return parser.parse_args()


def main():
    """Returns the length of a chromosome sequence.

    For command line help, run with the '-h' flag.

    Prints:
        The length of the sequence.
    """
    args = get_args()
    dataset, chrom, _ = common.process_args(args)
    if not dataset: return

    length = 0
    with open(chrom.get_path(), 'r') as fd:
        for line in fd:
            if line.startswith('>'): continue
            length += len(line.strip())
    print length


if __name__ == '__main__':
    main()
