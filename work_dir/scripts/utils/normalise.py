"""
Normalises the output of the prediction tools.

This method of unifying the outputs of multiple tools is based on code provided
in the supplementary information of Bradford and Perrin (2019).
"""

import common


def get_args(raw_args):
    parser = common.default_parser(
        "Normalise the outputs of the supported prediction tools.")
    parser.add_argument(
        "-t", "--tool", default=None,
        help="Prefix of one of the tools' names.")
    parser.add_argument(
        "-i", "--input", default=None,
        help="Path of input file (if not provided, the default input CSV of "
             "the dataset will be used)")
    parser.add_argument(
        "-o", "--output", default=None,
        help="Path of output file (if not provided, the default normalised "
             "output file will be used)")
    return parser.parse_args(raw_args)


def main(raw_args=None):
    """Normalises the output of the various prediction tools.

    For command line help, run with the '-h' flag.
    """
    args = get_args(raw_args)
    dataset, chrom, tool = common.process_args(args, tool_name=args.tool)
    if not dataset or not tool: return

    input_path = args.input if args.input else\
        tool.get_csv_out_path(args.chr, args.dataset)
    output_path = args.output if args.output else\
        tool.get_normalised_out_path(args.chr, args.dataset)

    # Reads the unnormalised output.
    first_entry = 1 if tool.has_headers else 0
    with open(input_path, 'r') as in_fd:
        entries = [line.strip() for line in in_fd.readlines()[first_entry:]]

    # Writes the normalised output.
    with open(output_path, 'w') as out_fd:
        normalised = tool.normalise(entries)
        for record in normalised:
            out_fd.write(record)


if __name__ == "__main__":
    main()
