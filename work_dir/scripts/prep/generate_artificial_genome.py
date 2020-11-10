"""
Creates an artificial genome from a dataset file.

This method of producing a sequence of desired targets is based on code provided
in the supplementary information of Bradford and Perrin (2019).
"""

import consts
import common
from utils.datasets import DATASETS


# The default length of lines in the fasta file.
LINE_LEN = 50
# A template for the header of the fasta file.
HEADER = ">chr%d\n"


def get_args():
    parser = common.default_parser(
        "Creates an artificial genome from a dataset file.")
    parser.add_argument(
        "-l", "--line_len", type=int, default=LINE_LEN,
        help="Length of lines in the output fasta file")
    parser.add_argument(
        "-t", "--target_len", type=int, default=consts.TARGET_LEN,
        help="Length of target sequences")
    return parser.parse_args()


def break_lines(genome, line_len):
    """Breaks the sequence into short lines.

    Args:
        genome: The genome sequence.
        line_len: The length of the lines.

    Returns:
        A list of strings, each of length line_len (except for the last which
        could be shorter), that concatenated together form the genome.
    """
    return [genome[i:i+line_len] for i in range(0, len(genome), line_len)]


def main():
    """Creates an artificial genome from a dataset file.

    For command line help, run with the '-h' flag.

    Writes:
        An artificial genome fasta file along with a matching .txt file with the
        entire sequence in one line.
    """
    args = get_args()
    dataset = DATASETS[args.dataset]
    dataset.set_work_dir(args.path)
    chrom = dataset.get_chr(args.chr)

    targets = dataset.get_targets(args.chr, args.target_len)
    guide_gap = 'N' * consts.NUM_N_BETWEEN_GUIDES
    # Starts the sequence with padding, since some tools struggle with targets
    # at the beginning of the sequence.
    genome = guide_gap
    for target in targets:
        genome += target
        genome += guide_gap

    lines = break_lines(genome, args.line_len)
    # Writes the fasta file.
    with open(chrom.get_path(), 'w') as fd_fasta:
        fd_fasta.write(HEADER % args.chr)
        [fd_fasta.write(line + '\n') for line in lines]
    # Writes a txt file, with the entire sequence in one line.
    # The genome could be too long for writing in a single chunk, hence we use
    # the broken representation.
    with open(chrom.get_path('txt'), 'w') as fd_txt:
        [fd_txt.write(line) for line in lines]


if __name__ == '__main__':
    main()
