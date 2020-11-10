"""Produces a fake refGene file for an artificial genome."""

import re

import consts
import common
from utils.datasets import DATASETS

# Identifies a target sequence.
TARGET_PATTERN = '([ACGT]{%d})[N]{%d}'
# A template for an entry in a refGene file.
REFGENE_LINE = (
    "{idx}\t{gene_name1}\t{chr_name}\t{strand}\t"
    "{start}\t{end}\t{start}\t{end}\t1\t{start},\t{end},\t"
    "0\t{gene_name2}\tcmpl\tcmpl\t0\n"
)
# Templates for gene names.
GENE_NAME1 = "NM_%d%08d"
GENE_NAME2 = "NM_%d%05d"


def get_args(raw_args):
    parser = common.default_parser(
        "Produces a fake refGene file for an artificial genome.")
    parser.add_argument(
        "-r", "--refgene", type=str, default=consts.REFGENE_FNAME,
        help="Name of output refGene file, with a single %%d for chr number")
    parser.add_argument(
        "-t", "--target_len", type=int, default=consts.TARGET_LEN,
        help="Length of target sequences")

    return parser.parse_args(raw_args)


def write_refgene(chrom, target_len, fd_out):
    """Writes the contents of the refGene file.

    Args:
        chrom: A Chromosome instance.
        target_len: The length of target sequences.
        fd_out: A file descriptor for the output file.
    """
    chr_name = chrom.get_name()
    # Reads the chromosome sequence.
    with open(chrom.get_path('txt'), 'r') as chr_file:
        genome = chr_file.read()

    # Identifies all the target sequences in the chromosome.
    matcher = re.compile(
        TARGET_PATTERN % (target_len, consts.NUM_N_BETWEEN_GUIDES))

    # Adds an entry for each target.
    for i, m in enumerate(matcher.finditer(genome)):
        idx = i+1
        gene_name1 = GENE_NAME1 % (chrom.num, idx)
        gene_name2 = GENE_NAME2 % (chrom.num, idx)
        strand = '+'
        start = m.start()
        end = start + target_len
        line = REFGENE_LINE.format(**locals())
        fd_out.write(line)


def main(raw_args=None):
    """Produces a fake refGene file for an artificial genome.

    For command line help, run with the '-h' flag.

    Writes:
        A refGene file.
    """
    args = get_args(raw_args)
    dataset = DATASETS[args.dataset]
    dataset.set_work_dir(args.path)
    chrom = dataset.get_chr(args.chr)

    out_path = chrom.get_refgene(args.refgene)
    with open(out_path, 'w') as fd_out:
        write_refgene(chrom, args.target_len, fd_out)


if __name__ == '__main__':
    main()
