"""Common utilities for tool modules."""

import os

import consts
from utils.normalise import main as _normalise


def create_csv(chrom, tool, tsv_out_name):
    """Converts a TSV output to a CSV.

    Args:
        chrom: A Chromosome instance for which the output was generated.
        tool: A Tool instance for the tool which generated the output.
        tsv_out_name: The name of the original TSV output. Should have a single
            %d for the number of the chromosome.
    """
    out_path = tool.get_out_dir(chrom.dataset.name)
    tsv_out_path = os.path.join(out_path, tsv_out_name % chrom.name)
    csv_out_path = os.path.join(out_path, consts.RAW_OUT_CSV % chrom.num)
    with open(tsv_out_path, 'r') as tsv, open(csv_out_path, 'w') as csv:
        output = tsv.read()
        csv.write(output.replace('\t', ','))


def normalise(chrom, tool):
    """Runs the output normalisation script.

    Args:
        chrom: A Chromosome instance for which the output was generated.
        tool: A Tool instance for the tool which generated the output.

    Writes:
        The normalised output.
    """
    args = [
        "-d", chrom.dataset.name,
        "-c", str(chrom.num),
        "-t", tool.name,
        "-p", chrom.work_dir
    ]
    _normalise(args)
