"""Common utilities."""

import argparse
from os.path import join as join_path

import consts
from utils.tools import TOOLS
from utils.datasets import DATASETS


def default_parser(description, multiple_chrs=False):
    """Creates an argparse ArgumentParser with the standard arguments:
    - Chromosome number(s) with a default of 1 for the single chromosome case
    - Path to the working directory with a default of consts.WORK_DIR
    - Dataset name chosen from the keys of utils.datasets.DATASETS

    Args:
        description: A description for the parser.
        multiple_chrs: A bool indicating whether multiple chromosome numbers
            should be accepted.

    Returns:
        The parser.
    """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    if multiple_chrs:
        parser.add_argument(
            "-c", "--chrs", type=int, nargs='+',
            help="The numbers of the chromosomes")
    else:
        parser.add_argument(
            "-c", "--chr", type=int, default=1,
            help="The number of the chromosome")
    parser.add_argument(
        "-p", "--path", type=str, default=consts.WORK_DIR,
        help="Path to the main work directory")
    parser.add_argument(
        "-d", "--dataset", type=str, choices=DATASETS.keys(),
        help="Name of the dataset")
    return parser


def process_args(args, multiple_chrs=False, tool_name=""):
    """Processes parsed argparse arguments to retrieve standard objects.

    Args:
        args: The parsed arguments (an argparse.Namespace).
        multiple_chrs: A bool indicating whether multiple chromosome numbers
            were allowed.
        tool_name: A prefix of a tool name to retrieve. If none is provided, no
            tool will be retrieved (i.e. the tool will be None).

    Returns:
        The Dataset and Chromosome(s) instances. If tool_name was specified,
        also returns a Tool instance. If the dataset or tool were not
        identified, None will be returned for them.
    """
    # Retrieves the dataset.
    dataset, chrom, tool = None, None, None

    dataset = DATASETS.get(args.dataset, None)
    if not dataset:
        print "Unknown dataset %s." % args.dataset
        return dataset, chrom, tool
    dataset.set_work_dir(args.path)

    # Retreieves the Chromosome(s).
    if multiple_chrs:
        chrom = [dataset.get_chr(chr_num) for chr_num in args.chrs]
    else:
        chrom = dataset.get_chr(args.chr)

    # Retrieves the tool.
    if tool_name:
        full_name = get_tool(tool_name)
        if not full_name:
            return dataset, chrom, tool
        tool = TOOLS[full_name]
        tool.set_work_dir(args.path)

    return dataset, chrom, tool


def dataset_path(dataset, work_dir=consts.WORK_DIR):
    """Returns the path to the main directory of a dataset.

    Args:
        dataset: The name of the dataset.
        work_dir: A path to the main work directory.

    Returns:
        The directory path.
    """
    return join_path(work_dir, consts.DATA_DIR, dataset)


def out_path(dataset, work_dir=consts.WORK_DIR):
    """Returns the path to the outputs directory of a dataset.

    Args:
        dataset: The name of the dataset.
        work_dir: A path to the main work directory.

    Returns:
        The directory path.
    """
    return join_path(dataset_path(dataset, work_dir), consts.OUTPUT_DIR)


def indexes_path(dataset, work_dir=consts.WORK_DIR):
    """Returns the path to the indices directory of a dataset.

    Args:
        dataset: The name of the dataset.
        work_dir: A path to the main work directory.

    Returns:
        The directory path.
    """
    return join_path(dataset_path(dataset, work_dir), consts.INDEXES_DIR)


def get_tools_path(work_dir=consts.WORK_DIR):
    """Returns the path to the directory where all the tools are located.

    Args:
        work_dir: A path to the main work directory.

    Returns:
        The directory path.
    """
    return join_path(work_dir, consts.TOOLS_DIR)

def to_csv_line(*args):
    """Joins the args to a single CSV line (with a line break at the end)."""
    return ','.join([str(value) for value in args]) + '\n'

def from_csv_line(line):
    """Extracts the list of values from a CSV line."""
    return line.strip().split(',')

def get_tool(prefix):
    """Retrieves a Tool instance based on a prefix to the tool's name."""
    candidates = [tool for tool in TOOLS if tool.startswith(prefix.lower())]
    if not candidates:
        print("No tools match prefix %s. Available tools are: %s" % (prefix, TOOLS.keys()))
        return None
    if len(candidates) > 1:
        print("Ambiguous prefix %s (matches %s)" % (prefix, candidates))
        return None
    return candidates[0]
