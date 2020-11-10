"""Prepares the files and directories for a new dataset."""

import os
from subprocess import Popen, PIPE

import common
from fake_refgene import main as refgene
from utils.tools import TOOLS
from utils.datasets import DATASETS


# Configures the indices which will be generated.
BWA_EXTS = [".amb", ".ann", ".bwt", ".pac", ".sa"]
BOWTIE_VERSIONS = [1, 2]
INDEXES = ["bwa"] + ["bowtie%d" % version for version in BOWTIE_VERSIONS]


def get_args():
    parser = common.default_parser(
        "Prepares the files and directories for a new dataset.", True)
    return parser.parse_args()


def make_dir(path):
    """Creates a directory if it does not already exist."""
    if not os.path.exists(path):
        os.makedirs(path)


def make_dirs(dataset):
    """Creates all the required directories for the dataset."""
    dirs =  [os.path.join(dataset.get_out_path(), tool) for tool in TOOLS] +\
            [os.path.join(dataset.get_indexes_path(), idx) for idx in INDEXES]
    [make_dir(d) for d in dirs]


def make_refgenes(dataset, chrs):
    """Creates the refGene files for the artificial genomes of the dataset.

    Args:
        dataset: A Dataset instance of the new dataset.
        chrs: A list of chromosome numbers included in the dataset.

    Writes:
        All the refGene files.
    """
    for chr_num in chrs:
        refgene_args = [
            "-c", str(chr_num),
            "-p", dataset.work_dir,
            "-d", dataset.name,
        ]
        refgene(refgene_args)


def make_bowtie_indices(dataset, chrom):
    """Creates the Bowtie indices for the artificial genome of the chromosome.

    Args:
        dataset: A Dataset instance of the new dataset.
        chrom: A Chromosome instance for which to create the indices.

    Writes:
        All the Bowtie indices for the chromosome.
    """
    index_name = "%s_chr%d" % (dataset.name, chrom.num)

    for bowtie_version in BOWTIE_VERSIONS:
        index_type = "bowtie%d" % bowtie_version
        bowtie_cmd = [
            "bowtie%d-build" % bowtie_version,
            "-f", chrom.get_path(),
            os.path.join(dataset.get_indexes_path(), index_type, index_name)
        ]
        _ = Popen(bowtie_cmd).wait()


def make_bwa_index(dataset, chrom):
    """Creates the BWA index for the artificial genome of the chromosome.

        Args:
            dataset: A Dataset instance of the new dataset.
            chrom: A Chromosome instance for which to create the index.

        Writes:
            The BWA index for the chromosome.
        """
    bwa_cmd = ["bwa", "index", chrom.get_path()]
    _ = Popen(bwa_cmd).wait()
    for ext in BWA_EXTS:
        file_path = chrom.get_path() + ext
        file_nme = os.path.basename(file_path)
        os.rename(file_path,
                  os.path.join(dataset.get_indexes_path(), "bwa", file_nme))


def make_indices(dataset, chrs):
    """Creates the indices for the artificial genomes of the dataset.

    Args:
        dataset: A Dataset instance of the new dataset.
        chrs: A list of chromosome numbers included in the dataset.

    Writes:
        All the indices.
    """
    for chr_num in chrs:
        chrom = dataset.get_chr(chr_num)
        # Creates bowtie indices.
        make_bowtie_indices(dataset, chrom)
        # Creates BWA index.
        make_bwa_index(dataset, chrom)


def make_2bit_and_sizes(dataset, chrs):
    """Creates 2bit and sizes files for the artificial genomes of the dataset.

    Args:
        dataset: A Dataset instance of the new dataset.
        chrs: A list of chromosome numbers included in the dataset.

    Writes:
        All the 2bit and sizes files.
    """
    for chr_num in chrs:
        chrom = dataset.get_chr(chr_num)
        fa_file = chrom.get_path()
        twobit_file = chrom.get_path("2bit")
        sizes_file = chrom.get_path("sizes")

        _ = Popen(["faToTwoBit", fa_file, twobit_file]).wait()
        info = Popen(["twoBitInfo", twobit_file, "stdout"], stdout=PIPE)
        info.wait()
        sizes = Popen(["sort", "-k2rn"], stdin=info.stdout, stdout=PIPE)
        sizes.wait()
        with open(sizes_file, 'w') as fd:
            fd.write(sizes.stdout.read())


def main():
    """Prepares the files and directories for a new dataset.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset = DATASETS[args.dataset]
    dataset.set_work_dir(args.path)

    make_dirs(dataset)
    make_refgenes(dataset, args.chrs)
    make_indices(dataset, args.chrs)
    make_2bit_and_sizes(dataset, args.chrs)


if __name__ == '__main__':
    main()
