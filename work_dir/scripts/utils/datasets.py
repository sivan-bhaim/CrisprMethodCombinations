"""
Represents datasets and chromosomes.

A dataset is usually composed of several chromosomes.
In the artificial genome sense, a chromosome is a specific part of the dataset.
For example, the Xu dataset is in fact a collection of multiple datasets, such
as ribosomal and nonribosomal. Thus, the Xu dataset here will have a chromosome
for the ribosomal set and a chromosome for the nonribosomal set.

Some of the logic in this code reflects the parsing methods of the raw datasets.
For the Xu and Doench datasets, the parsing methods are based on code provided
in the supplementary information of Bradford and Perrin (2019).
"""

import os
from os.path import join as join_path
from collections import OrderedDict

import consts


class Dataset(object):
    """Represents a dataset.

    Attributes:
        name: The name of the dataset (e.g. xu).
        sequence_idx: The index of the column with the target sequence in the
            dataset file.
        efficient_idx: The index of the column with the efficiency label in the
            dataset file.
        start_idx: The index of the column with the start position of the target
            within the genome (-1 indicates no such column exists).
        end_idx: The index of the column with the end position of the target
            within the genome (-1 indicates no such column exists).
        raw_activity_idx: The index of the column with the activity metric of
            the dataset.
        genome_label_idx: The index of the column with the GenomeCRISPR effect
            labels.
        target_pos: The position of the beginning of the target sequence within
            the full sequence provided in the dataset.
        chr_data: A dictionary mapping chromosome number to its name.
        work_dir: The main work directory.
    """

    def __init__(self, name, sequence_idx, efficient_idx, target_pos, chr_data,
                 start_idx=-1, end_idx=-1, raw_activity_idx=-1,
                 genome_label_idx=-1):
        self.name = name
        self.sequence_idx = sequence_idx
        self.efficient_idx = efficient_idx
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.raw_activity_idx = raw_activity_idx
        self.genome_label_idx = genome_label_idx
        self.target_pos = target_pos
        self.chr_data = chr_data
        self.work_dir = consts.WORK_DIR

    def set_work_dir(self, path):
        """Changes the working directory."""
        self.work_dir = path

    def get_data_path(self):
        """Returns the path to the data directory."""
        return join_path(self.work_dir, consts.DATA_DIR, self.name)

    def get_out_path(self):
        """Returns the path to the output directory."""
        return join_path(self.get_data_path(), consts.OUTPUT_DIR)

    def get_indexes_path(self):
        """Returns the path to the indexes directory (e.g. bowtie)."""
        return join_path(self.get_data_path(), consts.INDEXES_DIR)

    def get_raw_data_path(self, chr_num):
        """Returns the path to the chromosome's CSV file."""
        filename = consts.RAW_DATA_CSV % (chr_num, self.chr_data[chr_num])
        return os.path.join(self.get_data_path(), consts.RAW_DATA_DIR, filename)

    def get_target(self, values):
        """Extracts the target sequence from a CSV entry.

        Args:
            values: The list of values in the CSV entries.

        Returns:
            The target sequence.
        """
        start_pos = self.target_pos
        end_pos = start_pos + consts.TARGET_LEN
        sequence = values[self.sequence_idx]
        return sequence[start_pos:end_pos].upper()

    def get_targets(self, chr_num, length=consts.TARGET_LEN, condition=None):
        """Returns the list of targets in the dataset which pass the condition.

        Args:
            chr_num: The chromosome number.
            length: The length of the targets.
            condition: A condition which must be met by the target for it to be
                included. This should be a function with a single argument (the
                list of values in the entry for the target) which returns a
                bool. If no condition is provided, all the targets in the
                dataset are included.

        Returns:
            The list of target sequences.
        """
        targets = []

        with open(self.get_raw_data_path(chr_num), 'r') as fd:
            for line in fd.readlines()[1:]:
                fields = from_csv_line(line)
                if condition and not condition(fields): continue
                targets.append(self.get_target(fields))
        return targets


    def get_value(self, chr_num, target, column):
        """Returns a specific value for a target.

        Args:
            chr_num: The chromosome number.
            target: The target sequence.
            column: The index of the coloumn to retrieve.

        Returns:
            The desired value.
        """
        with open(self.get_raw_data_path(chr_num), 'r') as fd:
            for line in fd.readlines()[1:]:
                fields = from_csv_line(line)
                if self.get_target(fields) == target.upper():
                    return fields[column]

    def get_efficient_targets(self, chr_num, no_pam=False):
        """Returns a list of efficient targets.

        Args:
            chr_num: The chromosome number.
            no_pam: A boolean indicating whether to remove the PAM sequences.

        Returns:
            A list of efficient targets.
        """
        condition = lambda fields: int(fields[self.efficient_idx])
        targets = self.get_targets(chr_num, condition=condition)
        if no_pam:
            return remove_pam(targets)
        return targets

    def get_inefficient_targets(self, chr_num,  no_pam=False):
        """Returns a list of inefficient targets.

        Args:
            chr_num: The chromosome number.
            no_pam: A boolean indicating whether to remove the PAM sequences.

        Returns:
            A list of inefficient targets.
        """
        condition = lambda fields: not int(fields[self.efficient_idx])
        targets = self.get_targets(chr_num, condition=condition)
        if no_pam:
            return remove_pam(targets)
        return targets

    def get_chr(self, chr_num):
        """Returns a Chromosome instance."""
        return Chromosome(chr_num, self.chr_data[chr_num], self, self.work_dir)


class Chromosome(object):
    """Represents a specific chromosome in the dataset.

    Attributes:
        num: The number of the chromosome.
        name: A string representing the chromosome's number (chr<num>).
        description: A descriptive name of the chromosome (e.g. ribosomal).
        dataset: A Dataset instance of the dataset the chromosome belongs to.
        work_dir: The main work directory.
    """
    def __init__(self, num, description, dataset, work_dir=consts.WORK_DIR):
        self.num = num
        self.name = "chr%d" % num
        self.description = description
        self.dataset = dataset
        self.work_dir = work_dir

    def get_name(self):
        """Returns the name of the chromosome."""
        return self.name

    def get_path(self, extension="fa"):
        """Returns the path of the file of the chromosome's sequence.

        Args:
            extension: The extension of the file.

        Returns:
            The path to the sequence file.
        """
        file_name = "%s.%s" % (self.name, extension.strip('.'))
        return join_path(self.dataset.get_data_path(), file_name)

    def get_refgene(self, refgene=consts.REFGENE_FNAME):
        """Returns the path of the refGene file for the chromosome.

        Args:
            refgene: A template string of the refGene file name. Should have a
                single %d into which the chromosome's number will be placed.

        Returns:
            The path to the refGene file.
        """
        return join_path(self.dataset.get_data_path(), refgene % self.num)

    def add_values(self, new_headers, target_to_vals):
        """Adds values to the chromosome's CSV file.

        Args:
            new_headers: A list of the headers for the new columns.
            target_to_vals: A dictionary mapping each target sequence in the
                dataset to a list of the new values. The ordering of these
                values must match the ordering of the headers.
        """
        data = OrderedDict()
        with open(self.dataset.get_raw_data_path(self.num), 'r') as fd_in:
            lines = fd_in.readlines()
            headers = from_csv_line(lines[0]) + new_headers
            for line in lines[1:]:
                values = from_csv_line(line)
                target = self.dataset.get_target(values)
                data[target] = values + target_to_vals[target]

        with open(self.dataset.get_raw_data_path(self.num), 'w') as fd_out:
            fd_out.write(to_csv_line(*headers))
            for vals in data.values():
                fd_out.write(to_csv_line(*vals))


def remove_pam(targets):
    """Removes the PAM sequence from all the targets."""
    return [target[:-consts.PAM_LEN] for target in targets]


def to_csv_line(*args):
    """Joins all the values into a CSV line."""
    return ','.join([str(value) for value in args]) + '\n'


def from_csv_line(line):
    """Splits a CSV line into a list of values."""
    return line.strip().split(',')


XU_DATASET = Dataset(
    name="xu",
    sequence_idx=4,
    efficient_idx=9,
    target_pos=10,
    chr_data={
        1: "ribo",
        2: "nonribo",
    },
    raw_activity_idx=7,
    genome_label_idx=8,
)

DOENCH_DATASET = Dataset(
    name="doench",
    sequence_idx=1,
    efficient_idx=10,
    target_pos=4,
    chr_data={
        1: "train",
    },
    raw_activity_idx=7,
    genome_label_idx=9,
)

GENOMECRISPR_DATASET = Dataset(
    name="genome",
    sequence_idx=7,
    efficient_idx=16,
    target_pos=0,
    chr_data={
        1: "15k",
    },
    start_idx=0,
    end_idx=1,
    genome_label_idx=13,
)

CHARI_DATASET = Dataset(
    name="chari",
    sequence_idx=1,
    efficient_idx=6,
    target_pos=0,
    chr_data={
        1: "293t",
    },
    raw_activity_idx=2,
    genome_label_idx=5,
)

DATASETS = {
    "xu": XU_DATASET,
    "doench": DOENCH_DATASET,
    "genome": GENOMECRISPR_DATASET,
    "chari": CHARI_DATASET,
}
