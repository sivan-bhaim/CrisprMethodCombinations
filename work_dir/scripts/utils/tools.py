"""Represents prediction tools."""

from os import chdir
from os.path import join as join_path
from subprocess import Popen, PIPE
from collections import OrderedDict

import consts
import common


class Tool(object):
    """Represents prediction tools.

    Attributes:
        name: The name of the tool.
        seq_idx: The index of the sequence column in the tool's output.
        strand_idx: The index of the strand column in the tool's output.
        start_idx: The index of the sequence start position column in the tool's
            output.
        end_idx: The index of the sequence end position column in the tool's
            output.
        scores: A mapping from score column indices to their names.
        has_headers: A bool indicating whether the tool's output file has
            a headers line.
        work_dir: Path to the main work directory.
        bin_path: Path to the main binary or script of the tool, starting from
            the tool's main directory.
    """
    def __init__(self, name, seq_idx, strand_idx, start_idx, end_idx=-1, scores=None,
                 has_headers=True, work_dir=consts.WORK_DIR, bin_path=None):
        self.name = name
        self.seq_idx = seq_idx
        self.strand_idx = strand_idx
        self.start_idx = start_idx
        self.end_idx = start_idx if end_idx == -1 else end_idx
        self.scores = scores if scores else {}
        self.has_headers = has_headers
        self.work_dir = work_dir
        self.bin_path = bin_path

    def to_csv_line(self, *args):
        """Returns a CSV line with the args, prefixed by the tool's name."""
        return common.to_csv_line(self.name, *args)

    def normalise(self, entries):
        """Normalises the tool's output entries."""
        normalised = []
        for entry in entries:
            values = self.preprocess(entry.split(','))

            seq = self.extract_seq(values[self.seq_idx])
            if not seq: continue

            strand = self.extract_strand(values[self.strand_idx])
            start, end = self.extract_pos(
                values[self.start_idx], values[self.end_idx], strand)
            scores = self.extract_scores(values)

            normalised.append(
                self.to_csv_line(seq, start, end, strand, *scores))

        return normalised

    def preprocess(self, values):
        """Performs pre-processing."""
        return values

    def extract_seq(self, value):
        """Extracts the target sequence."""
        return value

    def extract_pos(self, value_start, value_end, strand):
        """Extracts the start and end positions of the target."""
        start = int(value_start)
        if value_start == value_end:
            end = start + consts.TARGET_LEN - 1
        else:
            end = int(value_end)
        return start, end

    def extract_strand(self, value):
        """Extracts the strand."""
        return value

    def extract_scores(self, values):
        """Extracts the list of scores."""
        return [values[idx] for idx in self.scores.keys()]

    def get_score_idx(self, name):
        """Returns the column index of the given score name."""
        for i, score in enumerate(self.scores.values()):
            if name == score:
                return i
        return -1

    def set_work_dir(self, path):
        """Changes the working directory."""
        self.work_dir = path

    def get_out_dir(self, dataset):
        """Returns the output directory of the tool.

        Args:
            dataset: The dataset for which the output is required.

        Returns:
             The path of the directory.
        """
        return join_path(common.out_path(dataset, self.work_dir), self.name)

    def _get_out_file_path(self, pattern, chr_num, dataset):
        """Returns the path of an output file.

        Args:
            pattern: A pattern for the name of the output file. Should have a
                single %d for the number of the chromosome.
            chr_num: The number of the chromosome for which the output is
                required.
            dataset: The name of the dataset for which the output is required.

        Returns:
            The path to the output file.
        """
        return join_path(self.get_out_dir(dataset), pattern % chr_num)

    def get_normalised_out_path(self, chr_num, dataset):
        """Returns the path of the normalised output file.

        Args:
            chr_num: The number of the chromosome for which the output is
                required.
            dataset: The name of the dataset for which the output is required.

        Returns:
            The path to the output file.
        """
        return self._get_out_file_path(consts.NORMALISED_OUT, chr_num, dataset)

    def get_csv_out_path(self, chr_num, dataset):
        """Returns the path of the raw CSV output file.

        Args:
            chr_num: The number of the chromosome for which the output is
                required.
            dataset: The name of the dataset for which the output is required.

        Returns:
            The path to the output file.
        """
        return self._get_out_file_path(consts.RAW_OUT_CSV, chr_num, dataset)

    def get_dir_path(self):
        """Returns the path of the tool's main directory."""
        return join_path(common.get_tools_path(self.work_dir), self.name)

    def get_bin_path(self):
        """Returns the path of the tool's main binary or script."""
        return join_path(self.get_dir_path(), self.bin_path)

    def chdir(self):
        """Changes directory to the main directory of the tool."""
        chdir(self.get_dir_path())

    def _run(self, full_params):
        """Runs the command given by full_params."""
        print("Executing: %s" % (' '.join(full_params),))
        process = Popen(full_params, stdout=PIPE)
        return process.communicate()

    def run_bash(self, params, bin_path=None):
        """Runs the tool with the given params as a bash command.

        Args:
            params: Command line parameters for the tool.
            bin_path: Path to the main binary or script of the tool. If not
                provided, the default path is used.

        Returns:
            A tuple (stdout_data, stderr_data).
        """
        if not bin_path: bin_path = self.bin_path
        full_path = join_path(self.get_dir_path(), bin_path)
        return self._run([full_path] + params)

    def run_python(self, params, script_path=None):
        """Runs the tool with the given params as a Python script.

        Args:
            params: Command line parameters for the tool.
            script_path: Path to the main script of the tool. If not provided,
                the default path is used.

        Returns:
            A tuple (stdout_data, stderr_data).
        """
        if not script_path: script_path = self.bin_path
        full_path = join_path(self.get_dir_path(), script_path)
        return self._run(["python", full_path] + params)


# Following are the configurations of the supported prediction tools.

class Chopchop(Tool):
    def __init__(self):
        scores = OrderedDict([
            (4, "gc_content"),
            (5, "self_complementarity"),
            (10, "xu_2015"),
            (11, "doench_2014"),
            (13, "moreno_mateos_2015"),
            (15, "g20"),
        ])
        super(Chopchop, self).__init__(
            name="chopchop",
            seq_idx=1,
            strand_idx=3,
            start_idx=2,
            scores=scores,
            bin_path="chopchop.py"
        )

    def extract_pos(self, value_start, value_end, strand):
        start = int(value_start.split(':')[1])
        end = start + consts.TARGET_LEN - 1
        return start, end

    def extract_scores(self, values):
        scores = [values[idx] for idx in self.scores.keys()]
        # The first two scores require scaling.
        normalised = [float(scores[0])/100, float(scores[1])/5]
        return normalised + scores[2:]


class FlashFry(Tool):
    def __init__(self):
        scores = OrderedDict([
            (7, "doench_2014"),
            (8, "moreno_mateos_2015"),
        ])
        super(FlashFry, self).__init__(
            name="flashfry",
            seq_idx=3,
            strand_idx=6,
            start_idx=1,
            end_idx=2,
            scores=scores,
            bin_path="flashfry.jar"
        )

    def extract_strand(self, value):
        return {'FWD': '+', 'RVS': '-'}[value]

    def extract_pos(self, value_start, value_end, strand):
        start = int(value_start) + 1
        end = int(value_end)
        return start, end


class mm10db(Tool):
    def __init__(self):
        scores = OrderedDict([
            (4, "decision"),
        ])
        super(mm10db, self).__init__(
            name="mm10db",
            seq_idx=0,
            strand_idx=8,
            start_idx=6,
            end_idx=7,
            has_headers=False,
            scores=scores
        )

    def extract_seq(self, value):
        return "%s..." % value

    def extract_pos(self, value_start, value_end, strand):
        start = int(value_start)
        end = int(value_end)
        if strand == '+': end += 3
        else: start -= 3
        return start, end

    def extract_scores(self, values):
        score = float(values[self.scores.keys()[0]])
        if score < 0:
            return [score]
        return [1]

class phytoCRISPEX(Tool):
    def __init__(self):
        scores = OrderedDict([
            (-1, "decision"),
        ])
        super(phytoCRISPEX, self).__init__(
            name="phytocrispex",
            seq_idx=4,
            strand_idx=1,
            start_idx=2,
            end_idx=3,
            scores=scores,
            bin_path="bin/phytoCRISPex"
        )

    def preprocess(self, values):
        return values[0].split('_')

    def extract_scores(self, values):
        # A target that appears in the output is accepted and gets the score 1.
        return [1]


class sgRNAScorer(Tool):
    def __init__(self):
        scores = OrderedDict([
            (2, "score"),
        ])
        super(sgRNAScorer, self).__init__(
            name="sgrnascorer",
            seq_idx=1,
            strand_idx=0,
            start_idx=0,
            scores=scores,
            bin_path="identifyAndScore.py"
        )

    def extract_strand(self, value):
        strand = value.split('_')[1]
        return {'Plus': '+', 'Minus': '-'}[strand]

    def extract_pos(self, value_start, value_end, strand):
        start = int(value_start.split('_')[2]) + 1
        end = start + consts.TARGET_LEN - 1
        return start, end

    def extract_scores(self, values):
        score = float(values[self.scores.keys()[0]])
        # We scale the scores by a factor 5 so they roughly fit in the range [-1,1].
        return [score / 5.0]


class SSC(Tool):
    def __init__(self):
        scores = OrderedDict([
            (5, "decision"),
        ])
        super(SSC, self).__init__(
            name="ssc",
            seq_idx=0,
            strand_idx=3,
            start_idx=1,
            has_headers=False,
            scores=scores,
            bin_path="bin/SSC"
        )

    def extract_seq(self, value):
        seq = value[:consts.TARGET_LEN]
        if 'N' in seq: return ""
        return seq

    def extract_pos(self, value_start, value_end, strand):
        start = int(value_start) + 1
        end = start + consts.TARGET_LEN - 1
        if strand == '-':
            start += 7
            end += 7
        return start, end


TOOLS = OrderedDict([
    ("chopchop", Chopchop()),
    ("flashfry", FlashFry()),
    ("mm10db", mm10db()),
    ("phytocrispex", phytoCRISPEX()),
    ("sgrnascorer", sgRNAScorer()),
    ("ssc", SSC()),
])
