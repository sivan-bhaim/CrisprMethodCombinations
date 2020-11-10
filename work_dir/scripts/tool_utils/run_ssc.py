"""Runs SSC."""

import os

import common
import tools_common


TOOL = "ssc"
MATRIX_PATH = "matrix/human_mouse_CRISPR_KO_30bp.matrix"
SPACER_PATH = "bin/Fasta2Spacer"
RAW_OUT_NAME = "%s_out.txt"


def get_args():
    parser = common.default_parser("Runs SSC.")
    return parser.parse_args()


def main():
    """Runs SSC.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset, chrom, tool = common.process_args(args, tool_name=TOOL)
    if not dataset or not tool:
        return

    tool_path = tool.get_dir_path()
    out_path = tool.get_out_dir(dataset.name)

    spacer = os.path.join(tool_path, SPACER_PATH)
    spacer_out = os.path.join(out_path, "%s_spacer.txt" % chrom.name)
    tool.run_bash(["-i", chrom.get_path(), "-o", spacer_out], spacer)

    ssc = tool.get_bin_path()
    matrix = os.path.join(tool_path, MATRIX_PATH)
    ssc_out = os.path.join(out_path, RAW_OUT_NAME % chrom.name)
    stdout, stderr = tool.run_bash(
        ["-i", spacer_out, "-o", ssc_out, "-m", matrix, "-l", "30"], ssc)
    if stderr: print stderr
    if stdout: print stdout

    tools_common.create_csv(chrom, tool, RAW_OUT_NAME)
    tools_common.normalise(chrom, tool)


if __name__ == "__main__":
    main()
