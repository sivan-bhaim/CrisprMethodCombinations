"""Runs sgRNA Scorer 2.0."""

import os
import subprocess

import common
import tools_common


RAW_OUT_NAME = "%s_out.txt"
TOOL = "sgrnascorer"


def get_args():
    parser = common.default_parser("Runs sgRNA Scorer 2.0.")
    return parser.parse_args()


def main():
    """Runs sgRNA Scorer 2.0.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset, chrom, tool = common.process_args(args, tool_name=TOOL)
    if not dataset or not tool:
        return
    tool.chdir()

    script_path = os.path.join(tool.get_dir_path(), tool.bin_path)
    params = [
        "-d", args.dataset,
        "-c", str(args.chr),
        "-o", RAW_OUT_NAME % chrom.name,
        "-p", "3", "-s", "20", "-l", "NGG",
    ]
    process = subprocess.Popen(
        ["python", script_path] + params,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = process.communicate()
    print(outs)
    print(errs)

    tools_common.create_csv(chrom, tool, RAW_OUT_NAME)
    tools_common.normalise(chrom, tool)


if __name__ == "__main__":
    main()
