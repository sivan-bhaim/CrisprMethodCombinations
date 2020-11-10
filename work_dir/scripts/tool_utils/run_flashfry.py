"""Runs FlashFry."""

import os
import subprocess

import common
import tools_common


TOOL = "flashfry"
INDEX_NAME = "%s_index.bin"
DISCOVER_NAME = "%s_discover.txt"
RAW_OUT_NAME = "%s_out.txt"


def get_args():
    parser = common.default_parser("Runs FlashFry.")
    return parser.parse_args()


def main():
    """Runs FlashFry.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset, chrom, tool = common.process_args(args, tool_name=TOOL)
    if not dataset or not tool:
        return

    chr_path = chrom.get_path()
    out_path = tool.get_out_dir(dataset.name)
    bin_path = tool.get_bin_path()

    chr_name = chrom.get_name()
    index_path = os.path.join(out_path, INDEX_NAME % chr_name)
    discover_path = os.path.join(out_path, DISCOVER_NAME % chr_name)
    targets_path = os.path.join(out_path, RAW_OUT_NAME % chr_name)

    cmd_index = [
        "java", "-Xmx4g",
        "-jar", bin_path,
        "index",
        "--tmpLocation", out_path,
        "--database", index_path,
        "--reference", chr_path,
        "--enzyme", "spcas9ngg",
    ]
    cmd_discover = [
        "java", "-Xmx4g",
        "-jar", bin_path,
        "discover",
        "--database", index_path,
        "--fasta", chr_path,
        "--output", discover_path,
    ]
    cmd_score = [
        "java", "-Xmx4g",
        "-jar", bin_path,
        "score",
        "--input", discover_path,
        "--output", targets_path,
        "--scoringMetrics", "doench2014ontarget,moreno2015",
        "--database", index_path,
    ]

    tool.chdir()
    _ = subprocess.Popen(cmd_index).wait()
    _ = subprocess.Popen(cmd_discover).wait()
    _ = subprocess.Popen(cmd_score).wait()

    tools_common.create_csv(chrom, tool, RAW_OUT_NAME)
    tools_common.normalise(chrom, tool)


if __name__ == "__main__":
    main()
