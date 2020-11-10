"""Runs phytoCRISP-Ex."""

import os
import shutil

import common
import tools_common


TOOL = "phytocrispex"
SRC_OUT_NAME = "%s-NGGpam.csv"
DST_OUT_NAME = "%s_out.csv"


def get_args():
    parser = common.default_parser("Runs phytoCRISP-Ex.")
    return parser.parse_args()


def main():
    """Runs phytoCRISP-Ex.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset, chrom, tool = common.process_args(args, tool_name=TOOL)
    if not dataset or not tool:
        return

    params = ["-g", args.dataset, "-c", chrom.name, "NGG", "G"]
    tool.run_bash(params)

    shutil.move(
        os.path.join(dataset.get_data_path(), SRC_OUT_NAME % chrom.name),
        os.path.join(tool.get_out_dir(dataset.name), DST_OUT_NAME % chrom.name)
    )
    tools_common.normalise(chrom, tool)


if __name__ == "__main__":
    main()