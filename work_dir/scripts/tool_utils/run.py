"""Runs one of the prediction tools."""

import os
import subprocess

import common


SCRIPT_NAME = "run_%s.py"


def get_args(raw_args):
    parser = common.default_parser("Runs one of the prediction tools.")
    parser.add_argument(
        "-t", "--tool", type=str,
        help="Prefix of one of the tools' names.")
    return parser.parse_args(raw_args)


def main(raw_args=None):
    """Runs one of the prediction tools.

    The available tools are:
    - CHOPCHOP
    - FlashFry
    - mm10db
    - phytoCRISP-Ex
    - sgRNA Scorer 2.0
    - SSC

    For command line help, run with the '-h' flag.
    """
    # Identifies the tool from the prefix.
    args = get_args(raw_args)
    tool = common.get_tool(args.tool)
    if not tool:
        return
    print("--Running %s--" % tool)

    # Runs the script for the specific tool.
    script_dir = os.path.dirname(os.path.realpath(__file__))
    script_path = os.path.join(script_dir, SCRIPT_NAME % tool)
    params = [
        "-p", args.path,
        "-d", args.dataset,
        "-c", str(args.chr),
    ]
    process = subprocess.Popen(
        ["python", script_path] + params,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = process.communicate()
    if outs: print(outs)
    if errs: print(errs)


if __name__ == '__main__':
    main()
