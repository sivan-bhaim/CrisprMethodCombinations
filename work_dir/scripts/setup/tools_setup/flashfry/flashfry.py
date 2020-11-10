"""
A wrapper script which operates FlashFry.
"""

import os
import argparse
import subprocess


WORK_DIR = os.environ['WORK_DIR']
TOOL_NAME = "flashfry.jar"


def get_args():
    parser = argparse.ArgumentParser(description="Run FlashFry.")
    parser.add_argument(
        "-p", "--path", type=str, default=WORK_DIR,
        help="Path to the main work directory")
    parser.add_argument(
        "-c", "--chr", type=str,
        help="The number of the chromosome")
    parser.add_argument(
        "-d", "--dataset", type=str,
        help="The name of the dataset")
    return parser.parse_args()


def main():
    args = get_args()
    
    data_dir = os.path.join(args.path, "data")
    data_path = os.path.join(data_dir, args.dataset)
    chr_path = os.path.join(data_path, "chr%s.fa" % args.chr)
    output_path = os.path.join(data_path, "flashfry")
    
    index_path = os.path.join(output_path, "index.bin")
    discover_path = os.path.join(output_path, "discover.txt")
    targets_path = os.path.join(output_path, "targets.txt")
    
    cmd_index = [
        "java", "-Xmx4g",
        "-jar", TOOL_NAME,
        "index",
        "--tmpLocation", output_path,
        "--database", index_path,
        "--reference", chr_path,
        "--enzyme", "spcas9ngg",
    ]
    cmd_discover = [
        "java", "-Xmx4g",
        "-jar", TOOL_NAME,
        "discover",
        "--database", index_path,
        "--fasta", chr_path,
        "--output", discover_path,
    ]
    cmd_score = [
        "java", "-Xmx4g",
        "-jar", TOOL_NAME,
        "score",
        "--input", discover_path,
        "--output", targets_path,
        "--scoringMetrics", "doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot",
        "--database", index_path,
    ]

    _ = subprocess.Popen(cmd_index).wait()
    _ = subprocess.Popen(cmd_discover).wait()
    _ = subprocess.Popen(cmd_score).wait()


if __name__ == "__main__":
    main()
