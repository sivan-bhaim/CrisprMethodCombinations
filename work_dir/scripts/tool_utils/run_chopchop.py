"""Runs CHOPCHOP."""

import os

import common
import tools_common


RAW_OUT_NAME = "%s_out.txt"
TOOL = "chopchop"
CONFIG_FILE = "config.json"
PLACE_HOLDER = "DATASET_DIR"
CONFIGURATION = """
{{
  "PATH": {{
    "PRIMER3": "./primer3_core",
    "BOWTIE": "bowtie/bowtie",
    "TWOBITTOFA": "./twoBitToFa",
    "TWOBIT_INDEX_DIR": "{PLACE_HOLDER}",
    "BOWTIE_INDEX_DIR": "{PLACE_HOLDER}/indexes/bowtie1",
    "ISOFORMS_INDEX_DIR": "/your/full/path/to/ebwt_transcriptome_folder_and_2bit_of_genome",
    "ISOFORMS_MT_DIR": "/your/full/path/to/vienna_MT_folder",
    "GENE_TABLE_INDEX_DIR": "{PLACE_HOLDER}/chopchop"
  }},
  "THREADS": 1
}}
""".format(**locals())


def config(dataset_dir, config_path):
    """Writes the configuration for CHOPCHOP.

    Args:
        dataset_dir: The directory with the chromosome files of the dataset.
        config_path: The path to the configuration file for CHOPCHOP.
    """
    configuration = CONFIGURATION.replace(PLACE_HOLDER, dataset_dir)
    with open(config_path, 'w') as fd_write:
        fd_write.write(configuration)


def unconfig(config_path):
    """Reverts the configuration file."""
    with open(config_path, 'w') as fd_write:
        fd_write.write(CONFIGURATION)


def get_args():
    parser = common.default_parser("Runs CHOPCHOP.")
    return parser.parse_args()


def main():
    """Runs CHOPCHOP.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset, chrom, tool = common.process_args(args, tool_name=TOOL)
    if not dataset or not tool:
        return
    tool.chdir()

    chr_name = chrom.get_name()
    out_path = tool.get_out_dir(dataset.name)
    config_path = os.path.join(tool.get_dir_path(), CONFIG_FILE)

    # Configures and runs CHOPCHOP, and then reverts the configuration file.
    config(dataset.get_data_path(), config_path)
    params = [
        "-G", dataset.name,
        "-o", os.path.join(out_path, chr_name),
        "-F", "--targets", chrom.get_path(),
        "--chr", chr_name,
        "--scoringMethod", "ALL",
    ]
    output, _ = tool.run_python(params)
    unconfig(config_path)

    # Writes the output.
    with open(os.path.join(out_path, RAW_OUT_NAME % chr_name), 'w') as txt:
       txt.write(output)
    tools_common.create_csv(chrom, tool, RAW_OUT_NAME)
    tools_common.normalise(chrom, tool)


if __name__ == "__main__":
    main()