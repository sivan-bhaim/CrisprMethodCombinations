"""Runs mm10db."""

import os

import consts
import common
import tools_common


TOOL = "mm10db"
# Files created by mm10db
RAW_OUT_NAME = "%s_accepted_targets.txt"
REJECT_FILE_NAME = "chr%d_rejected_targets.txt"
# Scripts required for running mm10db, in order.
SCRIPTS = [
    "createListExons.py",
    "prepareExonSequences.py",
    "prepareListOfftargetSites.py",
    "target_identitification_viaC.py",
]
# Rejections reasons.
AT = "AT%"
MULTI_MATCH = "Multiple matches in exons"
STRUCT_ENERGY = "Secondary structure or energy"
REVERSE_PRIMER = "Too close to reverse primer"
TTTT = "TTTT"
OFF_TARGET = "Off-target score"
# Encoding of rejection reasons as error codes. These will appear in the output file.
REJECTIONS = {
    AT: -1,
    MULTI_MATCH: -2,
    STRUCT_ENERGY: -3,
    REVERSE_PRIMER: -4,
    TTTT: -5,
    OFF_TARGET: -6,
}


def get_args():
    parser = common.default_parser("Runs mm10db.")
    return parser.parse_args()


def add_rejected(chrom, out_path):
    """Adds the rejected targets to the output, with the rejection error codes.

    Args:
        chrom: A Chromosome instance.
        out_path: Path to the directory of the tool's output for the processed
            dataset.

    Writes:
        The complete output CSV file.
    """
    # Collects rejected spacers with their rejection reasons.
    rejections = []
    rejection_file_path = os.path.join(out_path, REJECT_FILE_NAME % chrom.num)
    with open(rejection_file_path, 'r') as fd:
        for line in fd:
            seq, reason = line.strip().split('\t')
            if reason not in REJECTIONS:
                # Stops if there is an unknown rejection reason.
                print "Unknown rejection reason: %s" % reason
                return
            rejections.append((seq[:-consts.PAM_LEN], REJECTIONS[reason]))

    # Adds the rejected targets to the output file.
    out_file_path = os.path.join(out_path, consts.RAW_OUT_CSV % chrom.num)
    with open(out_file_path, 'a') as fd:
        for seq, reason in rejections:
            line = common.to_csv_line(seq, "", "", "", reason, chrom.name,
                                      -1, -1, "x", "", "", "", "")
            fd.write(line)


def main():
    """Runs mm10db.

    For command line help, run with the '-h' flag.
    """
    args = get_args()
    dataset, chrom, tool = common.process_args(args, tool_name=TOOL)
    if not dataset or not tool:
        return
    tool.chdir()

    script_params = ["-g", args.dataset, "-c", chrom.name]
    out_path = tool.get_out_dir(dataset.name)

    # Runs the mm10db pipeline.
    for script in SCRIPTS:
        tool.run_python(script_params, script)

    tools_common.create_csv(chrom, tool, RAW_OUT_NAME)
    add_rejected(chrom, out_path)
    tools_common.normalise(chrom, tool)



if __name__ == "__main__":
    main()
