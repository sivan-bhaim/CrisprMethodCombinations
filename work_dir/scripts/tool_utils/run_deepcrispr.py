"""
Runs DeepCRISPR on the desired dataset.
This script should be run inside the DeepCRISPR Docker container, available from
the Docker repository:
docker pull michaelchuai/deepcrispr:latest

Specifically, to import the deepcrispr module, the script should run from
/root/DeepCRISPR (inside the container).
We recommend running the container with a volume, to allow files to be moved in
and out easily:
sudo docker run --name deep -v <local path>:<container path> --rm -i -t michaelchuai/deepcrispr bash
"""

import re
import os
import numpy as np
import pickle
import argparse
import tensorflow as tf
from deepcrispr import DCModelOntar


# Identifies a target sequence inside the artificial genome.
TARGET_PATTERN = '([ACGT]{%d})[N]{50}'
# Path to the DeepCRISPR models directory.
MODELS_DIR = "/root/DeepCRISPR/trained_models"
# Maps DeepCRISPR model options to the model's name.
# The first key signifies whether epigenetic features are used. The second key
# signifies whether the model uses regression.
MODELS = {
    True: {
        True: "ontar_pt_cnn_reg",
        False: "ontar_ptaug_cnn",
    },
    False: {
        True: "ontar_cnn_reg_seq",
    }
}


def get_model(is_epigenetics, is_regression):
    """Returns the path to the required model.

    Args:
        is_epigenetics: Should the model include epigenetic features.
        is_regression: Should the model be a regressuion model.

    Returns:
        Full path to the required model.
    """
    model_name = MODELS[is_epigenetics][is_regression]
    return  os.path.join(MODELS_DIR, model_name)

def encode(seq, channels):
    """Encodes a target sequence according to the DeepCRISPR input format.

    Args:
        seq: Target sequence.
        channels: The number of required channels.

    Returns:
        The encoded target sequence.
    """
    map = {'A':0, 'C':1, 'G': 2, 'T':3}
    out = [[[]] for _ in range(channels)]
    for c in seq:
        for i in range(channels):
            if i == map[c]: out[i][0].append(1)
            else: out[i][0].append(0)
    return out


def get_args():
    parser = argparse.ArgumentParser(
        description="Run DeepCRISPR.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-s", "--sequence", type=str,
        help=("Path to a single-line sequence file from which the targets will "
              "be extracted. Mutually exclusive with --targets."))
    parser.add_argument(
        "-t", "--targets", type=str,
        help=("Path to a file all of the targets, each on in a separate line. "
              "Mutually exclusive with --sequence."))
    parser.add_argument(
        "-o", "--out_dir", type=str,
        help=("Path to the main DeepCRISPR output directory. The complete "
              "output path will then depend on the variant used and the name "
              "of the dataset."))
    parser.add_argument(
        "-d", "--dataset", type=str,
        help=("The name of the dataset. Used to choose an appropriate output "
              "file name."))
    parser.add_argument(
        "-r", "--regression", action="store_true",
        help="Use a regression model (otherwise classification is used.")
    parser.add_argument(
        "-e", "--epigenetics", action="store_true",
        help="Use an epigenetics model (otherwise sequence-only is used.)")
    parser.add_argument(
        "-l", "--target_len", type=int, default=23,
        help=("Length of target sequences. Not necessary when --targets is "
              "used."))
    return parser.parse_args()


def get_targets(path, target_len, is_sequence):
    """Returns the list of targets.

    Args:
        path: Path to the file with the targets.
        target_len: The length of the targets.
        is_sequence: Is the file a single-sequence of the artificial genome.

    Returns:
        The list of targets.
    """
    with open(path) as fd:
        if not is_sequence:
            return [target.strip() for target in fd]
        sequence = fd.read()
    # Identifies all the target sequences in the chromosome.
    return re.findall(TARGET_PATTERN % (target_len,), sequence)


def main():
    """Runs DeepCRISPR.

    For additional instructions, see comment at the top of this file.
    For command line help, run with the '-h' flag.

    Writes:
        A Pickle file, which holds a mapping from target sequences to their
        DeepCRISPR scores.
        The output directory will be:
        <out_dir>/<regression or classification>/<epigenetics or seq_only>
        And the file name will be:
        <dataset name>_results.pkl
    """
    args = get_args()
    if args.sequence:
        path = args.sequence
        is_sequence = True
    else:
        path = args.targets
        is_sequence = False

    # Creates the input for the model.
    targets = get_targets(path, args.target_len, is_sequence)
    channels = 8 if args.epigenetics else 4
    x_on_target = np.array(
        [encode(t, channels) for t in targets]
    )

    # Runs the model.
    sess = tf.InteractiveSession()
    model_path = get_model(args.epigenetics, args.regression)
    dcmodel = DCModelOntar(
        sess, model_path, args.regression, not args.epigenetics)
    predicted_on_target = dcmodel.ontar_predict(x_on_target)

    # Maps targets to their scores.
    results = {}
    for i, target in enumerate(targets):
        score = predicted_on_target[i]
        results[target] = score

    # Saves the results.
    out_path = os.path.join(
        args.out_dir, "data", "deepcrispr",
        "regression" if args.regression else "classification",
        "epigenetics" if args.epigenetics else "seq_only",
        "%s_results.pkl" % args.dataset,
    )

    with open(out_path, 'wb') as fd:
        pickle.dump(results, fd)


if __name__ == '__main__':
    main()
