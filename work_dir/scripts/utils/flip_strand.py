"""Flips a DNA sequence (complement and reverse)."""

import argparse


COMPLEMENT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'N': 'N',
}


def get_args():
    parser = argparse.ArgumentParser(
        description="Complements and reverses sequence.")
    parser.add_argument(
        "sequence", type=str,
        help="The sequence to flip to the other strand")
    return parser.parse_args()


def main():
    """Flips a DNA sequence (complement and reverse).

    Usage:
        flip_strand.py <sequence>

    Prints:
        The flipped sequence.
    """
    sequence = get_args().sequence
    flipped = ''.join([COMPLEMENT[b] for b in sequence])[::-1]
    print flipped


if __name__ == '__main__':
    main()
