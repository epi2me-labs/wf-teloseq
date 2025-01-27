"""Extend telomere in FASTA file."""
import pathlib
import sys

import pysam

from .util import wf_parser  # noqa: ABS101

# TODO: This script should be folded into its producer/consumer


def main(args):
    """Run the entry point."""
    with pysam.FastxFile(args.input_file) as infile:
        for entry in infile:
            # Trim the sequence to the first occurrence of 'repeat'. Then extend it by
            # adding repeats to the start of the sequence.
            trim_index = entry.sequence.find(args.repeat_sequence)
            if trim_index != -1:
                # Trim to the first occurrence of the repeat
                trimmed_seq = entry.sequence[trim_index:]
                # Add repeat to the start of the sequence
                extended_seq = args.repeat_sequence * args.repeat_count + trimmed_seq

            sys.stdout.write(f">{entry.name}\n{extended_seq}\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("ExtTelo")
    parser.add_argument(
        "input_file",
        type=pathlib.Path,
        help="FASTA file with telomere->subtelomere sequences",
    )
    parser.add_argument(
        "repeat_sequence",
        help="sequence of repeat to use for extending the telomere",
    )
    parser.add_argument(
        "repeat_count",
        type=int,
        help="number of repeats to use to extend the telomere",
    )
    return parser
