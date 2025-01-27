"""
This script extracts the sequence with the highest count.

Using a given input file and writes it to an output file.
"""

import re
import sys

import pysam

from .util import wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    with pysam.FastxFile(args.input_fasta) as f:
        max_count = 0
        max_seq = None

        for entry in f:
            (count,) = re.search(r";seqs=(\d+)", entry.name).groups()
            count = int(count)
            if count > max_count:
                max_count = count
                max_seq = entry.sequence

    if max_count >= args.min_seqs:
        cluster_name = args.input_fasta.split(".")[0]
        sys.stdout.write(f">{cluster_name}\n{max_seq}\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("ExtHigh")
    parser.add_argument(
        "input_fasta",
        help="VSEARCH output consensus FASTA file",
    )
    parser.add_argument(
        "min_seqs",
        type=int,
        help="Minimum number of sequences of a cluster to keep",
    )
    return parser
