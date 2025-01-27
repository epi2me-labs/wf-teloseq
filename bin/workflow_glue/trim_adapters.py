"""Trim adapter sequences from the 5' end of sequences in a FASTQ file."""

import re
import sys

import pysam

from .util import wf_parser  # noqa: ABS101

# TODO: This scripts should be folded into its consumer/producer


def main(args):
    """Run the entry point.

    Then writes the modified sequences to the output file.

    Parameters:
    fastq_path (str): Path to the input FASTQ file.
    output_path (str): Path to the output FASTQ file.
    ADAPTER_SEQUENCES (list): List of sequences to search for trimming.
    """
    with pysam.FastxFile(args.input) as f:
        for record in f:
            seq = record.sequence
            qual = record.quality
            earliest_end_pos = float("inf")

            # Consider only the 5' end and first 250 bp
            sequence_fragment = str(seq[:250])

            for search_seq in args.adapter_seqs.split(","):
                search_seq_clean = re.sub(r"\W+", "", search_seq)
                match = re.search(search_seq_clean, sequence_fragment)
                if match is not None and match.end() < earliest_end_pos:
                    earliest_end_pos = match.end()

            if earliest_end_pos != float("inf"):
                # Modify the sequence by removing from 1 to the (end position - 15)
                # Also, trim the quality scores accordingly
                trim_pos = max(0, earliest_end_pos - 14)
                record.sequence = seq[trim_pos:]
                record.quality = qual[trim_pos:]

                sys.stdout.write(str(record) + "\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("TrimAdap")
    parser.add_argument("input", help="FASTQ file with reads to trim")
    parser.add_argument("adapter_seqs", help="comma-separated list of adapters to trim")
    return parser
