"""Processes a BAM file to count indels above a certain threshold."""
from collections import defaultdict
import pathlib
import sys

import pysam

from .util import wf_parser  # noqa: ABS101

# TODO: This script should be folded into its producer/consumer


def main(args):
    """Run the entry point."""
    with pysam.AlignmentFile(args.bam_file, "rb") as bam:
        # results dict: {contig_name: {read_id: indel_count}}
        results = defaultdict(dict)

        for read in bam.fetch():
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue  # Skip unmapped reads

            contig_name = read.reference_name
            read_id = read.query_name

            # Count indels above the threshold in the current read
            indel_count = 0
            for operation, length in read.cigartuples:
                # 1: insertion, 2: deletion
                if operation in (1, 2) and length > args.indel_threshold:
                    indel_count += 1

            # Record the result for this read
            results[contig_name][read_id] = indel_count

    for contig, reads in results.items():
        for read_id, indel_count in reads.items():
            sys.stdout.write(f"{contig}\t{read_id}\t{indel_count}\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("CountInd")
    parser.add_argument(
        "bam_file",
        type=pathlib.Path,
        help="BAM file with indels to count",
    )
    parser.add_argument(
        "indel_threshold",
        type=int,
        help="Length of indels to consider",
    )
    return parser
