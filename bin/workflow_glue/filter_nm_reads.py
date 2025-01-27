"""Script: filter reads based upon proportion of NM should be below 4% user defined."""

import sys

import pysam

from .util import wf_parser  # noqa: ABS101

# Using pysam constants:
# pysam.CMATCH (0), pysam.CEQUAL (7), pysam.CDIFF (8)
MATCH_OPS = {pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF}


def argparser():
    """Argument parser for filter_NM_reads."""
    parser = wf_parser("FilterNm")
    parser.add_argument("input", help="Input BAM file")
    parser.add_argument("output", help="Output BAM file")
    parser.add_argument(
        "percent", type=float,
        help="Maximum allowed NM percentage (e.g., 4 for 4%)")
    return parser


def match_count(read):
    """Return the number of matched positions in the read."""
    return sum(length for op, length in (read.cigartuples or []) if op in MATCH_OPS)


def process_bam(input_bam, output_bam, nm_percent_threshold):
    """Process BAM file and filter reads based on NM percentage threshold."""
    try:
        bamfile = pysam.AlignmentFile(input_bam, "rb")
    except (OSError, ValueError) as e:
        sys.exit(f"Error opening input BAM file: {e}")

    try:
        filtered_bam = pysam.AlignmentFile(output_bam, "wb", header=bamfile.header)
    except (OSError, ValueError) as e:
        bamfile.close()
        sys.exit(f"Error creating output BAM file: {e}")

    total_reads = 0
    kept_reads = 0

    for read in bamfile:
        total_reads += 1
        if read.is_unmapped:
            continue

        try:
            nm_value = read.get_tag("NM")
            matched_length = match_count(read)
            if matched_length > 0:
                nm_percentage = (nm_value / matched_length) * 100
                if nm_percentage <= nm_percent_threshold:
                    filtered_bam.write(read)
                    kept_reads += 1
        except KeyError:
            continue

    bamfile.close()
    filtered_bam.close()


def main(args):
    """Process to filter bam file by NM."""
    process_bam(args.input, args.output, args.percent)
