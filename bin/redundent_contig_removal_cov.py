#!/usr/bin/env python
"""Script to filter contigs on coverage from a BAM file with coverage threshold."""

from collections import defaultdict
import sys

import pysam

# Read MAX_COVERAGE from coverage file
coverage_file_path = sys.argv[2]
with open(coverage_file_path, 'r') as file:
    coverage_value_str = file.readline().strip()  # first line and strip newline
    coverage_value = int(coverage_value_str)

# Convert the string value to float
MAX_COVERAGE = float(coverage_value)


def filter_contigs(bamfile):
    """
    Filter contigs based on coverage and write low coverage contigs to an output file.

    Args:
    bamfile (str): Path to BAM file.
    """
    bam = pysam.AlignmentFile(bamfile, "rb")

    primary_counts = defaultdict(int)
    for refid in bam.references:
        coverage = bam.count(
            contig=refid,
            read_callback=lambda x: (
                not x.is_secondary and not x.is_supplementary and x.mapping_quality > 0
            )
        )

        if coverage > 0:
            primary_counts[refid] = coverage
        else:
            primary_counts[refid] = 0

    low_coverage_contigs = set()
    for c in primary_counts:
        coverage = primary_counts[c]

        if coverage < MAX_COVERAGE:
            low_coverage_contigs.add(c)
            with open('idlisttoremove.txt', 'a') as outfile:
                outfile.write(c + '\n')


if __name__ == '__main__':
    bamfile = sys.argv[1]

    filter_contigs(bamfile)
