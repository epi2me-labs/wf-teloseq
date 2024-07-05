#!/usr/bin/env python
"""Script to filter contigs based on coverage from a BAM file."""

from collections import defaultdict
import sys

import pysam

# Exclude below 20% of max coverage done in coverage process
max_coverage_value = float(sys.argv[2])
MAX_COVERAGE = round(max_coverage_value)


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
            read_callback=lambda x: not x.is_secondary
            and not x.is_supplementary
            and x.mapping_quality > 0
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
