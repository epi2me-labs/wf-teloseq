#!/usr/bin/env python
"""
This module processes a BAM file to count indels.

(insertions and deletions) above a certain threshold.
"""

import pysam

# Define the threshold for indels
indel_threshold = 9

# Open the BAM file
bam_file = 'alignment.bam'  # Specify your BAM file path here
bam = pysam.AlignmentFile(bam_file, "rb")

# Dictionary to hold results: {contig_name: {read_id: indel_count}}
results = {}

for read in bam.fetch():
    if read.is_unmapped or read.is_supplementary or read.is_secondary:
        continue  # Skip unmapped reads

    contig_name = read.reference_name
    read_id = read.query_name

    # Initialize contig in results if not present
    if contig_name not in results:
        results[contig_name] = {}

    # Count indels above the threshold in the current read
    indel_count = 0
    for operation, length in read.cigartuples:
        # 1: insertion, 2: deletion
        if operation in (1, 2) and length > indel_threshold:
            indel_count += 1

    # Record the result for this read
    results[contig_name][read_id] = indel_count

bam.close()

# Example of how to print or save the results
with open('indel_counts.txt', 'w') as f:
    for contig, reads in results.items():
        for read_id, indel_count in reads.items():
            f.write(f"{contig}\t{read_id}\t{indel_count}\n")
