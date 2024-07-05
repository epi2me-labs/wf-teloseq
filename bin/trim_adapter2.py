#!/usr/bin/env python
"""This script trims adapter sequences from the 5' end of sequences in a FASTQ file."""

import logging
import re
import sys

from Bio import SeqIO


def find_earliest_match_and_modify(fastq_path, output_path, search_sequences):
    """
    Trims adapter sequences from the 5' end of sequences in the FASTQ file.

    Then writes the modified sequences to the output file.

    Parameters:
    fastq_path (str): Path to the input FASTQ file.
    output_path (str): Path to the output FASTQ file.
    search_sequences (list): List of sequences to search for trimming.
    """
    modified_records = []

    for record in SeqIO.parse(fastq_path, "fastq"):
        earliest_end_pos = float('inf')

        # Consider only the 5' end and first 250 bp
        sequence_fragment = str(record.seq[:250])

        for search_seq in search_sequences:
            search_seq_clean = re.sub(r'\W+', '', search_seq)
            match = re.search(search_seq_clean, sequence_fragment)
            if match and match.end() < earliest_end_pos:
                earliest_end_pos = match.end()

        if earliest_end_pos != float('inf'):
            # Modify the sequence by removing from 1 to the (end position - 15)
            # Also, trim the quality scores accordingly
            trim_pos = max(0, earliest_end_pos - 14)
            modified_seq = record.seq[trim_pos:]
            modified_qual = record.letter_annotations["phred_quality"][trim_pos:]

            modified_record = record[:]
            # Clear existing annotations before modifying the sequence
            modified_record.letter_annotations = {}
            modified_record.seq = modified_seq
            # Reassign the trimmed quality scores
            modified_record.letter_annotations["phred_quality"] = modified_qual

            modified_records.append(modified_record)

    # Write the modified sequences to an output file
    with open(output_path, "w") as output_handle:
        SeqIO.write(modified_records, output_handle, "fastq")


# Define the sequences to search for
search_sequences = [
    "CACCCTAA​CCCTAACCCTAACC",
    "CAACCCTA​ACCCTAACCCTAAC",
    "CCTAACCC​TAACCCTAACCCTA",
    "CCCTAACC​CTAACCCTAACCCT",
    "CTAACCCT​AACCCTAACCCTAA",
    "CCCCTAAC​CCTAACCCTAACCC"
]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        logging.error("Usage: python script.py <input_fastq> <output_fastq>")
        sys.exit(1)

    fastq_path = sys.argv[1]
    output_path = sys.argv[2]

    find_earliest_match_and_modify(fastq_path, output_path, search_sequences)
