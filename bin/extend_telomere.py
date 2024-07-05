#!/usr/bin/env python
"""Extend telomere in FASTA file."""
import sys


def trim_and_extend_sequence(seq, repeat="TAACCC", extension_count=4000):
    """
    Trim and extend sequence.

    Trim the sequence to the first occurrence of 'repeat'. Then extend it
    by adding 'repeat' to the start of the sequence.
    """
    trim_index = seq.find(repeat)
    if trim_index != -1:
        # Trim to the first occurrence of the repeat
        trimmed_seq = seq[trim_index:]
        # Add repeat to the start of the sequence
        extended_seq = repeat * extension_count + trimmed_seq
        return extended_seq
    else:
        # If repeat not found, return the original sequence
        return seq


def process_fasta(input_file, output_file):
    """Process the input FASTA file and save the output."""
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        sequence = ""
        header = ""  # Initialize header variable
        for line in infile:
            if line.startswith(">"):
                if sequence:
                    extended_seq = trim_and_extend_sequence(sequence)
                    outfile.write(header + "\n" + extended_seq + "\n")
                    sequence = ""
                header = line.strip()
            else:
                sequence += line.strip()
        # Process the last sequence in the file
        if sequence:
            extended_seq = trim_and_extend_sequence(sequence)
            outfile.write(header + "\n" + extended_seq + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stdout.write("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_fasta(input_file, output_file)
