#!/usr/bin/env python
"""Script for deduplicating multi fasta by sequence."""

import sys


def contains_homopolymer(sequence, length=5):
    """Check if the sequence contains any homopolymer of given length or more."""
    for i in range(len(sequence) - length + 1):
        if len(set(sequence[i:i + length])) == 1:
            return True
    return False


def remove_homopolymers(sequence, length=5):
    """
    Temporarily remove homopolymers of a specified length or more.

    The original sequence is not modified; this is for comparison only.
    """
    result = []
    i = 0
    while i < len(sequence):
        if i <= len(sequence) - length and len(set(sequence[i:i + length])) == 1:
            i += length
            while i < len(sequence) and sequence[i] == sequence[i - 1]:
                i += 1
        else:
            result.append(sequence[i])
            i += 1
    return ''.join(result)


def normalize_sequence(
        sequence, telomere_sequence="TAACCCTAACCCTAACCCTAACCCTAACCCTAACCC"):
    """Normalize sequence by removing telomere sequence if present at the end."""
    last_occurrence = sequence.rfind(telomere_sequence)
    if last_occurrence != -1:
        return sequence[last_occurrence + len(telomere_sequence):]
    return sequence


output_file = sys.argv[1]
sequences = []

current_header = None  # Store the current sequence header
current_seq = ''  # Initialize current_seq
output_lines = []  # Store output lines to write to file

for line in sys.stdin:
    line = line.strip()

    if line.startswith('>'):
        if current_header is not None:
            # Handle case where previous sequence had no duplicates but wasn't printed
            output_lines.append(current_header)
            output_lines.append(current_seq)
        current_header = line
        current_seq = ''  # Reset current sequence
    else:
        current_seq += line  # Accumulate sequence lines
        seq_trimmed = current_seq[80:-15]
        seq_trimmed_no_homopolymers = remove_homopolymers(seq_trimmed)
        seq_trimmed_no_homopolymers_notelomere = (
            normalize_sequence(seq_trimmed_no_homopolymers)
        )

        is_duplicate = False

        for saved_header, saved_seq in sequences:
            saved_seq_trimmed = saved_seq[80:-15]
            saved_seq_trimmed_no_homopolymers = remove_homopolymers(saved_seq_trimmed)
            saved_seq_trimmed_no_homopolymers_notelomere = (
                normalize_sequence(saved_seq_trimmed_no_homopolymers)
            )

            if (
                seq_trimmed_no_homopolymers_notelomere
                in saved_seq_trimmed_no_homopolymers_notelomere
                or saved_seq_trimmed_no_homopolymers_notelomere
                in seq_trimmed_no_homopolymers_notelomere
            ):

                is_duplicate = True
                break

        if not is_duplicate:
            if current_header is not None:
                # Save header if not duplicate
                output_lines.append(current_header)
                # Reset current header to avoid reprinting
                current_header = None
            # Save original sequence
            output_lines.append(current_seq)
            sequences.append((current_header, current_seq))

# Handle the last sequence in the file
if current_header is not None and not is_duplicate:
    output_lines.append(current_header)
    output_lines.append(current_seq)

with open(output_file, 'w') as f:
    for line in output_lines:
        f.write(f"{line}\n")
