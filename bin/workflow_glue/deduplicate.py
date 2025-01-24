"""Script for deduplicating multi fasta by sequence."""

import sys

from .util import wf_parser  # noqa: ABS101


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


def main(args):
    """Run the entry point."""
    sequences = []

    current_header = None  # Store the current sequence header
    current_seq = ''  # Initialize current_seq

    for line in sys.stdin:
        line = line.strip()

        if line.startswith('>'):
            if current_header is not None:
                # Handle case where previous sequence had no duplicates but
                # wasn't printed
                sys.stdout.write(current_header + "\n")
                sys.stdout.write(current_seq + "\n")
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
                saved_seq_trimmed_no_homopolymers = remove_homopolymers(
                    saved_seq_trimmed
                )
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
                    sys.stdout.write(current_header + "\n")
                    # Reset current header to avoid reprinting
                    current_header = None
                # Save original sequence
                sys.stdout.write(current_seq + "\n")
                sequences.append((current_header, current_seq))

    # Handle the last sequence in the file
    if current_header is not None and not is_duplicate:
        sys.stdout.write(current_header + "\n")
        sys.stdout.write(current_seq + "\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("deduplicate")
    # no args since we're reading from STDIN and writing to STDOUT
    return parser
