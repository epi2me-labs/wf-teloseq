"""Script to filter contigs based on sequence length, GC content, and telomere ends."""

import re

import pysam

from .util import wf_parser  # noqa: ABS101

# Constants
MIN_LENGTH_FOR_GC = 30000  # Minimum length for GC content check (30kbp)
TEL_PATTERN = re.compile(r"^(CTAACC)+|((CTAACC)+)$")  # Precompiled telomere regex
SEQ_PATTERN = re.compile(r'seqs=(\d+)')  # Precompiled seqs pattern regex
GC_LOWER = 41.1
GC_UPPER = 41.4


def is_telomere_sequence(sequence):
    """Check if the sequence is a repetitive telomere sequence."""
    telomere_pattern = "CTAACC"
    repeat_count = len(sequence) // len(telomere_pattern)
    return sequence == telomere_pattern * repeat_count


def calculate_gc_content(sequence):
    """Calculate GC content as a percentage."""
    gc_count = sequence.count('G') + sequence.count('C')
    total_length = len(sequence)
    return (gc_count / total_length) * 100 if total_length > 0 else 0


def main(args):
    """Run the entry point."""
    qualifying_contigs = []
    filtered_records = []

    # Variable for min_seqs value
    min_seqs = int(args.min_seqs)

    # Process input file using Pysam
    with pysam.FastxFile(args.input_file) as infile:
        for entry in infile:
            sequence = entry.sequence
            header = entry.name

            # Skip pure telomere sequences
            if not is_telomere_sequence(sequence):
                # Remove tel repeats from both ends
                processed_seq = TEL_PATTERN.sub('', sequence)
                sequence_length = len(processed_seq)

                # Extract seqs value from header
                match = SEQ_PATTERN.search(header)
                seqs_value = int(match.group(1)) if match else 0

                # Determine if the record qualifies based on length and GC content
                if sequence_length >= MIN_LENGTH_FOR_GC:
                    gc_content = calculate_gc_content(processed_seq)
                    if GC_LOWER <= gc_content <= GC_UPPER:
                        qualifying_contigs.append((header, sequence))
                        if seqs_value >= min_seqs:
                            filtered_records.append((header, sequence))
                    elif seqs_value >= min_seqs:
                        filtered_records.append((header, sequence))
                elif seqs_value >= min_seqs:
                    filtered_records.append((header, sequence))

    # Determine whether to apply min_seqs filter based on number of qualifying contigs
    apply_min_seqs = len(qualifying_contigs) > 2

    # Write output file
    with open(args.output_file, "w") as outfile:
        for header, sequence in filtered_records:
            if (
                not apply_min_seqs or
                (
                    'seqs=' in header and
                    int(SEQ_PATTERN.search(header).group(1)) >= min_seqs
                )
            ):
                # Write in FASTA format
                outfile.write(f">{header}\n{sequence}\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("FiltCont")
    parser.add_argument(
        "input_file",
        help="Input multi-FASTA file",
    )
    parser.add_argument(
        "output_file",
        help="Output file to save filtered contigs",
    )
    parser.add_argument(
        "min_seqs",
        help="Minimum seqs value",
    )
    return parser
