#!/usr/bin/env python
"""Extract reference from FASTA file."""
import argparse
import sys

from Bio import Seq
from Bio import SeqIO


def main(fasta_file, telomere_sequence, output_file, restriction_cut):
    """Read FASTA file and try to find telomeres."""
    with open(output_file, "w") as output:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq).upper()

            last_3000_bases = sequence[1000:]
            reverse_sequence = Seq.reverse_complement(sequence)
            reverse_last_3000_bases = reverse_sequence[1000:]

            forward_position = -1
            reverse_position = -1

            if telomere_sequence in last_3000_bases:
                forward_position = last_3000_bases.find(telomere_sequence)

            if telomere_sequence in reverse_last_3000_bases:
                reverse_position = reverse_last_3000_bases.find(telomere_sequence)

            forward_telomere_pos = forward_position
            reverse_telomere_pos = reverse_position

            if forward_telomere_pos != -1:
                gatatc_pos = sequence.find(restriction_cut)

                if gatatc_pos != -1:
                    end_pos = min(gatatc_pos + 300, len(sequence))
                    if record.id[-2:].lower() not in ["_q", "_p"]:
                        corrected_header = record.id + "_p"
                    else:
                        corrected_header = record.id
                    extracted_sequence = sequence[1:end_pos]
                    output.write(f">{corrected_header}\n{extracted_sequence}\n")
            else:
                sys.stdout.write(f"No forward telomere found for {record.id}.")
            if reverse_telomere_pos != -1:
                gatatc_pos = reverse_sequence.find(restriction_cut)

                if gatatc_pos != -1:
                    end_pos = min(gatatc_pos + 300, len(reverse_sequence))
                    if record.id[-2:].lower() not in ["_q", "_p"]:
                        corrected_header = record.id + "_q"
                    else:
                        corrected_header = record.id
                    extracted_sequence = reverse_sequence[1:end_pos]
                    output.write(f">{corrected_header}\n{extracted_sequence}\n")
            else:
                sys.stdout.write(f"No reverse telomere found for {record.id}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file", help="Input reference FASTA file")
    parser.add_argument("restriction_cut", help="restriction_cut_site")
    args = parser.parse_args()

    fasta_file = args.fasta_file
    restriction_cut = args.restriction_cut

    telomere_sequence = "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC"
    output_file = "reference2.fasta"  # Replace with the desired output file name
    main(fasta_file, telomere_sequence, output_file, restriction_cut)
