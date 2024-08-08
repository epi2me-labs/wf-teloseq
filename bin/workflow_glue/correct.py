"""
This script reads a VCF file and a reference FASTA file.

Applies the variants to create consensus sequences,
and writes the output to a specified file.
It processes the VCF file for 1/1 genotypes and applies
the changes to all contigs.
"""

# TODO: looks like we could do this with `bcftools consensus`?

import vcf

from .util import wf_parser  # noqa: ABS101


def apply_variant(ref_seq, position, ref_allele, alt_allele):
    """
    Apply a variant to the reference sequence.

    Parameters:
    ref_seq (str): The reference sequence.
    position (int): The position of the variant.
    ref_allele (str): The reference allele.
    alt_allele (str): The alternate allele.

    Returns:
    str: The sequence with the variant applied.
    """
    if len(ref_allele) == len(alt_allele):
        # SNP
        return ref_seq[:position] + alt_allele + ref_seq[position + len(ref_allele):]
    elif len(alt_allele) > len(ref_allele):
        # Insertion
        insertion_seq = alt_allele[len(ref_allele):]
        new_seq = ref_seq[:position + len(ref_allele)] + insertion_seq
        return new_seq + ref_seq[position + len(ref_allele):]
    elif len(ref_allele) > len(alt_allele):
        # Deletion
        return ref_seq[:position] + ref_seq[position + len(ref_allele):]
    return ref_seq


def main(args):
    """Run the entry point."""
    # Reading the FASTA file
    with open(args.fasta_file) as f:
        fasta_data = f.read().split(">")[1:]
        reference_dict = {
            seq.split("\n", 1)[0]: seq.split("\n", 1)[1].replace("\n", "")
            for seq in fasta_data
        }

    # Processing the VCF file for 1/1 genotypes
    with open(args.vcf_file, 'r') as infile:
        vcf_reader = vcf.Reader(infile)

        # Apply 1/1 genotype changes to all contigs
        for record in vcf_reader:
            ref_name = record.CHROM
            if ref_name in reference_dict:
                ref_seq = reference_dict[ref_name]
                for sample in record.samples:
                    if sample.data.GT == '1/1':
                        position = record.POS - 1
                        ref_allele = str(record.REF)
                        alt_allele = str(record.ALT[0])
                        ref_seq = apply_variant(
                            ref_seq, position, ref_allele, alt_allele
                        )
                reference_dict[ref_name] = ref_seq

    # Writing the consensus sequences to the output file, including unchanged contigs
    with open(args.output_file, 'w') as out_file:
        for ref_name, ref_seq in reference_dict.items():
            out_file.write(f">{ref_name}\n{ref_seq}\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("correct")
    parser.add_argument(
        "vcf_file",
        help="VCF with variants to apply",
    )
    parser.add_argument(
        "fasta_file",
        help="FASTA file with sequence which the variants should be applied to",
    )
    parser.add_argument(
        "output_file",
        help="output FASTA filename",
    )
    return parser
