"""
Script: BAM Analyzer for Large Indels and Mismatch Clusters.

Analyzes BAM files to identify contigs with well-supported large indels (≥40bp)
and extracts reads from contigs where indels have support threshold or clusters
of mismatches for later reclustering by other scripts.
"""

import logging
from pathlib import Path

import numpy as np
import pysam

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entry point."""
    parser = wf_parser("Mismatch")
    parser.add_argument(
        "bam_file",
        help="Input BAM file"
    )
    parser.add_argument(
        "reference_file",
        help="Reference FASTA file"
    )
    parser.add_argument(
        "output_dir",
        help="Output directory"
    )
    parser.add_argument(
        "support_value",
        help="Minimum support value (single integer)"
    )
    parser.add_argument(
        "min_indel_length",
        help="Minimum indel length to consider",
        type=int,
        default=40
    )
    parser.add_argument(
        "min_mismatches",
        help="Minimum mismatches in window",
        type=int,
        default=40
    )
    parser.add_argument(
        "mismatch_window",
        help="Window size for mismatch clustering",
        type=int,
        default=1000
    )
    parser.add_argument(
        "secondary_min_mismatches",
        help="Minimum mismatches in secondary window",
        type=int,
        default=30
    )
    parser.add_argument(
        "secondary_mismatch_window",
        help="Secondary window size for mismatch clustering",
        type=int,
        default=600
    )
    return parser


def main(args):
    """Run the entry point."""
    analyze_bam(
        args.bam_file,
        args.reference_file,
        args.output_dir,
        min_support=int(args.support_value),
        min_indel_length=args.min_indel_length,
        min_mismatches=args.min_mismatches,
        mismatch_window=args.mismatch_window,
        secondary_min_mismatches=args.secondary_min_mismatches,
        secondary_mismatch_window=args.secondary_mismatch_window,
    )


def analyze_bam(
    bam_file,
    reference_file,
    output_dir,
    min_support,
    min_indel_length=40,
    min_mismatches=40,
    mismatch_window=1000,
    secondary_min_mismatches=30,
    secondary_mismatch_window=600,
):
    """Analyzes BAM file to identify large indels and mismatch clusters."""
    logger = logging.getLogger(__name__)

    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting analysis...")

    summary_records, contigs_with_support = analyze_indels_and_mismatches(
        bam_file,
        reference_file,
        min_support,
        min_indel_length,
        min_mismatches,
        mismatch_window,
        secondary_min_mismatches,
        secondary_mismatch_window,
        output_dir,
        logger
    )

    logger.info(f"Found {len(contigs_with_support)} contigs supported")

    # Write summary
    write_summary(
        output_dir,
        min_indel_length,
        min_support,
        min_mismatches,
        mismatch_window,
        secondary_min_mismatches,
        secondary_mismatch_window,
        summary_records,
        contigs_with_support,
    )

    logger.info("Analysis complete!")


def analyze_indels_and_mismatches(
    bam_file,
    reference_file,
    min_support,
    min_indel_length,
    min_mismatches,
    mismatch_window,
    secondary_min_mismatches,
    secondary_mismatch_window,
    output_dir,
    logger
):
    """Process all contigs to find indels and mismatch clusters."""
    logger.info("Analyzing contigs for large indels and mismatch clusters...")

    bam = pysam.AlignmentFile(bam_file, "rb")
    reference = pysam.FastaFile(reference_file)
    contigs = list(bam.references)

    # Simple accumulators:
    summary_records = []
    contigs_with_support = set()

    for contig in contigs:
        try:
            (
                contig_records,
                contig_support,
                reads_to_extract
            ) = find_indels_and_mismatches_in_contig(
                bam,
                reference,
                contig,
                min_support,
                min_indel_length,
                min_mismatches,
                mismatch_window,
                secondary_min_mismatches,
                secondary_mismatch_window,
                logger
            )
            # If the contig has support, extract relevant reads
            if contig_support:
                contigs_with_support.add(contig)
                extract_reads_from_contig(
                    bam, output_dir, contig, reads_to_extract, logger
                )
            # Accumulate summary records
            summary_records.extend(contig_records)

        except Exception as e:
            logger.error(f"Error processing contig {contig}: {e}")

    return summary_records, contigs_with_support


def find_indels_and_mismatches_in_contig(
    bam,
    reference,
    contig,
    min_support,
    min_indel_length,
    min_mismatches,
    mismatch_window,
    secondary_min_mismatches,
    secondary_mismatch_window,
    logger
):
    """Identify large indels and mismatch clusters in a contig."""
    from pysam import (
        CMATCH,
        CINS,
        CDEL,
        CREF_SKIP,
        CSOFT_CLIP,
        CEQUAL,
        CDIFF,
    )

    ref_seq = reference.fetch(contig)

    # Temporary structures to hold indel and mismatch cluster info for this contig
    primary_mismatch_supporting_reads = set()
    secondary_mismatch_supporting_reads = set()
    large_indel_reads = set()

    # Predefine sets for quick membership checks
    match_ops = {CMATCH, CEQUAL, CDIFF}

    for read in bam.fetch(contig):
        if read.is_secondary or read.is_supplementary:
            continue

        if read.reference_end is None:
            continue

        ref_start = read.reference_start
        ref_end = read.reference_end
        ref_length = ref_end - ref_start
        if ref_length <= 0:
            continue

        # Collect errors
        errors = np.zeros(ref_length, dtype=int)
        read_pos = 0
        ref_pos = ref_start

        # Track large indels supported by this read
        this_read_large_indels = False

        for operation, length in (read.cigartuples or []):
            if operation in match_ops:
                read_seq = read.query_sequence[read_pos:read_pos + length]
                ref_seq_segment = ref_seq[ref_pos:ref_pos + length]

                if operation == CEQUAL:
                    # Perfect matches, no errors
                    pass
                elif operation == CDIFF:
                    # All mismatches
                    errors[(ref_pos - ref_start):(ref_pos - ref_start + length)] = 1
                else:  # CMATCH
                    # Check base by base
                    for i, (r_base, q_base) in enumerate(
                        zip(ref_seq_segment, read_seq)
                    ):
                        if r_base.upper() != q_base.upper():
                            errors[(ref_pos - ref_start) + i] = 1

                read_pos += length
                ref_pos += length

            elif operation == CINS:
                # Insertion in read
                if length >= min_indel_length:
                    this_read_large_indels = True
                    # Mark an error at the current ref position
                if ref_start <= ref_pos < ref_end:
                    errors[ref_pos - ref_start] = 1
                read_pos += length

            elif operation == CDEL:
                # Deletion from reference
                if length >= min_indel_length:
                    this_read_large_indels = True
                # Mark deleted region as errors
                errors[(ref_pos - ref_start):(ref_pos - ref_start + length)] = 1
                ref_pos += length

            elif operation == CSOFT_CLIP:
                # Soft clip: no ref advance, just move read pos
                read_pos += length

            elif operation == CREF_SKIP:
                # Skipped region on the reference
                ref_pos += length

        # Check clusters after processing the read
        if has_cluster_of_errors(errors, mismatch_window, min_mismatches):
            primary_mismatch_supporting_reads.add(read.query_name)

        if has_cluster_of_errors(
            errors,
            secondary_mismatch_window,
            secondary_min_mismatches
        ):
            secondary_mismatch_supporting_reads.add(read.query_name)

        if this_read_large_indels:
            large_indel_reads.add(read.query_name)

    # Determine support for this contig
    # If contig has ≥min_support reads with large indels, or mismatch clusters
    contig_has_support = False
    if len(large_indel_reads) >= min_support:
        contig_has_support = True
    if len(primary_mismatch_supporting_reads) >= min_support:
        contig_has_support = True
    if len(secondary_mismatch_supporting_reads) >= min_support:
        contig_has_support = True

    # Prepare summary records
    contig_records = []

    # Indel summary (if contig is supported)
    if contig_has_support and len(large_indel_reads) > 0:
        contig_records.append({
            "Contig": contig,
            "Position": "NA",
            "Type": "Indel_Cluster",
            "Length": "NA",
            "Supporting_Reads": len(large_indel_reads),
            "Cluster_Type": "Indel"
        })

    # Mismatch cluster summary
    if contig_has_support and len(primary_mismatch_supporting_reads) > 0:
        contig_records.append({
            "Contig": contig,
            "Position": "NA",
            "Type": "Mismatch_Cluster",
            "Length": "NA",
            "Supporting_Reads": len(primary_mismatch_supporting_reads),
            "Cluster_Type": "primary"
        })

    if contig_has_support and len(secondary_mismatch_supporting_reads) > 0:
        contig_records.append({
            "Contig": contig,
            "Position": "NA",
            "Type": "Mismatch_Cluster",
            "Length": "NA",
            "Supporting_Reads": len(secondary_mismatch_supporting_reads),
            "Cluster_Type": "secondary"
        })

    # Reads to extract if contig is supported
    reads_to_extract = (
        large_indel_reads
        .union(primary_mismatch_supporting_reads)
        .union(secondary_mismatch_supporting_reads)
    )

    # Write read IDs if supported
    if contig_has_support:
        reads_file = Path(reference.filename.decode()).parent / f"{contig}_read_ids.txt"
        with open(reads_file, 'w') as f_reads:
            for read_id in sorted(reads_to_extract):
                f_reads.write(f"{read_id}\n")

    return contig_records, contig_has_support, reads_to_extract


def has_cluster_of_errors(errors, window_length, cluster_size):
    """Check if there is at least one window of size window_length.

    with at least cluster_size errors using convolution.
    """
    if len(errors) == 0:
        return False
    if window_length > len(errors):
        # If window is larger than the read span, just sum everything
        return np.sum(errors) >= cluster_size
    window_sum = np.convolve(errors, np.ones(window_length, dtype=int), mode='valid')
    return np.any(window_sum >= cluster_size)


def extract_reads_from_contig(bam, output_dir, contig, reads_to_extract, logger):
    """Extract reads from contig with well-supported variants."""
    logger.info(f"Extracting reads from contig {contig}...")

    output_path = output_dir / f"{contig}_reads.bam"
    with pysam.AlignmentFile(output_path, "wb", template=bam) as out_bam:
        try:
            for read in bam.fetch(contig):
                if (
                    read.query_name in reads_to_extract and
                    not (read.is_secondary or read.is_supplementary)
                ):
                    out_bam.write(read)
        except ValueError as e:
            logger.error(f"Error extracting reads from contig {contig}: {e}")

    # Index the contig BAM file
    pysam.index(str(output_path))


def write_summary(
    output_dir,
    min_indel_length,
    min_support,
    min_mismatches,
    mismatch_window,
    secondary_min_mismatches,
    secondary_mismatch_window,
    summary_records,
    contigs_with_support,
):
    """Write analysis summary and read IDs to files."""
    summary_path = output_dir / "indel_summary.csv"

    with open(summary_path, 'w') as f:
        # Write header comments
        f.write("# Large Indel and Mismatch Cluster Analysis Summary\n")
        f.write(f"# Minimum indel length: {min_indel_length}bp\n")
        f.write(f"# Minimum support threshold: {min_support}\n")
        f.write(
            f"# Minimum mismatches in primary window "
            f"({mismatch_window}bp): {min_mismatches}\n")
        f.write(
            f"# Minimum mismatches in secondary window "
            f"({secondary_mismatch_window}bp): "
            f"{secondary_min_mismatches}\n")
        f.write(
            f"# Contigs with well-supported variants: "
            f"{len(contigs_with_support)}\n")

        # Write CSV header
        f.write("Contig,Position,Type,Length,Supporting_Reads,Cluster_Type\n")

        for record in summary_records:
            f.write(
                f"{record['Contig']},"
                f"{record['Position']},"
                f"{record['Type']},"
                f"{record['Length']},"
                f"{record['Supporting_Reads']},"
                f"{record['Cluster_Type']}\n"
            )


if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)
