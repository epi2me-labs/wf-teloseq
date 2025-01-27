"""Script: filter bam and reference by minimum primary coverage."""

import logging
from pathlib import Path

import pysam

from .util import wf_parser  # noqa: ABS101


def validate_input_files(bam_file, reference_file):
    """Validate that input files exist."""
    if not bam_file.exists() or not reference_file.exists():
        raise FileNotFoundError("Input BAM or reference file does not exist")


def setup_logging():
    """Set up logging for the script."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def filter_contigs(bam_file, min_primary_count, logger):
    """
    Filter contigs based on primary read count.

    Uses direct fetching instead of iterating through all reads.
    """
    contigs_to_keep = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        if not bam.header.get('HD', {}).get('SO') == 'coordinate':
            raise ValueError(
                "BAM file is not coordinate sorted. Please provide a \
                coordinate-sorted BAM file."
            )

        for contig in bam.references:
            primary_reads = sum(
                1 for read in bam.fetch(contig)
                if (not read.is_secondary and
                    not read.is_supplementary and
                    read.mapping_quality > 0)
            )

            if primary_reads >= min_primary_count:
                contigs_to_keep.append(contig)

    return contigs_to_keep


def create_filtered_bam(bam_file, contigs_to_keep, output_dir, logger):
    """Create a filtered BAM file using region-based filtering."""
    output_bam_path = output_dir / "filtered.bam"

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Create a new header with only the contigs to keep
        new_header = bam.header.to_dict()
        new_header['SQ'] = [
            sq for sq in new_header['SQ'] if sq['SN'] in contigs_to_keep
        ]

        # Create a mapping from original reference IDs to new reference IDs
        new_ref_name_to_id = {sq['SN']: idx for idx, sq in enumerate(new_header['SQ'])}

        with pysam.AlignmentFile(output_bam_path, "wb", header=new_header) as out_bam:
            for contig in contigs_to_keep:
                for read in bam.fetch(contig):
                    # Update the read's reference_id to match the new header
                    read.reference_id = new_ref_name_to_id[read.reference_name]
                    out_bam.write(read)

    pysam.index(str(output_bam_path))
    logger.info(f"Filtered BAM file created: {output_bam_path}")


def create_filtered_reference(reference_file, contigs_to_keep, output_dir, logger):
    """Create a filtered reference file."""
    output_reference_path = output_dir / "filtered_reference.fasta"

    with (
        pysam.FastaFile(reference_file) as ref_in,
        open(output_reference_path, "w") as ref_out
    ):
        for contig in contigs_to_keep:
            sequence = ref_in.fetch(contig)
            ref_out.write(f">{contig}\n{sequence}\n")

    logger.info(f"Filtered reference file created: {output_reference_path}")


def run(args):
    """Run the contig filtering process."""
    bam_file = Path(args.bam_file)
    reference_file = Path(args.reference_file)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging()

    validate_input_files(bam_file, reference_file)
    min_primary_count = int(args.min_primary_count)

    contigs_to_keep = filter_contigs(bam_file, min_primary_count, logger)
    create_filtered_bam(bam_file, contigs_to_keep, output_dir, logger)
    create_filtered_reference(reference_file, contigs_to_keep, output_dir, logger)

    logger.info("Contig filtering complete!")


def argparser():
    """Argument parser for entry point."""
    parser = wf_parser("FiltCont")
    parser.add_argument(
        "bam_file",
        help="Input BAM file to filter",
    )
    parser.add_argument(
        "reference_file",
        help="Reference FASTA file",
    )
    parser.add_argument(
        "output_dir",
        help="Directory for output files",
    )
    parser.add_argument(
        "min_primary_count",
        type=int,
        help="Value with minimum primary read count to keep a contig",
    )
    parser.add_argument(
        "threads",
        type=int,
        default=None,
        help="Number of threads to use",
    )
    return parser


def main(args):
    """Run the entry point."""
    run(args)
