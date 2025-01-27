"""Add multiple QC tags to BAM file reads."""

import sys

import pandas as pd
import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def tag_bam_reads(input_bam, output_bam, tag_files):
    """
    Add multiple QC tags to reads in a BAM file based on multiple tag files.

    Args:
        input_bam (str): Path to input BAM file
        output_bam (str): Path to output BAM file
        tag_files (list): List of tuples (tag_file_path, tag_value)
    """
    logger = get_named_logger("bamMultiTagging")

    # Combine tags from multiple files
    tags_dict = {}
    try:
        for tag_file, tag_value in tag_files:
            current_tags_df = pd.read_csv(
                tag_file, sep='\t', header=None, names=['read_id', 'existing_tag']
            )

            # Group by read_id and aggregate tags
            for read_id, group in current_tags_df.groupby('read_id'):
                if read_id not in tags_dict:
                    tags_dict[read_id] = set()
                tags_dict[read_id].add(tag_value)
    except Exception as e:
        logger.error(f"Error reading tag files: {e}")
        sys.exit(1)

    # Open input and output BAM files
    try:
        with (
            pysam.AlignmentFile(input_bam, 'rb') as infile,
            pysam.AlignmentFile(output_bam, 'wb', template=infile) as outfile
        ):

            for read in infile:
                # Add multiple QC tags if read ID exists in tags dictionary
                if read.query_name in tags_dict:
                    # Combine multiple tags
                    combined_tags = ','.join(sorted(tags_dict[read.query_name]))
                    read.set_tag('ZQ', combined_tags)

                outfile.write(read)

        # Index the output BAM file
        pysam.index(output_bam)
        logger.info(f"Successfully tagged reads and indexed {output_bam}")

    except Exception as e:
        logger.error(f"Error processing BAM file: {e}")
        sys.exit(1)


def main(args):
    """Run the entry point."""
    # Prepare tag files with their corresponding tag values
    tag_files = [
        (args.strict_tags, 'HS'),      # High Confidence
        (args.lenient_tags, 'LS'),     # Lenient Confidence
        (args.no_filter_tags, 'NS')    # No Confidence
    ]

    tag_bam_reads(
        input_bam=args.input_bam,
        output_bam=args.output_bam,
        tag_files=tag_files
    )


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("tag_bam_reads")
    parser.add_argument(
        "input_bam",
        help="Input BAM file to tag reads",
    )
    parser.add_argument(
        "output_bam",
        help="Output BAM file with added tags",
    )
    parser.add_argument(
        "strict_tags",
        help="TSV file with read IDs for strict filtering",
    )
    parser.add_argument(
        "lenient_tags",
        help="TSV file with read IDs for lenient filtering",
    )
    parser.add_argument(
        "no_filter_tags",
        help="TSV file with read IDs for no filtering",
    )
    return parser
