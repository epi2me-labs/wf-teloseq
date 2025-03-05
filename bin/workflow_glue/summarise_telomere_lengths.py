"""
Script to analyze telomere lengths and coverage data.

It can operate in two modes:
1. 'alignment' mode: Processes telomere lengths along with alignment data.
2. 'raw' mode: Processes telomere lengths from raw data without alignment.

Use the '--mode' command-line argument to select the mode.
"""

import pandas as pd

from .util import wf_parser, get_named_logger  # noqa: ABS101

logger = get_named_logger(__name__)


def parse_seqkit(fname):
    """Parse seqkit output file and return a DataFrame."""
    cols = {
        "Read": str,
        "Ref": str,
        "Acc": float,
        "Strand": int,
        "IsSec": bool,
        "IsSup": bool,
    }
    df = (
        pd.read_csv(fname, sep="\t", dtype=cols, usecols=cols.keys(), index_col="Read")
        .query(
            # only keep primary alignments on + strand
            "not IsSec and not IsSup and Strand == 1"
        )
        .drop(columns=["IsSec", "IsSup", "Strand"])
    )
    return df


def process_alignment_mode(args):
    """Process data in alignment mode with tag filtering, pat/mat arm separation."""
    # Parse input files
    alignment_df = parse_seqkit(args.seqkit_bam)

    # Load telomere lengths with tags
    lengths_df = pd.read_csv(
        args.telomere_lengths,
        sep="\t",
        header=None,
        names=["Read", "Telomere_length", "Tags"],
        dtype={"Read": str, "Telomere_length": float, "Tags": str}
    )

    # Merge alignment data with telomere lengths
    merged_df = alignment_df.merge(lengths_df, left_index=True, right_on="Read")
    merged_df.set_index("Read", inplace=True)

    # Handle missing values in the Tags column
    merged_df["Tags"] = merged_df["Tags"].fillna("")

    # Load unaligned telomere lengths if provided
    unaligned_df = None
    if args.unaligned:
        unaligned_summary = pd.read_csv(
            args.unaligned,
            sep=',',
            dtype={
                'Read count': int,
                'Telomere length mean': float,
                'Telomere length SD': float,
                'Telomere length max': int,
                'Telomere length N50': int,
                'Telomere length CV': float
            }
        )

        # Augment the summary with the tag information
        unaligned_summary['Filter'] = ['Unaligned']
        unaligned_summary['Ref_Type'] = ['UH']

        # Add this to the combined summary
        if 'combined_summary' not in locals():
            combined_summary = []
        combined_summary.append(unaligned_summary)

    # Define tag filters
    tag_filters = ["HS", "LS", "NS"]

    # Define paternal and maternal chromosome/arm identifiers
    def classify_ref(ref):
        ref_lower = ref.lower()
        if 'pat' in ref_lower:
            ref_type = 'Pat'
        elif 'mat' in ref_lower:
            ref_type = 'Mat'
        else:
            ref_type = 'UH'
        return ref_type

    merged_df['Ref_Type'] = merged_df['Ref'].apply(classify_ref)

    combined_summary = []  # Store individual summaries to combine later

    # Process each tag filter
    for tag in tag_filters:
        # Filter the dataframe for reads that contain the specific tag independently
        filtered_df = merged_df[merged_df["Tags"].apply(lambda x: tag in x.split(','))]

        # If unaligned data exists, combine unaligned reads with the same tag
        if unaligned_df is not None:
            filter_condition = unaligned_df["Tags"].apply(lambda x: tag in x.split(','))
            unaligned_tag_df = unaligned_df[filter_condition]

            # Combine aligned and unaligned dataframes
            # Add a column to distinguish between aligned and unaligned
            filtered_df['Alignment_Status'] = 'Aligned'
            unaligned_tag_df['Alignment_Status'] = 'Unaligned'
            filtered_df = pd.concat([filtered_df, unaligned_tag_df], ignore_index=True)

        if filtered_df.empty:
            continue

        # First, generate global summary across UH references
        global_summary = None

        # For the aligned data only, generate per-reference statistics
        if 'Alignment_Status' in filtered_df.columns:
            aligned_df = filtered_df[filtered_df['Alignment_Status'] == 'Aligned']
        else:
            aligned_df = filtered_df

        # Group by reference and Ref_Type and aggregate statistics
        per_ref_summary = aligned_df.groupby(["Ref", "Ref_Type"]).agg(
            **{
                "Coverage": ("Acc", "count"),
                "avg_accuracy": ("Acc", "mean"),
                "SD_accuracy": ("Acc", "std"),
                "Telomere_mean": ("Telomere_length", "mean"),
                "Telomere_sd": ("Telomere_length", "std"),
                "Telomere_max": ("Telomere_length", "max"),
            }
        ).reset_index()

        # Round float columns to two decimals
        columns_to_round = [
            "avg_accuracy",
            "SD_accuracy",
            "Telomere_mean",
            "Telomere_sd"
        ]
        per_ref_summary[columns_to_round] = per_ref_summary[columns_to_round].round(2)

        # Apply minimum coverage filter
        if args.min_coverage is not None:
            coverage_filter = per_ref_summary['Coverage'] >= args.min_coverage
            per_ref_summary = per_ref_summary[coverage_filter]

        # Re-order summary to chr order
        per_ref_summary = per_ref_summary.copy()
        # per_ref_summary['sort_order'] = per_ref_summary['Ref'].apply(sort_key)
        per_ref_summary = per_ref_summary.sort_values(['sort_order', 'Ref']).drop(
            'sort_order', axis=1
        )

        # Save filtered summary to CSV
        per_ref_summary.to_csv(
            f"{args.output_prefix}_{tag}_chr_arm_coverage.csv",
            index=False
        )

        # Reset index to ensure 'Read' is a column, then reorder columns
        if 'Read' not in filtered_df.columns:
            filtered_df = filtered_df.reset_index()

        # Put 'Read' as the first column
        cols = ['Read'] + [c for c in filtered_df.columns if c != 'Read']
        filtered_df = filtered_df[cols]

        # Save per-read information including telomere length to CSV
        filtered_df.to_csv(
            f"{args.output_prefix}_{tag}_per_read_telomere_length.csv",
            index=False
        )

        # Compute summaries for each ref type (paternal, maternal, other)
        type_summaries = [global_summary]  # Start with global summary
        for ref_type in ['Pat', 'Mat', 'UH']:
            ref_type_df = filtered_df[filtered_df['Ref_Type'] == ref_type]

            if ref_type_df.empty:
                continue

            # Generate summary for this specific ref type
            # type_summary = generate_summary(ref_type_df, tag, ref_type)
            # type_summaries.append(type_summary)

        # Combine ref type summaries and save
        if type_summaries:
            type_summary_df = pd.concat(type_summaries, ignore_index=True)
            type_summary_df.to_csv(f"{args.output_prefix}_{tag}.csv", index=False)

            # Only add HC and LC tags to combined summary
            if tag in ['HS', 'LS']:
                # Remove duplicate rows based on all columns - temp fix
                type_summary_df = type_summary_df.drop_duplicates()
                combined_summary.append(type_summary_df)

    # Combine UH individual summaries into one DataFrame
    if combined_summary:
        # Ensure unaligned data is included
        if unaligned_summary is not None:
            combined_summary.append(unaligned_summary)

        combined_summary_df = pd.concat(combined_summary, ignore_index=True)
        combined_summary_df.to_csv(
            f"{args.output_prefix}_combined_summary.csv", index=False
        )


def main(args):
    """Run the entry point."""
    process_alignment_mode(args)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("TeloLen")
    parser.add_argument(
        "--telomere_lengths",
        required=True,
        help="Tab-separated file with read IDs, telomere lengths, and tags (no header)",
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="Prefix for output CSV files",
    )
    parser.add_argument(
        "--unaligned",
        help="Raw telomere coverage data",
    )
    # Arguments specific to 'alignment' mode
    parser.add_argument(
        "--seqkit_bam",
        help="Seqkit BAM output TSV file (required for 'alignment' mode)",
    )
    parser.add_argument(
        "--min_coverage",
        type=int,
        help="Minimum coverage to include contigs (required for 'alignment' mode)",
    )
    return parser
