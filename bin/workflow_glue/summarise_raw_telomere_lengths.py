"""
Script to process telomere length coverage of raw data.

Then generate summaries and plots.
"""
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    # Load the datasets
    df = pd.read_csv(
        args.telomere_lengths, sep="\t", header=None, names=["Read", "Telomere_length"]
    )

    # Define aggregation functions
    agg_functions = {"Telomere_length": ["mean", "std", "max"]}

    # Aggregate statistics
    summary_stats = df.agg(agg_functions).round(0)
    summary_stats["Read_count"] = len(df)
    read_count = len(df)

    # Create a new summary DataFrame with the desired format
    summary_formatted = pd.DataFrame(
        {
            "Read count": [read_count],
            "Telomere length mean": [summary_stats["Telomere_length"]["mean"]],
            "Telomere length SD": [summary_stats["Telomere_length"]["std"]],
            "Telomere length max": [int(summary_stats["Telomere_length"]["max"])],
        }
    )

    # Save the formatted summary DataFrame to a CSV
    summary_formatted.to_csv(f"{args.output_prefix}_raw_Coverage.csv", index=False)

    # Save per-Read information including telomere length
    df.to_csv(f"{args.output_prefix}_raw_Per_Read_telomere_length.csv", index=False)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("telomere_length_coverage_raw")
    parser.add_argument(
        "telomere_lengths",
        help="tab-separated file with read ids and telomere lengths (no header line)",
    )
    parser.add_argument("output_prefix", help="prefix for output file names")
    return parser
