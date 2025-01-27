"""
Script to process telomere length coverage of raw data.

Then generate summaries and plots.
"""
import numpy as np
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

    # Calculate coefficient of variation (CV)
    telomere_mean = round(summary_stats["Telomere_length"]["mean"], 0)
    telomere_std = round(summary_stats["Telomere_length"]["std"], 0)
    telomere_cv = (telomere_std / telomere_mean) if telomere_mean != 0 else 0
    telomere_max = summary_stats["Telomere_length"]["max"]

    # Calculate N50 for telomere length using original dataset
    sorted_lengths = df["Telomere_length"].sort_values(ascending=False).values
    cumulative_sum = np.cumsum(sorted_lengths)
    half_sum = cumulative_sum[-1] / 2
    telomere_n50 = sorted_lengths[np.where(cumulative_sum >= half_sum)[0][0]]

    # Create a new summary DataFrame with the desired format
    phenotype = "TERT+" if telomere_cv < 0.7 else "ALT+"

    # Create a new summary DataFrame with the desired format
    summary_formatted = pd.DataFrame(
        {
            "Read count": [read_count],
            "Telomere length mean": telomere_mean,
            "Telomere length SD": telomere_std,
            "Telomere length max": [int(telomere_max)],
            "Telomere length N50": [telomere_n50],
            "Telomere length CV": [round(telomere_cv, 2)],
            "Phenotype": [phenotype],
        }
    )

    # Save the formatted summary DataFrame to a CSV
    summary_formatted.to_csv(f"{args.output_prefix}_raw_coverage.csv", index=False)

    # Save per-Read information including telomere length
    df.to_csv(f"{args.output_prefix}_raw_per_read_telomere_length.csv", index=False)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("telomere_length_coverage_raw")
    parser.add_argument(
        "telomere_lengths",
        help="tab-separated file with read ids and telomere lengths (no header line)",
    )
    parser.add_argument("output_prefix", help="prefix for output file names")
    return parser
