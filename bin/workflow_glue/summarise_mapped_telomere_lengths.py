#!/usr/bin/env python
"""Script to analyze telomere lengths and coverage data."""
import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


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


def main(args):
    """Run the entry point."""
    # Parse input files
    alignment_df = parse_seqkit(args.seqkit_bam)
    lengths = pd.read_csv(
        args.telomere_lengths,
        sep="\t",
        header=None,
        index_col="Read",
        names=["Read", "Telomere_length"],
    )

    alignment_df["Telomere_length"] = lengths

    # Group by reference and aggregate statistics
    per_ref_summary = alignment_df.groupby("Ref").agg(
        **{
            "Coverage": ("Acc", "count"),
            "avg_accuracy": ("Acc", "mean"),
            "SD_accuracy": ("Acc", "std"),
            "Telomere_mean": ("Telomere_length", "mean"),
            "Telomere_sd": ("Telomere_length", "std"),
            "Telomere_max": ("Telomere_length", "max"),
        }
    )

    # Round float columns to two decimals
    per_ref_summary[["avg_accuracy", "SD_accuracy", "Telomere_mean", "Telomere_sd"]] = (
        per_ref_summary[
            ["avg_accuracy", "SD_accuracy", "Telomere_mean", "Telomere_sd"]
        ].round(2)
    )

    # drop refs with low coverage from summary (and the corresponding reads)
    per_ref_summary = per_ref_summary[per_ref_summary["Coverage"] >= args.min_coverage]
    alignment_df = alignment_df.loc[alignment_df["Ref"].isin(per_ref_summary.index)]

    # read ref `.fai` index to get order of reference sequences
    refs = list(pd.read_csv(args.ref_fai, sep="\t", usecols=[0], header=None).squeeze())

    # re-order summary to match ref FASTA
    per_ref_summary = per_ref_summary.loc[
        [ref for ref in refs if ref in per_ref_summary.index]
    ]

    # Save filtered summary to CSV
    per_ref_summary.to_csv(f"{args.output_prefix}_chr_arm_coverage.csv")

    # Save per-read information including telomere length to CSV
    alignment_df.to_csv(f"{args.output_prefix}_per_read_telomere_length.csv")

    # Calculate coefficient of variation (CV)
    telomere_mean = round(alignment_df["Telomere_length"].mean(), 0)
    telomere_std = round(alignment_df["Telomere_length"].std(), 0)
    telomere_cv = (telomere_std / telomere_mean) if telomere_mean != 0 else 0
    telomere_max = int(alignment_df["Telomere_length"].max())

    # Calculate N50 for telomere length
    sorted_lengths = alignment_df["Telomere_length"].sort_values(ascending=False).values
    cumulative_sum = np.cumsum(sorted_lengths)
    half_sum = cumulative_sum[-1] / 2
    telomere_n50 = sorted_lengths[np.where(cumulative_sum >= half_sum)[0][0]]

    # Create a new summary DataFrame with the desired format
    phenotype = "TERT+" if telomere_cv < 0.7 else "ALT+"

    # Create global summary DataFrame
    global_summary = pd.DataFrame(
        {
            "Read count": alignment_df.shape[0],
            # "global mean / std" means the mean / std of all reads post-filtering
            "Telomere length mean": telomere_mean,
            "Telomere length SD": telomere_std,
            "Telomere length max": telomere_max,
            "Telomere length N50": [telomere_n50],
            "Telomere length CV": [round(telomere_cv, 2)],
            "Phenotype": [phenotype],

            # get mean / std of per-ref means (i.e. the mean across all chromosome arms)
            "Telomere length per-contig mean": round(
                per_ref_summary["Telomere_mean"].mean(), 2
            ),
            "Telomere length per-contig SD": round(
                per_ref_summary["Telomere_mean"].std(), 2
            ),
            # alignment accuracy mean and std
            "Accuracy mean": round(alignment_df["Acc"].mean(), 2),
            "Accuracy SD": round(alignment_df["Acc"].std(), 2),
        },
        index=[0],
    )

    # Save the summary DataFrame to CSV
    global_summary.to_csv(f"{args.output_prefix}.csv", index=False)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("summarise_lengths_results")
    parser.add_argument(
        "--seqkit-bam",
        required=True,
        help="seqkit bam output TSV file",
    )
    parser.add_argument(
        "--ref-fai",
        required=True,
        help="reference FASTA .fai index file",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        required=True,
        help="drop contigs with fewer reads than this",
    )
    parser.add_argument(
        "--telomere-lengths",
        required=True,
        help="tab-separated read IDs and telomere lengths (no header line)",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="prefix for output CSV files",
    )
    return parser
