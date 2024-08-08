#!/usr/bin/env python
"""Script to analyze telomere lengths and coverage data."""
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
    alignments_df = parse_seqkit(args.seqkit_bam)
    lengths = pd.read_csv(
        args.telomere_lengths,
        sep="\t",
        header=None,
        index_col="Read",
        names=["Read", "Telomere_length"],
    )

    alignments_df["Telomere_length"] = lengths

    # Group by reference and aggregate statistics
    per_ref_summary = alignments_df.groupby("Ref").agg(
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
    alignments_df = alignments_df.loc[alignments_df["Ref"].isin(per_ref_summary.index)]

    # read ref `.fai` index to get order of reference sequences
    refs = list(pd.read_csv(args.ref_fai, sep="\t", usecols=[0], header=None).squeeze())

    # re-order summary to match ref FASTA
    per_ref_summary = per_ref_summary.loc[
        [ref for ref in refs if ref in per_ref_summary.index]
    ]

    # Save filtered summary to CSV
    per_ref_summary.to_csv(f"{args.output_prefix}_chr_arm_Coverage.csv")

    # Save per-read information including telomere length to CSV
    alignments_df.to_csv(f"{args.output_prefix}_Per_Read_telomere_length.csv")

    # Create global summary DataFrame
    global_summary = pd.DataFrame(
        {
            "Read count": alignments_df.shape[0],
            # "global mean / std" means the mean / std of all reads post-filtering
            "Telomere length mean": round(alignments_df["Telomere_length"].mean(), 2),
            "Telomere length SD": round(alignments_df["Telomere_length"].std(), 2),
            "Telomere length max": int(alignments_df["Telomere_length"].max()),
            # get mean / std of per-ref means (i.e. the mean across all chromosome arms)
            "Telomere length per-contig mean": round(
                per_ref_summary["Telomere_mean"].mean(), 2
            ),
            "Telomere length per-contig SD": round(
                per_ref_summary["Telomere_mean"].std(), 2
            ),
            # alignment accuracy mean and std
            "Accuracy mean": round(alignments_df["Acc"].mean(), 2),
            "Accuracy SD": round(alignments_df["Acc"].std(), 2),
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
