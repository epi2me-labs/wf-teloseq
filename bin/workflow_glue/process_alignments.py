"""Filter the aligned telomeric sequences, tagging them."""

from functools import partial
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101,E501


def assign_haplotype(ref):
    """Assign haplotype to given reference. High tech."""
    ref_lower = ref.lower()
    if "pat" in ref_lower:
        ref_type = "Pat"
    elif "mat" in ref_lower:
        ref_type = "Mat"
    else:
        ref_type = "UH"
    return ref_type


def calculate_chr_box_stats(df):
    """Calculate read statistics for a dataframe of read stats.

    Returns a transformed dataframe, containing the min, Q1, median,
    Q3 and max for each contig.
    """
    # Box plot requires the format
    # [Min, Q1, Median, Q3, Max]
    # Define partials for quantiles
    quantile_25 = partial(pd.Series.quantile, q=0.25)
    quantile_75 = partial(pd.Series.quantile, q=0.75)
    numeric_cols = ["query_length", "telomere_length", "identity"]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric)
    chr_box_plots_data = (
        df[(df["tp"] == "P") & (df["qc"] == "Good")]
        .groupby("reference_name")["telomere_length"]
        .agg((np.min, quantile_25, np.median, quantile_75, np.max))
    )
    return chr_box_plots_data


def main(args):
    """Entry point."""
    logger = get_named_logger(__name__)
    # (contig, assigned haplotype) -> raw data about record
    records = []

    with (
        pysam.AlignmentFile(args.input_bam, "r") as bam_in,
        pysam.AlignmentFile(args.output_bam, "wb", template=bam_in) as bam_out,
    ):
        for i, record in enumerate(bam_in):
            if i and not i % 100:
                logger.info(f"Analysed {i} records.")

            record.set_tag(
                "HP",
                assign_haplotype(record.reference_name),
                value_type="Z"
            )
            # de tag is the gap compressed identity as assigned by minimap2
            # see https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity  # noqa: E501
            identity = 1 - record.get_tag("de")
            if identity < args.identity_threshold:
                record.set_tag("qc", "BadAlign", value_type="Z")

            bam_out.write(record)
            # Add line to stats generation
            records.append(
                (
                    record.query_length,
                    record.get_tag("tl"),
                    record.reference_name,
                    record.get_tag("HP"),
                    record.get_tag("de"),
                    record.get_tag("qc"),
                    record.get_tag("tp"),
                )
            )
    # Make all reads into a dataframe.
    df = pd.DataFrame(
        records,
        columns=(
            "query_length",
            "telomere_length",
            "reference_name",
            "haplotype",
            "identity",
            "qc",
            "tp",
        ),
    )

    chr_box_plot_data = calculate_chr_box_stats(df)
    chr_box_plot_data.to_csv(args.boxplot_csv_name)

    pysam.index("--csi", str(args.output_bam))


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("tag_bam_reads")

    parser.add_argument(
            "--input-bam", type=Path, help="Input BAM file.", default=sys.stdin
    )
    parser.add_argument(
        "--output-bam", default=sys.stdout, type=Path,
        help="Output BAM filename.",
    )
    parser.add_argument(
        "--boxplot-csv-name", type=Path, default="chr_box_plot_data.csv",
        help="Name of stats CSV file to output.",
    )
    parser.add_argument(
        "--identity-threshold", type=float, default=0.8,
        help="Minimum gap compressed identity for an alignment to be considered.",
    )
    return parser
