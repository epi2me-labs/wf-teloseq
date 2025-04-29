"""Filter the aligned telomeric sequences, tagging them."""

from pathlib import Path
import sys

import natsort
import numpy as np
import pandas as pd
import pysam
from workflow_glue import ts_utils

from .util import get_named_logger, wf_parser  # noqa: ABS101,E501

NOT_GOOD_PRIMARY = (
        pysam.FUNMAP | pysam.FSECONDARY | pysam.FQCFAIL | pysam.FSUPPLEMENTARY)

QUAL_TO_PROB = np.fromiter(
    (10 ** (-q / 10) for q in range(0, 100)),
    dtype=np.float64, count=100)


def determine_haplotype(ref):
    """Determine the haplotype of a given reference. High tech."""
    ref_type = "UH"
    if ref is not None:
        ref_lower = ref.lower()
        if "pat" in ref_lower:
            ref_type = "Pat"
        elif "mat" in ref_lower:
            ref_type = "Mat"
    return ref_type


def get_tag_with_default(record, tag, default=None):
    """Get a tag from an `AlignedSegment`, returning default if not present."""
    try:
        return record.get_tag(tag)
    except KeyError:
        return default


def mean_quality(record, trim=60):
    """Calculate the mean quality of a read.

    :param trim: number of quality scores to trim from the start of the read,
        the default is to be equivalent to dorado.
    """
    try:
        qual = record.get_tag("qs")
    except KeyError:
        if record.query_qualities is None:
            return np.nan
    else:
        return qual

    probs = np.fromiter(
        (QUAL_TO_PROB[q] for q in record.query_qualities),
        dtype=np.float64, count=len(record.query_qualities))
    try:
        probs = probs[trim:]
    except IndexError:
        return np.nan
    return -10 * np.log10(np.mean(probs))


def main(args):
    """Entry point."""
    logger = get_named_logger(__name__)
    records = []

    with (
        pysam.AlignmentFile(args.input_bam, "r") as bam_in,
        pysam.AlignmentFile(args.output_bam, "wb", template=bam_in) as bam_out,
    ):
        for i, record in enumerate(bam_in):
            if i and not i % 100:
                logger.info(f"Analysed {i} records.")

            record.set_tag(
                "HP", determine_haplotype(record.reference_name), value_type="Z"
            )
            # Only want to analyse good primary mappings
            analyse = False
            # de tag is the gap compressed identity as assigned by minimap2
            # see https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity  # noqa: E501
            # Only check identity for Primary alignments
            if not (record.flag & NOT_GOOD_PRIMARY):
                identity = 1 - record.get_tag("de")
                if (
                    identity < args.identity_threshold or
                    record.mapping_quality < args.mapq_threshold
                ):
                    record.set_tag("qc", "BadAlign", value_type="Z")
                else:
                    # Consider this record in dataframe aggregations
                    analyse = True

            bam_out.write(record)
            # Add line to stats generation
            records.append(
                (
                    record.query_name,
                    record.query_length,
                    record.get_tag("tl"),
                    record.reference_name,
                    record.get_tag("HP"),
                    1 - get_tag_with_default(record, "de", 1),
                    record.get_tag("qc"),
                    mean_quality(record),
                    analyse
                )
            )
    # Make all reads into a dataframe.
    df = pd.DataFrame(
        records,
        columns=(
            "read_id",
            "query_length",
            "telomere_length",
            "reference_name",
            "haplotype",
            "identity",
            "qc",
            "median_qual",
            "analyse"
        ),
    )

    # Compute box plot stats per contig on primary alignments which passed
    # all pipeline filtering
    df_good_reads = df[(df["analyse"]) & (df["qc"] == "Good")]
    boxplot_data = (
        df_good_reads.groupby("reference_name")["telomere_length"]
        .agg(ts_utils.boxplot_stats)
        .apply(pd.Series)
    )
    boxplot_data = boxplot_data.loc[natsort.natsorted(boxplot_data.index)]
    boxplot_data.to_csv(args.boxplot_tsv_name, sep="\t")
    # Compute summary stats for all data
    summary_data = ts_utils.process_telomere_stats(df_good_reads["telomere_length"])
    if summary_data is not None:
        summary_data["Sample"] = args.sample
    else:
        summary_data = pd.DataFrame([])
    summary_data.to_csv(
        args.summary_tsv_name, index=False, sep="\t",
        float_format="%.2f"
    )

    # Compute summary stats for per contig data.
    per_contig_summary_df = df_good_reads.groupby(
        "reference_name", as_index=False
    )["telomere_length"].apply(ts_utils.process_telomere_stats)
    # Handle case where there are no stats
    if per_contig_summary_df.empty:
        pd.DataFrame([]).to_csv(args.contig_summary_tsv_name, index=False, sep="\t")
    else:
        per_contig_summary_df.rename(
            columns={"reference_name": "Contig name"}, inplace=True
        )
        per_contig_summary_df.sort_values(
            "Contig name", inplace=True, key=natsort.natsort_key
        )
        per_contig_summary_df.to_csv(
            args.contig_summary_tsv_name, index=False, sep="\t"
        )

    # QC modes for reads - what about if there is no data
    qc_df = df.set_index("read_id")
    qc_df = (
        qc_df.loc[qc_df.index.drop_duplicates()]
        .groupby("qc", as_index=False)
        .agg(
            {
                "query_length": ["size", "median"],
                "median_qual": "median",
                "identity": "median",
            }
        )
    )
    # Round median
    qc_df[("query_length", "median")] = qc_df[
        ("query_length", "median")
    ].round(0).astype(int)
    # Order the types of errors into the order the filters are applied
    application_order = [
        "TooShort",
        "TooFewRepeats",
        "StartNotRepeats",
        "TooCloseStart"
        "TooCloseEnd",
        "LowSubTeloQual",
        "TelomereOnly",
        "TooErrorful",
        "BadAlign",
        "Good"
    ]
    qc_df["qc"] = qc_df["qc"].astype(
        pd.CategoricalDtype(categories=application_order, ordered=True)
    )
    # Order the df
    qc_df = qc_df.sort_values(by="qc")
    qc_df.to_csv(
        args.qc_tsv_name,
        index=False,
        header=[
            "Status", "Total reads", "Median read length",
            "Median quality", "Median identity"
        ],
        float_format="%.2f",
        sep="\t",
    )
    # Index new BAM
    pysam.index("--csi", str(args.output_bam))


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("tag_bam_reads")

    parser.add_argument(
        "sample", type=str, help="Sample name.",
    )
    parser.add_argument(
        "input_bam", help="Input BAM file. Use - for stdin."
    )
    parser.add_argument(
        "--output-bam", default=sys.stdout, help="Output BAM filename.",
    )
    parser.add_argument(
        "--summary-tsv-name", type=Path, default="summary_stats.tsv",
        help="Name of telomere length summary stats TSV file to output.",
    )
    parser.add_argument(
        "--boxplot-tsv-name", type=Path, default="chr_box_plot_data.tsv",
        help="Name of boxplot stats TSV file to output.",
    )
    parser.add_argument(
        "--qc-tsv-name", type=Path, default="qc_mode_stats.tsv",
        help="Name of assigned qc stats TSV file to output.",
    )
    parser.add_argument(
        "--contig-summary-tsv-name", type=Path, default="contig_summary.tsv",
        help="Name of telomere stats per contig TSV file to output.",
    )
    parser.add_argument(
        "--identity-threshold", type=float, default=0.8,
        help="Minimum gap compressed identity for an alignment to be considered.",
    )
    parser.add_argument(
        "--mapq-threshold", type=int, default=20,
        help="Minimum map quality for an alignment to be considered.",
    )
    return parser
