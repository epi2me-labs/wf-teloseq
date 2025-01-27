"""
Filter and extract data from input files.

Aligned reads are filtered / split into three categories:

* no filter: primary alignments
* lenient: the alignment extends beyond the telomere boundary + an extra margin
  (`--beyond_telomere_margin`)
* strict: the read extends (almost) all the way to the restriction cut site (within
  `--close_to_cutsite_margin`)
"""

import pandas as pd

from .util import wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    # Read input files
    seqkit_bam_df = pd.read_csv(args.seqkit_bam_out, sep="\t")
    cutsites_df = pd.read_csv(
        args.cut_sites_bed,
        sep="\t",
        header=None,
        names=["Ref", "Start", "cut_site"],
        usecols=["Ref", "cut_site"],
    )
    telomere_ends_df = pd.read_csv(
        args.telomere_ends_bed,
        sep="\t",
        header=None,
        names=["Ref", "telomere_end"],
        usecols=["Ref", "telomere_end"],
    )

    # Merge data
    merged_data = pd.merge(seqkit_bam_df, cutsites_df, on="Ref", how="inner")
    merged_data = pd.merge(merged_data, telomere_ends_df, on="Ref", how="inner")

    # Drop secondary and supplementary alignments
    primary_alignments_df = merged_data.query("IsSec == 0 and IsSup == 0")

    # Dynamically determine lenient filter thresholds per `Ref`
    primary_alignments_df["lenient_threshold"] = primary_alignments_df.apply(
        lambda row: min(
            row["cut_site"] - args.close_to_cutsite_margin,  # Cut site - margin
            row["telomere_end"] + args.beyond_telomere_margin  # Telomere + margin
        ),
        axis=1,
    )

    # Apply lenient filter
    beyond_telomere_boundary_ids = primary_alignments_df.query(
        "EndPos >= lenient_threshold"
    )["Read"]

    # Strict filter (unchanged logic)
    close_to_cutsite_ids = primary_alignments_df.query(
        "EndPos >= cut_site - @args.close_to_cutsite_margin"
    )["Read"]

    chr21_p_reads = primary_alignments_df[
        primary_alignments_df['Ref'].str.contains('chr21') &
        (primary_alignments_df['Ref'].str.contains('_p') |
         primary_alignments_df['Ref'].str.contains('_P'))]["Read"]

    # df with chr21_p included again as the cut site is too far away to use.
    combined_df = pd.concat(
        [close_to_cutsite_ids, chr21_p_reads],
        ignore_index=True
    ).drop_duplicates()

    # Write IDs of "no filter" reads
    with open(args.output_no_filter, "w") as outfile:
        for read_id in primary_alignments_df["Read"]:
            outfile.write(read_id + "\n")

    # Write IDs of "lenient filter" reads
    with open(args.output_lenient_filter, "w") as outfile:
        for read_id in set(beyond_telomere_boundary_ids).union(
            set(close_to_cutsite_ids)
        ):
            outfile.write(read_id + "\n")

    # Write IDs of "strict filter" reads
    with open(args.output_strict_filter, "w") as outfile:
        for read_id in set(combined_df):
            outfile.write(read_id + "\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("FiltBam")
    parser.add_argument(
        "seqkit_bam_out",
        help="Seqkit bam TSV file.",
    )
    parser.add_argument(
        "cut_sites_bed",
        help="BED file with cutsite coordinates.",
    )
    parser.add_argument(
        "telomere_ends_bed",
        help="BED file with telomere end coordinates.",
    )
    parser.add_argument(
        "beyond_telomere_margin",
        type=int,
        help=(
            "How far beyond the telomere boundary a read has to extend to pass "
            "the lenient filter."
        ),
    )
    parser.add_argument(
        "close_to_cutsite_margin",
        type=int,
        help="How close to the cut site a read has to reach to pass the strict filter.",
    )
    parser.add_argument(
        "output_no_filter",
        help="Output file for unfiltered mapped read IDs.",
    )
    parser.add_argument(
        "output_lenient_filter",
        help="Output file for IDs of mapped reads that passed the lenient filter.",
    )
    parser.add_argument(
        "output_strict_filter",
        help="Output file for IDs of mapped reads that passed the strict filter.",
    )
    return parser
