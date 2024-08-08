"""
Filter and extract data from input files.

# TODO: should the terminology of "lenient" and "strict" change?

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
        usecols=[0, 1],
        names=["Ref", "cut_site"],
    )
    telomere_ends_df = pd.read_csv(
        args.telomere_ends_bed,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["Ref", "telomere_end"],
    )

    # Merge data
    merged_data = pd.merge(seqkit_bam_df, cutsites_df, on="Ref", how="inner")
    merged_data = pd.merge(merged_data, telomere_ends_df, on="Ref", how="inner")

    # drop secondary and supplementary alignments
    primary_alignments_df = merged_data.query("IsSec == 0 and IsSup == 0")

    # lenient filter (read extends beyond the telomere boundary + extra margin)
    beyond_telomere_boundary_ids = primary_alignments_df.query(
        "EndPos >= telomere_end + @args.beyond_telomere_margin"
    )["Read"]

    # strict filter (read end pos is at or almost at cut site); note that the
    # restriction site might be fairly far inwards of the telomere boundary, leading to
    # very long fragments (e.g. 100 kb or more for chr21 and EcoRV); for these
    # chromosome arms no reads will pass the filter usually
    close_to_cutsite_ids = primary_alignments_df.query(
        ("EndPos >= cut_site - @args.close_to_cutsite_margin")
    )["Read"]

    # write IDs of "no filter" reads
    with open(args.output_no_filter, "w") as outfile:
        for read_id in primary_alignments_df["Read"]:
            outfile.write(read_id + "\n")

    # write IDs of "lenient filter" reads
    with open(args.output_lenient_filter, "w") as outfile:
        for read_id in set(beyond_telomere_boundary_ids).union(
            set(close_to_cutsite_ids)
        ):
            outfile.write(read_id + "\n")

    # write IDs of "strict filter" reads
    with open(args.output_strict_filter, "w") as outfile:
        for read_id in set(close_to_cutsite_ids):
            outfile.write(read_id + "\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("filter_bam_reads_output")
    parser.add_argument(
        "--seqkit_bam_out",
        required=True,
        help="Seqkit bam TSV file.",
    )
    parser.add_argument(
        "--cut_sites_bed",
        required=True,
        help="BED file with cutsite coordinates.",
    )
    parser.add_argument(
        "--telomere_ends_bed",
        required=True,
        help="BED file with telomere end coordinates.",
    )
    parser.add_argument(
        "--beyond_telomere_margin",
        required=True,
        type=int,
        help=(
            "How far beyond the telomere boundary a read has to extend to pass "
            "the lenient filter."
        ),
    )
    parser.add_argument(
        "--close_to_cutsite_margin",
        required=True,
        type=int,
        help="How close to the cut site a read has to reach to pass the strict filter.",
    )
    parser.add_argument(
        "--output_no_filter",
        required=True,
        help="Output file for unfiltered mapped read IDs.",
    )
    parser.add_argument(
        "--output_lenient_filter",
        required=True,
        help="Output file for IDs of mapped reads that passed the lenient filter.",
    )
    parser.add_argument(
        "--output_strict_filter",
        required=True,
        help="Output file for IDs of mapped reads that passed the strict filter.",
    )
    return parser
