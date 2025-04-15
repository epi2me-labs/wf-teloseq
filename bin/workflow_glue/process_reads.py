"""
Process input reads, returning high quality reads spanning the telomere repeat boundary.

Module applies preflight checks to reads, and detects the end boundary of
telomeric repeats.
Reads passing all checks are emitted to stdout.
"""

from collections import Counter
import enum
from pathlib import Path
import re
import sys

import edlib
import numpy as np
import pandas as pd
import pysam
from scipy import ndimage
from workflow_glue import ts_utils

from .util import wf_parser  # noqa: ABS101

# the star of the show, actually a reverse complement and permutation
# of the usual TTAGGG
TELOMERE_MOTIF = "TAACCC"

# Known telomeric variants to catch
VARIANT_MOTIFS = {
    "CACCCT",
    "ACCCCT",
    "CCCAAA",
    "CCCCGA",
    "CCCTGA",
    "CCCTCA",
    "CCCTAC",
    "CCCTAT",
}


# known pathological error motifs
ERROR_MOTIFS = {
    "GTATAG",
    "CGCGCGCG",
    "CCACCG",
    "AGCGACAG",
    "ATAAGT",
    "CCTCGTCC",
    "TATAGT",
    "AGTACT",
    "GAGTCC",
    "TATACA",
    "TGGTCC",
    "CTCTCCTCT",
}
#  Regex pattern for error detection
ERROR_MOTIF_REGEX = re.compile("(?=(" + "|".join(ERROR_MOTIFS) + "))")

# we would hope these have been removed by the basecaller, but we
# allow the option to do it here.
# the structure of reads is: adapter-barcode-telomere
# note that there is an extra C on the end of these compared to their
# usual representation in e.g. dorado source code
# there is a meaningful order here, so not a set
BARCODES = [
    "CACAAAGACACCGACAACTTTCTTC",
    "AAGGTTACACAAACCCTGGACAAGC",
    "AAGGATTCATTCCCACGGTAACACC",
    "GAGAGGACAAAGGTTTCAACGCTTC",
    "TCCGATTCTGCTTCTTTCTACCTGC",
    "AGAACGACTTCCATACTCGTGTGAC",
    "CGTCAACTGACAGTGGTTCGTACTC",
    "CCAAACCCAACAACCTAGATAGGCC",
    "CCAGTAGAAGTCCGACAACGTCATC",
    "GGAGTTCGTCCAGAGAAGTACACGC",
    "CTTTCGTTGTTGACTCGACGGTAGC",
    "CATCTGGAACGTGGTACACCTGTAC",
]
# this sequence is jammed on the end of the adapter-barcode, its complementary
# to telomere motif TTAGGG (plus another extra C) and in the frame found to be
# the most common
BARCODES = [bc + "CCTAACC" for bc in BARCODES]

# I am likely to try and break this, so check that the tag is valid SAM format
# NB. lower case is not reserved by SAM
SAM_TAG_REGEX = re.compile(r"^[a-zA-Z]{2}:[AifZHB]:.+$")

# Compile regex to detect any known telomeric repeats in a read
TELOMERE_MOTIF_REGEX = re.compile(
    r"(" + "|".join(VARIANT_MOTIFS.union({TELOMERE_MOTIF})) + ")"
)

# CCC only motif
CCC_REGEX = re.compile(
    r"CCC"
)


class BoundaryFinder(enum.Enum):
    """Enumeration of possible boundary finder states."""

    Good = "Good"
    # The read was too short for analysis
    TooShort = "TooShort"
    # Too few telomere repeats in the whole read
    TooFewRepeats = "TooFewRepeats"
    # Unable to determine the boundary for whatever reason
    FailedAnalysis = "FailedAnalysis"
    # Read start (where repeats should be) did not have enough repeats
    StartNotRepeats = "StartNotRepeats"
    # The telomere boundary is too close to the end of the read, likely false positive
    TooCloseToEnd = "TooCloseEnd"
    # Too many error kmers clustered together
    TooErrorful = "TooErrorful"
    # Basecalls after boundary have low Q score
    LowQuality = "LowSubTeloQual"
    # Sequence after boundary is rich in CCC, so is likely just extra telomere
    TelomericOnly = "TelomereOnly"


def find_telo_boundary(
    record,
    motif,
    filter_width=8,
    min_repeats=20,
    start_window=0.3,
    start_repeats=0.8,
    min_qual_non_telo=9,
    post_boundary_ccc_threshold=0.25
):
    """Identify telomere boundary in the provided read, or None if read fails checks.

    :return: The detected telomere boundary position and classification status.
    :rtype: tuple[int | None, BoundaryFinder]
    """
    # find motif
    motifs = np.zeros(len(record.query_sequence), dtype=int)
    matches = 0
    for match in TELOMERE_MOTIF_REGEX.finditer(record.query_sequence):
        matches += 1
        motifs[match.start(): match.end()] = 1

    # quick return for few repeats (absolute or relative to filter width)
    if matches < max(min_repeats, filter_width):
        return None, BoundaryFinder.TooFewRepeats

    # common error is to have few repeats at the start of the read
    # require the start of read to be composed of repeats
    # TODO: evaluate this heuristic more
    start = int(len(record.query_sequence) * start_window)
    if np.sum(motifs[:start]) < start_repeats * start:
        return None, BoundaryFinder.StartNotRepeats

    # make an edge filter
    width = len(motif) * filter_width
    edge_filter = np.ones(width + 1, dtype=int)
    edge_filter[width // 2] = 0
    edge_filter[width // 2 + 1:] = -1

    # find edges, the median filter is a sharpening filter
    # that remove artefacts from inexact matches not detected
    # by the motif detection
    motifs = ndimage.median_filter(motifs, size=width)
    edges = np.convolve(motifs, edge_filter, mode="valid")

    # boundary is defined as last drop of large magnitude
    # "large" here is currently a complete half-window
    boundary = None
    min_value = np.min(edges)
    if min_value < filter_width * len(motif):
        boundary = np.where(edges == min_value)[0][-1] + width
    else:
        return None, BoundaryFinder.CannotAnalyse

    # a boundary within the filter width is more likely a false positive
    if len(record.query_sequence) - boundary < width:
        return None, BoundaryFinder.TooCloseToEnd

    # if there's a low quality region after the telomere region,
    # it's likely the basecaller went down the wrong track
    # and we've misidentified the boundary
    if quals := record.query_qualities:  # as we process references too
        if np.median(quals[boundary:]) < min_qual_non_telo:
            return None, BoundaryFinder.LowQuality

    # Check for CCC-rich motifs in post-boundary region
    post_seq = record.query_sequence[boundary:]
    # 3 base motif, so sum 3
    ccc_count = sum(3 for _ in CCC_REGEX.finditer(post_seq))
    ccc_proportion = ccc_count / len(post_seq)
    # If the sequence after the boundary is composed of more than
    # `post_boundary_CCC_threshold` CCC motifs, we deem it likely to be
    # just more telomere, so we tag it to be filtered out
    if ccc_proportion > post_boundary_ccc_threshold:
        return None, BoundaryFinder.TelomericOnly

    return boundary, BoundaryFinder.Good


def largest_error_cluster(sequence, last_position, distance=500):
    """Retrieve size of largest error motif clusters in the given sequence.

    Only bases before `last_position` are searched.

    :return: The size of the largest detected error cluster.
    """
    errors = np.fromiter(
        ((m.start() + 1) for m in ERROR_MOTIF_REGEX.finditer(sequence[:last_position])),
        dtype=int,
    )

    if len(errors) == 0:
        return 0

    # calculate distance matrix, and filter to "close" errors
    # find the largest neighbour count
    diff_matrix = np.abs(errors[:, None] - errors)
    counts = (diff_matrix <= distance).sum(axis=1)
    return np.max(counts)


def trim_adapters(
    record, adapters, prefix=200,
    max_errors=3, fallback_errors=1, adapter_motif="CCTAACC"
):
    """Trim adapters and barcode from read, updating qualities length as well.

    Doesn't trim if there is no barcode match.
    """
    seq_ = record.query_sequence[:prefix]
    trim = None
    for _ibc, bc in enumerate(adapters, start=1):
        hits = edlib.align(bc, seq_, mode="HW", task="path")
        if hits["editDistance"] <= max_errors:
            trim = hits["locations"][0][1] + 1
            break  # exceedingly unlikely to have multiple hits
    # if we didn't find a barcode, we will trim up to the
    # first telomere motif, allowing up to `fallback_error` mismatches
    else:
        hits = edlib.align(adapter_motif, seq_, mode="HW", task="path")
        # Short sequence, default of 1 mismatch
        if hits["editDistance"] <= fallback_errors:
            trim = hits["locations"][0][0]
    if trim is not None:
        trimmed_quals = record.query_qualities[trim:]
        record.query_sequence = record.query_sequence[trim:]
        record.query_qualities = trimmed_quals
    return record


def main(args):
    """Process input file to find telomere boundaries."""
    barcodes = BARCODES
    if args.barcode is not None:
        barcodes = [BARCODES[args.barcode - 1]]  # humans count from 1
    # Track a count of how we do on these reads
    boundary_result_count = Counter()
    boundaries = []
    qualities = []

    with (
        pysam.AlignmentFile(args.input_bam, mode="r", check_sq=False) as bam_in,
        pysam.AlignmentFile(args.output_bam, mode="w", template=bam_in) as bam_out
    ):
        for i, record in enumerate(bam_in):
            if i and i % 25000 == 0:
                sys.stderr.write(f"Processed {i} reads\n")

            # As we are writing everything out now - we need to trim
            # adapters for all reads
            if not args.skip_trimming:
                record = trim_adapters(record, barcodes)
                # In the event the basecaller isn't used to trim (why?)
                # There may be some reads are adapter+barcode, which are
                # trimmed in their entirety
                # This is a special case, so we will just skip the read
                # entirely
                if record.query_sequence is None:
                    sys.stderr.write(
                        f"Skipping read {record.query_name}, as it has no sequence left\n"  # noqa:E501
                    )
                    continue

            # Immediately tag the record with -1 for boundary coordinate
            # Only updated later if a boundary is detected
            record.set_tag("tl", -1, value_type="i")
            # maybe some other preflight checks
            if len(record.query_sequence) < 2 * args.filter_width * len(TELOMERE_MOTIF):
                boundary_result_count[BoundaryFinder.TooShort] += 1
                record.set_tag("qc", BoundaryFinder.TooShort.value, value_type="Z")
                bam_out.write(record)
                continue

            # find telomere boundary
            boundary, classification = find_telo_boundary(
                record,
                TELOMERE_MOTIF,
                filter_width=args.filter_width,
                min_repeats=args.min_repeats,
                start_window=args.start_window,
                start_repeats=args.start_repeats,
                min_qual_non_telo=args.min_qual_non_telo,
                post_boundary_ccc_threshold=args.post_boundary_ccc_threshold
            )

            if boundary is not None:
                record.set_tag("tl", boundary, value_type="i")

            if classification != BoundaryFinder.Good:
                boundary_result_count[classification] += 1
                record.set_tag("qc", classification.value, value_type="Z")
                bam_out.write(record)
                continue

            # remove reads with pathological errors
            # this test is expensive, do it last
            largest_cluster = largest_error_cluster(
                record.query_sequence, boundary, distance=args.error_distance
            )
            if largest_cluster > args.max_errors:
                boundary_result_count[BoundaryFinder.TooErrorful] += 1
                record.set_tag("qc", BoundaryFinder.TooErrorful.value, value_type="Z")
                bam_out.write(record)
                continue

            # we have a good read
            boundary_result_count[classification] += 1
            record.set_tag("qc", classification.value, value_type="Z")
            # Add the telomere boundary location to the comment
            bam_out.write(record)

            # Gather boundary and mean quality of good reads
            boundaries.append(boundary)
            qualities.append(np.mean(record.query_qualities))

    # Write the summary metrics out for use in report.
    summary_data = ts_utils.process_telomere_stats(pd.Series(boundaries))
    if summary_data is not None:
        summary_data.insert(0, "Sample", args.sample)
        summary_data.to_csv(
            args.summary_tsv_name,
            index=False,
            sep="\t",
            float_format="%.2f",
        )
    else:
        # No reads passed filtering, so to include in final report table
        # We write out a TSV which indicates that.
        pd.DataFrame([(args.sample, 0, 0, 0, 0, 0, 0, 0)]).to_csv(
            args.summary_tsv_name, sep="\t", index=False,
            header=[
                "Sample",
                "Read count",
                "Min length",
                "Q1",
                "Median length",
                "Q3",
                "Max length",
                "CV",
            ],
        )

    # Emit a short summary to stderr
    pretty_str = ", ".join(
        f"{key.name}: {value}" for key, value in boundary_result_count.items()
    )
    sys.stderr.write(f"{pretty_str}\n")


def argparser():
    """Argument parser for entry point."""
    parser = wf_parser("ProReads")
    parser.add_argument("sample", help="Sample name.")
    parser.add_argument("input_bam", help="Input BAM file. Use - for stdin.")
    parser.add_argument("--output-bam", default=sys.stdout, help="Output BAM file.")
    parser.add_argument(
        "--summary-tsv-name", type=Path, default="unaligned_summary_stats.tsv",
        help="Name of stats TSV file to output.",
    )

    # Motif and read filtering
    grp = parser.add_argument_group(
        "Motif Detection and basic read filtering",
        "Parameters for telomere motif detection.",
    )
    grp.add_argument(
        "--min-repeats", type=int, default=20,
        help="Minimum number of motif repeats to keep read.",
    )
    grp.add_argument(
        "--min-qual-non-telo", type=int, default=9,
        help="Minimum median qscore of non-telomeric region.",
    )
    grp.add_argument(
        "--post_boundary_ccc_threshold", type=float, default=0.25,
        help="Maximum threshold for proportion of post telomere boundary read"
        " that is composed of CCC motifs.",
    )
    grp.add_argument(
        "--barcode", type=int, default=None,
        help="If provided, trim this numbered barcode (and preceeding adapter) "
        "from the read. "
        "If not provided, read will be searched for all possible barcodes. "
        "Note that this program does not perform demultiplexing.",
    )
    grp.add_argument(
        "--skip-trimming", action="store_true", default=False,
        help="Skip trimming adapters and barcodes off of reads.",
    )

    # Repeat filtering
    grp = parser.add_argument_group(
        "Start Repeat Filtering", "Filtering reads without repeats at the start."
    )
    grp.add_argument(
        "--start-window", type=float, default=0.3,
        help="Fraction of read to consider for start repeats.",
    )
    grp.add_argument(
        "--start-repeats", type=float, default=0.8,
        help="Fraction of start window to require to be repeats.",
    )

    # Edge detection
    grp = parser.add_argument_group(
        "Edge Detection", "Parameters for telomere edge detection."
    )
    grp.add_argument(
        "--filter-width", type=int, default=10,
        help="Width of edge filter window, as a multiple of motif length.",
    )

    # Error filtering
    grp = parser.add_argument_group(
        "Error Filtering", "Filtering reads with pathological errors."
    )
    grp.add_argument(
        "--max-errors", type=int, default=5,
        help="Maximum number of co-located errors allowed in a read.",
    )
    grp.add_argument(
        "--error-distance", type=int, default=500,
        help="Window size to search for co-localized errors.",
    )
    return parser
