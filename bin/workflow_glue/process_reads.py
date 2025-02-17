"""
Process input reads, returning high quality reads spanning the telomere repeat boundary.

Module applies preflight checks to reads, and detects the end boundary of
telomeric repeats.
Reads passing all checks are emitted to stdout.
"""

from collections import Counter
import csv
import enum
import re
import sys

import edlib
import numpy as np
import pandas as pd
import pysam
from scipy import ndimage

from .util import wf_parser  # noqa: ABS101

# the star of the show, actually a reverse complement and permutation
# of the usual TTAGGG
TELOMERE_MOTIF = "TAACCC"

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


class BoundaryFinder(enum.Enum):
    """Enumeration of possible boundary finder states."""

    Good = 0
    # The read was too short for analysis
    TooShort = 1
    # Too few telomere repeats in the whole read
    TooFewRepeats = 2
    # Unable to determine the boundary for whatever reason
    FailedAnalysis = 3
    # Read start (where repeats should be) did not have enough repeats
    StartNotRepeats = 4
    # The telomere boundary is too close to the end of the read, likely false positive
    TooCloseToEnd = 5
    # Too many error kmers clustered together
    TooErrorful = 6
    # Basecalls after boundary have low Q score
    LowQuality = 7


def find_telo_boundary(
    record,
    motif,
    filter_width=8,
    min_repeats=20,
    start_window=0.3,
    start_repeats=0.8,
    min_qual_non_telo=9,
):
    """
    Return the telomere boundary in the provided read, or None if read fails checks.

    Returns
    -------
    :return: The detected telomere boundary position and classification status.
    :rtype: tuple[int | None, BoundaryFinder]
    """
    # find motif
    motifs = np.zeros(len(record.sequence), dtype=int)
    matches = 0
    for match in re.finditer(motif, record.sequence):
        matches += 1
        motifs[match.start(): match.end()] = 1

    # quick return for few repeats (absolute or relative to filter width)
    if matches < max(min_repeats, filter_width):
        return None, BoundaryFinder.TooFewRepeats

    # common error is to have few repeats at the start of the read
    # require the start of read to be composed of repeats
    # TODO: evaluate this heuristic more
    start = int(len(record.sequence) * start_window)
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
        boundary = np.where(edges == min_value)[0][-1]
    else:
        return None, BoundaryFinder.CannotAnalyse

    # a boundary within the filter width is more likely a false positive
    if len(record.sequence) - boundary < width:
        return None, BoundaryFinder.TooCloseToEnd

    # if theres a low quality region ofter the telomere region,
    # thats likely the basecaller when down the wrong track
    # and we've misidentified the boundary
    quals = record.get_quality_array()[boundary:]
    if np.median(quals) < min_qual_non_telo:
        return None, BoundaryFinder.LowQuality

    return boundary, BoundaryFinder.Good


def largest_error_cluster(sequence, last_position, distance=500):
    """
    Retrieve size of largest error motif clusters in the given sequence.

    Only bases before `last_position` are searched.

    Returns
    -------
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


def trim_adapters(record, adapters, prefix=200, max_errors=3):
    """Trim adapters and barcode from read, updating qualities length as well.

    Doesn't trim if there is no barcode match.
    """
    seq_ = record.sequence[:prefix]
    trim = None
    for _ibc, bc in enumerate(adapters, start=1):
        hits = edlib.align(bc, seq_, mode="HW", task="path")
        if hits["editDistance"] <= max_errors:
            trim = hits["locations"][0][1] + 1
            break  # exceedingly unlikely to have multiple hits
    if trim is not None:
        record.set_sequence(record.sequence[trim:], quality=record.quality[trim:])
    return record


def calculate_n50(series):
    """Calculate N50 for a pandas series using np.searchsorted."""
    if not pd.api.types.is_numeric_dtype(series):
        raise ValueError("Input series is non numeric")
    sorted_lengths = series.sort_values(ascending=False).values
    cumulative_sum = np.cumsum(sorted_lengths)
    half_sum = cumulative_sum[-1] / 2
    index = np.searchsorted(cumulative_sum, half_sum)
    return sorted_lengths[index]


def calculate_cv(series, ddof=1):
    """
    Calculate the coefficient of variance for a pandas series.

    Uses the std deviation with 1 as the degrees of freedom,
    as this is likely to be sample of a larger population.

    Parameters
    ----------
    :param ddof: Degrees of freedom for calculating stddev.

    """
    if not pd.api.types.is_numeric_dtype(series):
        raise ValueError("Input series is non numeric")
    series_mean = np.mean(series)
    series_std = np.std(series, ddof=ddof)
    return np.nan_to_num(series_std / series_mean)


def dump_to_csv(filename, data, header=None):
    """Dump elements in the provided list to a csv.

    Parameters
    ----------
    :param filename: Filename to write to.
    :param data: list of sequences to be written. Each sequence must be equal length
    :param header: Header of the csv. If None, no header is written
    """
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        if header:
            writer.writerow(header)
        writer.writerows(data)


def process_telomere_stats(telomere_lengths):
    """
    Process telomere length statistics and output summary metric and per record CSVs.

    Compute various summary statistics, returning the results and the input as
      dataframes.

    Parameters
    ----------
    :param telomere_lengths: List of tuples containing read IDs and telomere lengths.
    :type telomere_lengths: list[tuple[str, int]]

    Returns
    -------
    :param summary_df: Summary dataframe containing agregated stats for
        the provided read lengths
    :param telomere_length_df: The raw lengths of the data frame.
    """
    telomere_length_df = pd.DataFrame(
        telomere_lengths, columns=["Read_ID", "Telomere_length"]
    )
    if telomere_length_df.empty:
        return None

    # Define aggregation functions
    agg_functions = {
        "Read count": "size",
        "Telomere length mean": "mean",
        "Telomere length SD": "std",
        "Telomere length max": "max",
        "Telomere length N50": calculate_n50,
    }

    # Aggregate statistics
    summary_stats = telomere_length_df["Telomere_length"].agg(agg_functions)
    # Wierdly the CV function throws an error if included in the agg functions,
    # so it's run separately
    # See https://stackoverflow.com/questions/68091853/python-cannot-perform-both-aggregation-and-transformation-operations-simultaneo  # noqa: E501
    summary_stats["Telomere length CV"] = calculate_cv(
        telomere_length_df["Telomere_length"]
    )
    summary_df = summary_stats.to_frame().T
    # Round to human friendly DP
    summary_df.round(
        {
            "Telomere length mean": 0,
            "Telomere length SD": 2,
            "Telomere length max": 0,
            "Telomere length N50": 2,
        }
    )
    # These make the most sense as whole ints
    summary_df[["Telomere length mean", "Telomere length N50"]] = summary_df[
        ["Telomere length mean", "Telomere length N50"]
    ].astype(int)
    return summary_df


def main(args):
    """Process input file to find telomere boundaries."""
    barcodes = BARCODES
    if args.barcode is not None:
        barcodes = [BARCODES[args.barcode - 1]]  # humans count from 1
    # Track a count of how we do on these reads
    boundary_result_count = Counter()
    boundaries = []
    with pysam.FastxFile(args.input) as fastq_in:
        for i, record in enumerate(fastq_in):
            if i % 25000 == 0:
                sys.stderr.write(f"Processed {i} reads\n")

            # maybe some other preflight checks
            if len(record.sequence) < 2 * args.filter_width * len(TELOMERE_MOTIF):
                boundary_result_count[BoundaryFinder.TooShort] += 1
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
            )

            if classification != BoundaryFinder.Good:
                boundary_result_count[classification] += 1
                continue

            # remove reads with pathological errors
            # this test is expensive, do it last
            largest_cluster = largest_error_cluster(
                record.sequence, boundary, distance=args.error_distance
            )
            if largest_cluster > args.max_errors:
                boundary_result_count[BoundaryFinder.TooErrorful] += 1
                continue

            # trim adapters + barcodes
            if not args.skip_trimming:
                record = trim_adapters(
                    record, barcodes
                )
            # we have a good read
            boundary_result_count[BoundaryFinder.Good] += 1
            sys.stdout.write(f"{str(record)}\n")
            boundaries.append((record.name, boundary))
    summary_df = process_telomere_stats(
        boundaries
    )
    summary_df.to_csv("sample_raw_coverage.csv", index=False)
    # Wrire telomeric sequence lengths to CSV
    dump_to_csv(
        "sample_raw_per_read_telomere_length.csv",
        boundaries,
        ("Read ID", "Telomere_length"),
    )
    pretty_str = ", ".join(
        f"{key.name}: {value}" for key, value in boundary_result_count.items()
    )
    sys.stderr.write(f"{pretty_str}\n")


def argparser():
    """Argument parser for entry point."""
    parser = wf_parser("ProReads")
    parser.add_argument("--input", help="Input FASTQ/FASTA file")
    parser.add_argument(
        "--output", default=sys.stdout, type=str, help="Output BAM file"
    )
    # Motif and read filtering options
    grp = parser.add_argument_group(
        "Motif Detection and basic read filtering",
        "Parameters for telomere motif detection.",
    )
    grp.add_argument(
        "--min-repeats", type=int, default=20,
        help="Minimum number of motif repeats to keep read.",
    )
    grp.add_argument(
        "--min-qual-non-telo", type=int, default=10,
        help="Minimum median qscore of non-telomeric region.",
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

    # Repeat filtering options
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
    # Edge detection options
    grp = parser.add_argument_group(
        "Edge Detection", "Parameters for telomere edge detection."
    )
    grp.add_argument(
        "--filter-width", type=int, default=10,
        help="Width of edge filter window, as a multiple of motif length.",
    )
    # Error filtering options
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
