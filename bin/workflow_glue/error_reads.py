"""Script to process error reads data and print ID results to a file."""

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    error_locs_df = pd.read_csv(
        args.error_locations,
        delimiter="\t",
        dtype={"seqID": str, "start": int},
    ).rename(columns={"start": "error_pos"})

    # Import reference fai as dataframe for list of chr that should be in list
    telomere_locs_df = pd.read_csv(
        args.telomere_locations,
        sep="\t",
        dtype={"seqID": str, "start": int},
    ).rename(columns={"start": "last_telomere_pos"})

    # Merge by seqID
    merged_df = error_locs_df.merge(telomere_locs_df, on="seqID")

    # Remove rows where with errors after the telomere boundary
    merged_df = merged_df[merged_df["error_pos"] <= merged_df["last_telomere_pos"]]

    # Prepare to write output to a file
    with open(args.output, "w") as outfile:
        # Iterate through groups (i.e. through reads) and count occurrences of erroneous
        # k-mers
        for name, group in merged_df.groupby("seqID"):
            err_positions = group["error_pos"].to_numpy()
            for err_pos in err_positions:
                # for each error position we check how many others are within
                # `args.window_size`; if more than `args.max_count`, this read should be
                # removed
                count = (np.abs(err_positions - err_pos) <= args.window_size).sum()
                if count >= args.max_count:
                    outfile.write(f"{name}\n")
                    break


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("ErrReads")
    parser.add_argument(
        "telomere_locations",
        help="`seqkit locate` output file of telomeres",
    )
    parser.add_argument(
        "error_locations",
        help="`seqkit locate` output file of erroneous k-mers",
    )
    parser.add_argument(
        "output",
        help="file to write read IDs to remove",
    )
    parser.add_argument(
        "window_size",
        type=int,
        help="size of sliding window",
    )
    parser.add_argument(
        "max_count",
        type=int,
        help="maximum number of occurrences of erroneous k-mers to filter out a read",
    )
    return parser
