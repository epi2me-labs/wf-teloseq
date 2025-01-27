"""Script to combine statistics CSV."""
import pathlib

import pandas as pd

from .util import wf_parser  # noqa: ABS101

# TODO: this seems like it could be done elsewhere?


def main(args):
    """Run the entry point."""
    output_row_names = {
        args.raw: "Raw Reads",
        args.no_filter: "Mapped: No Filter",
        args.low_filter: "Mapped: Lenient Filter",
        args.high_filter: "Mapped: Strict Filter",
    }

    dfs = []
    for infile, row_name in output_row_names.items():
        df = pd.read_csv(infile)
        df.index = [row_name]
        dfs.append(df)

    comb = pd.concat(dfs)
    comb.to_csv(args.output)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("CombiCsv")
    parser.add_argument(
        "raw",
        type=pathlib.Path,
        help="CSV file with telomere stats of raw reads",
    )
    parser.add_argument(
        "no_filter",
        type=pathlib.Path,
        help="CSV file with telomere stats of unfiltered reads",
    )
    parser.add_argument(
        "low_filter",
        type=pathlib.Path,
        help="CSV file with telomere stats of low-filtered reads",
    )
    parser.add_argument(
        "high_filter",
        type=pathlib.Path,
        help="CSV file with telomere stats of high-filtered reads",
    )
    parser.add_argument(
        "output",
        type=pathlib.Path,
        help="filename for output CSV",
    )
    return parser
