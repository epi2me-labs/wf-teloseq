"""Shared code for teloseq wf-glue."""

import numpy as np
import pandas as pd


def calculate_cv(data, ddof=1):
    """Calculate the coefficient of variance for a sequence.

    Uses the std deviation with 1 as the degrees of freedom,
    as this is likely to be sample of a larger population.

    :param ddof: Degrees of freedom for calculating stddev.

    :raises ValueError: Raised if input data is non numeric or contains NaN values.
    """
    data = np.array(data)
    if not np.issubdtype(data.dtype, np.number) or np.any(np.isnan(data)):
        raise ValueError("Input sequence is non numeric or contains NaN values.")
    series_mean = np.mean(data)
    series_std = np.std(data, ddof=ddof)
    return np.nan_to_num(series_std / series_mean)


def process_telomere_stats(series):
    """Process telomere length statistics, output summary metrics and CSVs.

    Compute various summary statistics, returning the results and the input as
      dataframes.

    :param telomere_lengths: List of tuples containing read IDs and telomere lengths.
    :type telomere_lengths: list[tuple[str, int]]

    :return: Summary dataframe containing aggregated stats for
        the provided read lengths.
    """
    if series.empty:
        return None

    # Define aggregation functions
    agg_functions = {
        "Read count": "size",
        "Min length": "min",
        "Mean length": "mean",
        "Max length": "max",
    }

    # Aggregate statistics
    summary_stats = series.agg(agg_functions)
    # Annoyingly the CV and N50 functions throw errors if included in the agg functions,
    # so run separately
    # See https://stackoverflow.com/questions/68091853/python-cannot-perform-both-aggregation-and-transformation-operations-simultaneo  # noqa: E501
    summary_stats["CV"] = calculate_cv(series.values)
    summary_df = summary_stats.to_frame().T
    # Round to human friendly DP
    summary_df = summary_df.round({"CV": 2})
    # These make the most sense as whole ints
    int_fields = [
        "Read count",
        "Min length",
        "Mean length",
        "Max length",
    ]
    summary_df[int_fields] = (
        summary_df[int_fields]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .astype(int)
    )

    return summary_df
