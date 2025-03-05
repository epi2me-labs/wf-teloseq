"""Shared code for teloseq wf-glue."""

import numpy as np


def calculate_n50(data):
    """
    Calculate N50 for a sequence.

    Raises
    ------
    :raises ValueError: Raised if input data is non numeric or contains NaN values.
    """
    data = np.array(data)
    if not np.issubdtype(data.dtype, np.number) or np.any(np.isnan(data)):
        raise ValueError("Input sequence is non numeric or contains NaN values.")
    # Empty input
    if not data.size:
        return None
    sorted_lengths = np.sort(data)
    cumulative_sum = np.cumsum(sorted_lengths)
    half_sum = cumulative_sum[-1] / 2
    index = np.searchsorted(cumulative_sum, half_sum)
    return sorted_lengths[index]


def calculate_cv(data, ddof=1):
    """
    Calculate the coefficient of variance for a sequence.

    Uses the std deviation with 1 as the degrees of freedom,
    as this is likely to be sample of a larger population.

    Parameters
    ----------
    :param ddof: Degrees of freedom for calculating stddev.

    Raises
    ------
    :raises ValueError: Raised if input data is non numeric or contains NaN values.
    """
    data = np.array(data)
    if not np.issubdtype(data.dtype, np.number) or np.any(np.isnan(data)):
        raise ValueError("Input sequence is non numeric or contains NaN values.")
    series_mean = np.mean(data)
    series_std = np.std(data, ddof=ddof)
    return np.nan_to_num(series_std / series_mean)
