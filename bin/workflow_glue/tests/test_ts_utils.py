"""Tests for shared utils."""

import numpy as np
import pandas as pd
import pytest
from workflow_glue.ts_utils import calculate_cv, calculate_n50


@pytest.mark.parametrize(
    "series, expected",
    [
        (pd.Series([1, 3, 3, 7, 10]), 7),
        (pd.Series([1, 2, 3, 4, 5, 6]), 5),
        (pd.Series([100, 200, 300, 400, 500]), 400),
        (np.array([]), None),
        ([1, 2, 3, 4, 5, 6], 5),
    ],
)
def test_calculate_n50(series, expected):
    """Test N50 gives expected results."""
    assert calculate_n50(series) == expected


@pytest.mark.parametrize(
    "series", [pd.Series(["a", 7, 9]), (np.array([np.nan, np.nan, np.nan, np.nan]))]
)
def test_calculate_n50_invalid(series):
    """Should error."""
    with pytest.raises(
        ValueError, match="Input sequence is non numeric or contains NaN values."
    ):
        calculate_n50(series)


@pytest.mark.parametrize(
    "array, expected",
    [
        (np.array([10, 20, 30, 40, 50]), np.float64(0.52704627)),
        (pd.Series([5, 5, 5, 5]), 0),
        (np.array([1, 100, 1000, 10000]), np.float64(1.74305)),
        (np.array([42]), 0),
    ],
)
def test_calculate_cv(array, expected):
    """Test CV gives expected results."""
    assert np.isclose(calculate_cv(array), expected)


@pytest.mark.parametrize(
    "array", [np.array(["x", "y", "z"]), [np.nan, np.nan, np.nan, np.nan], ("hello")]
)
def test_calculate_cv_invalid(array):
    """Should error."""
    with pytest.raises(
        ValueError, match="Input sequence is non numeric or contains NaN values."
    ):
        calculate_cv(array)
