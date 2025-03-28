"""Tests for shared utils."""

import numpy as np
import pandas as pd
import pytest
from workflow_glue.ts_utils import (
    calculate_cv,
    process_telomere_stats,
)


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


@pytest.mark.parametrize(
    "series, expected_read_count, expected_min, expected_q1, expected_median, expected_q3, expected_max, expected_cv",  # noqa: E501
    [
        (
            pd.Series([], dtype=float),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ),  # Empty series
        (
            pd.Series([100, 200, 300, 400, 500]),
            5,  # Read count
            100,  # Min length
            200,  # Q1
            300,  # Median length
            400,  # Q3
            500,  # Max length
            round(calculate_cv(np.array([100, 200, 300, 400, 500])), 2),  # CV
        ),  # Normal case
        (
            pd.Series([150]),
            1,  # Read count
            150,  # Min length
            150,  # Q1
            150,  # Median length
            150,  # Q3
            150,  # Max length
            0.0,  # CV (should be 0 for single value)
        ),
        (
            pd.Series([200, 200, 200, 200, 200]),
            5,  # Read count
            200,  # Min length
            200,  # Q1
            200,  # Median length
            200,  # Q3
            200,  # Max length
            0.0,  # CV (should be 0 for identical values)
        ),
    ],
)
def test_process_telomere_stats(
    series,
    expected_read_count,
    expected_min,
    expected_q1,
    expected_median,
    expected_q3,
    expected_max,
    expected_cv,
):
    """Parameterized test for process_telomere_stats function."""
    result = process_telomere_stats(series)

    if series.empty:
        assert result is None
    else:
        assert isinstance(result, pd.DataFrame)
        assert result["Read count"].iloc[0] == expected_read_count, (
            "Read count incorrect"
        )
        assert result["Min length"].iloc[0] == expected_min, "Min length incorrect"
        assert result["Q1"].iloc[0] == expected_q1, "Q1 incorrect"
        assert result["Median length"].iloc[0] == expected_median, (
            "Mean length incorrect"
        )
        assert result["Q3"].iloc[0] == expected_q3, "Q3 incorrect"
        assert result["Max length"].iloc[0] == expected_max, "Max length incorrect"
        assert np.isclose(result["CV"].iloc[0], expected_cv, equal_nan=True), (
            "CV incorrect"
        )
