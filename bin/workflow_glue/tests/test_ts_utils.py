"""Tests for shared utils."""

import numpy as np
import pandas as pd
import pytest
from workflow_glue.ts_utils import (
    calculate_cv,
    calculate_n50,
    process_telomere_stats,
)


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


@pytest.mark.parametrize(
    "series, expected_read_count, expected_yield, expected_min, expected_mean, expected_max, expected_n50, expected_cv",  # noqa: E501
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
            1500,  # Yield (sum)
            100,  # Min length
            300,  # Mean length
            500,  # Max length
            int(calculate_n50(pd.Series([100, 200, 300, 400, 500]))),  # N50
            round(calculate_cv(np.array([100, 200, 300, 400, 500])), 2),  # CV
        ),  # Normal case
        (
            pd.Series([150]),
            1,  # Read count
            150,  # Yield (sum)
            150,  # Min length
            150,  # Mean length
            150,  # Max length
            150,  # N50
            0.0,  # CV (should be 0 for single value)
        ),
        (
            pd.Series([200, 200, 200, 200, 200]),
            5,  # Read count
            1000,  # Yield (sum)
            200,  # Min length
            200,  # Mean length
            200,  # Max length
            200,  # N50
            0.0,  # CV (should be 0 for identical values)
        ),
    ],
)
def test_process_telomere_stats(
    series,
    expected_read_count,
    expected_yield,
    expected_min,
    expected_mean,
    expected_max,
    expected_n50,
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
        assert result["Yield"].iloc[0] == expected_yield, "Yield incorrect"
        assert result["Min length"].iloc[0] == expected_min, "Min length incorrect"
        assert result["Mean length"].iloc[0] == expected_mean, "Mean length incorrect"
        assert result["Max length"].iloc[0] == expected_max, "Max length incorrect"
        assert result["N50"].iloc[0] == expected_n50, "N50 incorrect"
        assert np.isclose(result["CV"].iloc[0], expected_cv, equal_nan=True), (
            "CV incorrect"
        )
