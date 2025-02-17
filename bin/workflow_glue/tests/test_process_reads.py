"""Tests for the functionaility in workflow-glue process_reads."""

from collections import defaultdict
from pathlib import Path
import random

import pandas as pd
import pysam
import pytest
from workflow_glue.process_reads import (
    BoundaryFinder,
    find_telo_boundary,
    process_telomere_stats,
    TELOMERE_MOTIF,
    trim_adapters
)

SCRIPT_DIR = script_dir = Path(__file__).parent


# Function to create a mock FastxRecord
def create_fastx_record(seq, qual, name="test_read"):
    """Create a mock Pysam FastxRecord."""
    record = pysam.FastxRecord()
    record.set_name(name)
    record.set_sequence(seq, quality=qual)
    return record


@pytest.fixture
def records():
    """Open the test bam, yielding a dict of results."""
    records = defaultdict(dict)
    fields = ("query_sequence", "query_qualities")
    with pysam.AlignmentFile(
        SCRIPT_DIR / "static/test.bam", check_sq=False, threads=3
    ) as bam:
        for record in bam:
            records[record.query_name] = {att: getattr(record, att) for att in fields}
    return records


@pytest.mark.parametrize(
    "query_name, expected_boundary, expected_class",
    [
        ("0a0d6af3-1e97-4126-8ffb-7591c0e94cdf", 2678, BoundaryFinder.Good),
        ("0a0e2380-bc8e-454c-90b2-d53370c2c924", 2160, BoundaryFinder.Good),
        ("0a0faaca-299d-4925-83af-a6454c02e56f", 2632, BoundaryFinder.Good),
        ("0a0fb7db-c209-42a7-8d45-d3d04235f711", 2475, BoundaryFinder.Good),
        ("0a0fb381-f2ef-4bbf-a55c-c747019ceb4f", 3848, BoundaryFinder.Good),
        (
            "00a44e02-7303-4a87-89ac-2d349481afd0",
            None,
            BoundaryFinder.LowQuality,
        ),
        (
            "0ae826fd-4d3b-42e2-ab44-6385e3c121cd",
            None,
            BoundaryFinder.LowQuality,
        ),
        (
            "0b1a806a-dd36-473d-a1e0-289f6eb68713",
            None,
            BoundaryFinder.LowQuality,
        ),
        (
            "0a000a74-82f8-406f-ad3e-cf5cb1441efa",
            None,
            BoundaryFinder.StartNotRepeats,
        ),
        (
            "0a0c81d7-1624-4c29-a50c-326b0709ef95",
            None,
            BoundaryFinder.StartNotRepeats,
        ),
        (
            "0a0f89bd-c7d5-4a0b-8029-a36480ab9bee",
            None,
            BoundaryFinder.StartNotRepeats,
        ),
    ],
)
def test_find_telo_boundary(
    query_name, expected_boundary, expected_class, records
):
    """Test the boundary detection algorithm on a small known dataset."""
    sequence = records[query_name]["query_sequence"]
    # fastx record expects a ASCII set of qualities
    qualities = "".join(chr(q + 33) for q in records[query_name]["query_qualities"])
    record = create_fastx_record(seq=sequence, qual=qualities)

    boundary, classification = find_telo_boundary(record, TELOMERE_MOTIF)
    assert classification == expected_class
    assert boundary == expected_boundary


def random_dna(length):
    """Generate any lenght of random garbage nmucleotides."""
    return "".join(random.choices("ACGT", k=length))


# Define test sequences with embedded barcodes
BARCODES = ["AGCTTAGCTTAGCTTAGCTT", "CGTACGTACGTACGTACGTA"]
SEQUENCES = [
    random_dna(100) + barcode + random_dna(200)  # Barcode after 200 bases
    for barcode in BARCODES * 2  # Generate 4 sequences
]
QUALITY_STRS = ["I" * len(seq) for seq in SEQUENCES]  #


@pytest.mark.parametrize(
    "sequence, quality_str, adapters, prefix, max_errors, expected_trim, expected_length",  # noqa: E501
    [
        (SEQUENCES[i], QUALITY_STRS[i], BARCODES, 200, 3, True, 200)
        for i in range(len(SEQUENCES))  # Expect this to be 200 bases long
    ]
    + [
        (SEQUENCES[i], QUALITY_STRS[i], ["TGCATGCATGCATGCATGCA"], 200, 3, False, None)
        for i in range(len(SEQUENCES))  # Wrong barcode, should be full length
    ],
)
def test_trim_adapters(
    sequence, quality_str, adapters, prefix, max_errors, expected_trim, expected_length
):
    """Test that mock reads have the adapter and barcode correctly trimmed."""
    record = create_fastx_record(seq=sequence, qual=quality_str)
    record = trim_adapters(
        record, adapters, prefix, max_errors
    )
    trimmed_seq, trimmed_qual = record.sequence, record.quality
    if expected_trim:
        assert len(trimmed_seq) == expected_length, (
            f"Unexpected length of trimmed seq {len(trimmed_seq)}"
        )
        assert len(trimmed_qual) == expected_length, (
            f"Unexpected length of trimmed qual {len(trimmed_qual)}"
        )

    else:
        assert trimmed_seq == sequence, "Sequence should remain unchanged"
        assert trimmed_qual == quality_str, "Quality scores should remain unchanged"


@pytest.mark.parametrize(
    "input_data, expected_read_count, expected_mean, expected_max, expected_n50, expected_cv",  # noqa: E501
    [
        # Normal case: Various telomere lengths
        (
            [
                ("read1", 1000),
                ("read2", 1500),
                ("read3", 1200),
                ("read4", 1700),
                ("read5", 1300),
            ],
            5,
            1340,
            1700,
            1300,
            0.2,
        ),
        # Identical lengths: Ensures CV is 0
        (
            [("read1", 1000), ("read2", 1000), ("read3", 1000)],
            3,
            1000,
            1000,
            1000,
            0.0,
        ),
        # Zero mean case: Ensures division by zero is handled
        (
            [("read1", 0), ("read2", 0), ("read3", 0)],
            3,
            0,
            0,
            0,
            0,
        ),
        # Single read case: Ensures correct computation when only one read is available
        (
            [("read1", 1500)],
            1,
            1500,
            1500,
            1500,  # N50 should be the only read length
            0,  # CV should be 0 as there’s no variation
        ),
    ],
)
def test_process_telomere_stats(
    input_data,
    expected_read_count,
    expected_mean,
    expected_max,
    expected_n50,
    expected_cv,
):
    """Test calculating read length statistics."""
    summary_df = process_telomere_stats(input_data)

    assert isinstance(summary_df, pd.DataFrame)
    assert summary_df["Read count"].iloc[0] == expected_read_count
    assert summary_df["Telomere length mean"].iloc[0] == expected_mean
    assert summary_df["Telomere length max"].iloc[0] == expected_max
    assert summary_df["Telomere length N50"].iloc[0] == expected_n50
    assert round(summary_df["Telomere length CV"].iloc[0], 2) == round(
        expected_cv, 2
    )
