"""Tests for the functionaility in workflow-glue process_reads."""

from collections import defaultdict, namedtuple
import os
from pathlib import Path
import random

import pandas as pd
import pysam
import pytest
from workflow_glue.process_reads import (
    BoundaryFinder,
    find_telo_boundary,
    largest_error_cluster,
    main,
    process_telomere_stats,
    TELOMERE_MOTIF,
    trim_adapters,
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
    with pysam.FastxFile(
        SCRIPT_DIR / "static/test_telomere_boundary.fastq.gz"
    ) as fastq_file:
        for record in fastq_file:
            records[record.name] = record
    return records


@pytest.mark.parametrize(
    "query_name, expected_boundary, expected_class",
    [
        ("0f88584f-bfd2-4b5c-9a19-5ef2402a54bf", 4262, BoundaryFinder.Good),
        ("b4470370-96c3-41ab-b1e6-11ba396c5524", 3783, BoundaryFinder.Good),
        ("a03facf6-9418-4d92-b0f1-7550c1569f38", 1827, BoundaryFinder.Good),
        ("b8090f53-643c-48b2-bc74-30a2d233e9a7", 2510, BoundaryFinder.Good),
        (
            "62726334-0660-4759-a7d2-b30c966fc098",
            None,
            BoundaryFinder.LowQuality,
        ),
        (
            "502e8d1b-e5dc-4c3d-94f8-54a067e1b411",
            None,
            BoundaryFinder.LowQuality,
        ),
        (
            "b6e75c9f-eb90-4bb9-b809-ed141e474d9a",
            None,
            BoundaryFinder.LowQuality,
        ),
        (
            "669ee9fa-2206-4124-9d0c-d879519f31fc",
            None,
            BoundaryFinder.StartNotRepeats,
        ),
        (
            "7dbba504-6c22-4c53-8abf-6e68fd5bd8db",
            None,
            BoundaryFinder.StartNotRepeats,
        ),
        (
            "4546cfd2-f09d-4594-af29-32a480db025c",
            None,
            BoundaryFinder.StartNotRepeats,
        ),
        (
            "45a5e92c-6ee7-44a9-9b30-affa1581b645",
            None,
            BoundaryFinder.TooFewRepeats,
        ),
        (
            "ddb5ebba-c386-4ce8-91c5-31a1fca25922",
            None,
            BoundaryFinder.TooFewRepeats,
        ),
        (
            "23ecb278-2c5d-469a-90cc-8a8e7eabd747",
            None,
            BoundaryFinder.TooFewRepeats,
        ),
        (
            "too_close_to_end",
            None,
            BoundaryFinder.TooCloseToEnd,
        ),
    ],
)
def test_find_telo_boundary(query_name, expected_boundary, expected_class, records):
    """Test the boundary detection algorithm on a small known dataset."""
    record = records[query_name]
    boundary, classification = find_telo_boundary(record, TELOMERE_MOTIF)
    assert classification == expected_class
    assert boundary == expected_boundary


@pytest.mark.parametrize(
    "query_name",
    [
        ("too_errorful"),
    ],
)
def test_error_clustering(query_name, records):
    """Test the max error detection - read has at least 20 non overlapping errors in."""
    record = records[query_name]
    boundary, _classification = find_telo_boundary(record, TELOMERE_MOTIF)
    clust_size = largest_error_cluster(record.sequence, boundary)
    assert clust_size >= 20


def random_dna(length):
    """Generate any length of random garbage nucleotides."""
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
    record = trim_adapters(record, adapters, prefix, max_errors)
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
        # Empty table
        ([], *[None] * 5),
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
    # Special case empty input
    if not input_data:
        assert summary_df is None
        return
    assert isinstance(summary_df, pd.DataFrame)
    assert summary_df["Read count"].iloc[0] == expected_read_count
    assert summary_df["Telomere length mean"].iloc[0] == expected_mean
    assert summary_df["Telomere length max"].iloc[0] == expected_max
    assert summary_df["Telomere length N50"].iloc[0] == expected_n50
    assert round(summary_df["Telomere length CV"].iloc[0], 2) == round(expected_cv, 2)


def test_main_too_short(capsys, tmp_path):
    """Test main function.

    Input is at least one read which hits every BoundaryFinder variant.
    """
    # Change into the temporary path, so written out files get blasted into
    # oblivion afterwards
    os.chdir(tmp_path)
    Args = namedtuple(
        "Args",
        [
            "input",
            "filter_width",
            "min_repeats",
            "start_window",
            "start_repeats",
            "min_qual_non_telo",
            "barcode",
            "skip_trimming",
            "max_errors",
            "error_distance",
        ],
    )

    args = Args(
        input=str(
            (SCRIPT_DIR / "static" / "test_telomere_boundary.fastq.gz").resolve()
        ),
        filter_width=10,
        min_repeats=20,
        start_window=0.3,
        start_repeats=0.8,
        min_qual_non_telo=9,
        barcode=None,
        skip_trimming=False,
        max_errors=5,
        error_distance=500,
    )

    # Run main function
    main(args)
    captured = capsys.readouterr()
    assert "TooShort: 3, TooFewRepeats: 1, StartNotRepeats: 3, LowQuality: 4, Good: 3, TooErrorful: 2, TooCloseToEnd: 1" in captured.err  # noqa: E501
