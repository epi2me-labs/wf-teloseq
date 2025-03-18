"""Tests for the functionaility in workflow-glue process_reads."""

from collections import defaultdict
import os
from pathlib import Path
import random
from unittest.mock import Mock

import pysam
import pytest
from workflow_glue.process_reads import (
    BoundaryFinder,
    find_telo_boundary,
    largest_error_cluster,
    main,
    TELOMERE_MOTIF,
    trim_adapters,
)

SCRIPT_DIR = script_dir = Path(__file__).parent
TEST_BAM = SCRIPT_DIR / "static" / "test_telomere_boundary.bam"


@pytest.fixture
def bam_file():
    """Yield BAM file for tests."""
    with pysam.AlignmentFile(TEST_BAM, "rb", check_sq=False) as bam:
        yield bam


@pytest.fixture
def records(bam_file):
    """Open the test bam, yielding a dict of results."""
    records = defaultdict(dict)
    for record in bam_file:
        records[record.query_name] = record
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
    clust_size = largest_error_cluster(record.query_sequence, boundary)
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
    record = Mock()
    record.query_sequence = sequence
    record.query_qualities = quality_str
    record = trim_adapters(record, adapters, prefix, max_errors)
    trimmed_seq, trimmed_qual = record.query_sequence, record.query_qualities
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


def test_main(capsys, tmp_path):
    """Test main function.

    Input is at least one read which hits every BoundaryFinder variant.
    """
    # Change into the temporary path, so written out CSV stats files will
    # get blasted into
    # oblivion afterwards
    os.chdir(tmp_path)
    summary_tsv = tmp_path / "summary.tsv"

    args = Mock()

    args.input_bam = str(
        (SCRIPT_DIR / "static" / "test_telomere_boundary.bam").resolve()
    )
    args.output_bam = "/dev/null"
    args.summary_tsv_name = summary_tsv
    args.sample = "TEST"
    args.filter_width = 10
    args.min_repeats = 20
    args.start_window = 0.3
    args.start_repeats = 0.8
    args.min_qual_non_telo = 9
    args.barcode = None
    args.skip_trimming = False
    args.max_errors = 5
    args.error_distance = 500

    # Run main function
    main(args)
    captured = capsys.readouterr()
    assert "TooShort: 3, TooFewRepeats: 1, LowQuality: 5, StartNotRepeats: 2, Good: 3, TooErrorful: 2, TooCloseToEnd: 1" in captured.err  # noqa: E501
    assert summary_tsv.exists()
