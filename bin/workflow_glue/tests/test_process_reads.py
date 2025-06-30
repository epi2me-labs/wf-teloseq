"""Tests for the functionaility in workflow-glue process_reads."""

from collections import Counter, defaultdict
import csv
import os
from pathlib import Path
import random
from unittest.mock import Mock

import numpy as np
import pysam
import pytest
from workflow_glue.process_reads import (
    BoundaryFinder,
    calculate_kde,
    find_telo_boundary,
    largest_error_cluster,
    main,
    TELOMERE_MOTIF,
    trim_adapters,
)

random.seed(0)

SCRIPT_DIR = script_dir = Path(__file__).parent
TEST_BAM = SCRIPT_DIR / "static" / "test_telomere_boundary.bam"
TEST_SHORT_SUBTELO_BAM = SCRIPT_DIR / "static" / "short_subtelomere.bam"


@pytest.fixture
def bam_file():
    """Yield BAM file for tests."""
    with pysam.AlignmentFile(TEST_BAM, "rb", check_sq=False) as bam:
        yield bam


@pytest.fixture
def short_subtelo_bam():
    """Yield BAM file for tests."""
    with pysam.AlignmentFile(TEST_SHORT_SUBTELO_BAM, "rb", check_sq=False) as bam:
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
        ("0f88584f-bfd2-4b5c-9a19-5ef2402a54bf", 4301, BoundaryFinder.Good),
        ("b4470370-96c3-41ab-b1e6-11ba396c5524", 3795, BoundaryFinder.Good),
        ("b8090f53-643c-48b2-bc74-30a2d233e9a7", 2516, BoundaryFinder.Good),
        ("a03facf6-9418-4d92-b0f1-7550c1569f38", 1833, BoundaryFinder.TelomericOnly),
        ("62726334-0660-4759-a7d2-b30c966fc098", 1125, BoundaryFinder.LowQuality),
        ("502e8d1b-e5dc-4c3d-94f8-54a067e1b411", 1541, BoundaryFinder.LowQuality),
        ("b6e75c9f-eb90-4bb9-b809-ed141e474d9a", 2882, BoundaryFinder.LowQuality),
        ("669ee9fa-2206-4124-9d0c-d879519f31fc", 858, BoundaryFinder.StartNotRepeats),
        ("7dbba504-6c22-4c53-8abf-6e68fd5bd8db", 4663, BoundaryFinder.StartNotRepeats),
        ("4546cfd2-f09d-4594-af29-32a480db025c", 4945, BoundaryFinder.Good),
        ("45a5e92c-6ee7-44a9-9b30-affa1581b645", None, BoundaryFinder.TooFewRepeats),
        ("ddb5ebba-c386-4ce8-91c5-31a1fca25922", None, BoundaryFinder.TooFewRepeats),
        ("23ecb278-2c5d-469a-90cc-8a8e7eabd747", None, BoundaryFinder.TooFewRepeats),
        ("too_close_to_end", 599, BoundaryFinder.TooCloseEnd),
    ],
)
def test_find_telo_boundary(query_name, expected_boundary, expected_class, records):
    """Test the boundary detection algorithm on a small known dataset."""
    record = records[query_name]
    boundary, classification = find_telo_boundary(
        record,
        TELOMERE_MOTIF,
        filter_width=10,
        min_repeats=100,
        start_window=0.3,
        start_repeats=0.8,
        min_qual_non_telo=9,
        post_boundary_ccc_threshold=0.25,
    )
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
    boundary, _classification = find_telo_boundary(
        record,
        TELOMERE_MOTIF,
        filter_width=10,
        min_repeats=100,
        start_window=0.3,
        start_repeats=0.8,
        min_qual_non_telo=9,
        post_boundary_ccc_threshold=0.25,
    )
    clust_size = largest_error_cluster(record.query_sequence, boundary)
    assert clust_size >= 20


def random_dna(length):
    """Generate any length of random garbage nucleotides."""
    return "".join(random.choices("ACGT", k=length))


def mutate_barcode(sequence, n_mutations):
    """Mutate the barcode, for testing edlib alignment."""
    mutate_points = (3, 7, 9, 11)
    mutated_barcode = list(sequence)
    for i in range(n_mutations):
        char = mutated_barcode[mutate_points[i]]
        mutated_barcode[mutate_points[i]] = "T" if char != "T" else "A"
    return "".join(mutated_barcode)


# Define test sequences with embedded barcodes
BARCODES = ["AGCTTAGCTTAGCTTAGCTT", "CGTACGTACGTACGTACGTA"]
SEQUENCES = [
    random_dna(100) + barcode + random_dna(200)  # Barcode after 200 bases
    for barcode in BARCODES * 2  # Generate 4 sequences
]
# Barcode just outside prefix
BARCODE_AFTER_PREFIX = f"{random_dna(201)}{BARCODES[0]}{random_dna(100)}"
# Barcode with 2 mismatches within first 200bp, should match
BARCODE_MISMATCH_2 = (
    f"{random_dna(100)}{mutate_barcode(BARCODES[0], 2)}{random_dna(200)}"  # noqa: E501
)
# Barcode with 4 mismatches within first 200bp, shouldn't match
BARCODE_MISMATCH_4 = (
    f"{random_dna(100)}{mutate_barcode(BARCODES[0], 4)}{random_dna(200)}"  # noqa: E501
)


@pytest.mark.parametrize(
    "sequence, adapters, prefix, max_errors, expected_trim, expected_length",
    [
        # Exact barcode matches
        *[(SEQUENCES[i], BARCODES, 200, 3, True, 200) for i in range(4)],
        # Barcode outside prefix — should not match
        (BARCODE_AFTER_PREFIX, BARCODES, 200, 3, False, None),
        # Barcode with 2 mismatches — still within error tolerance
        (BARCODE_MISMATCH_2, BARCODES, 200, 3, True, 200),
        # Barcode with 4 mismatches — should not trim (max errors is 3)
        (BARCODE_MISMATCH_4, BARCODES, 200, 3, False, None),
        # No barcode match and no telomere motif match — should not trim
        (SEQUENCES[0], ["TTTTTTTTTTTTTTTTTT"], 200, 3, False, None),
    ],
)
def test_trim_adapters(
    sequence, adapters, prefix, max_errors, expected_trim, expected_length
):
    """Test that mock reads have the adapter and barcode correctly trimmed."""
    record = Mock()
    record.query_sequence = sequence
    record.query_qualities = "I" * len(sequence)
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
        assert trimmed_qual == "I" * len(sequence), (
            "Quality scores should remain unchanged"
        )


def test_trim_adapters_telomere_fallback():
    """Test trimming when barcode is absent but telomere motif is found."""
    # Telomere motif appears at position 50-56 (HW alignment), within prefix
    # We don't actually trim out the motif
    telomere_motif = "CCTAACC"
    sequence = "A" * 50 + telomere_motif + "G" * 250
    quality_str = "I" * len(sequence)

    record = Mock()
    record.query_sequence = sequence
    record.query_qualities = quality_str

    result = trim_adapters(
        record, adapters=["TTTTTTTTTTTTTTTTTT"], prefix=200, max_errors=3
    )
    expected_trim_position = 50  # start location of motif

    assert result.query_sequence == sequence[expected_trim_position:]
    assert result.query_qualities == quality_str[expected_trim_position:]


@pytest.mark.parametrize(
    "values, expected",
    [
        ([], ["length", "density"]),  # Should fail and create Empty TSV (empty input)  # noqa: E501
        ([420], ["length", "density"]),  # Should fail and create Empty TSV (Only 1 value)  # noqa: E501
        ([1, 1], ["length", "density"]),  # Should fail and create Empty TSV (1 unique value)  # noqa: E501
        (
            [1, 2],
            [
                'length',
                'density',
                1.000000000000000000e+00,
                4.306592395074496649e-01
            ]
        ),  # Should pass (2 unique values)
        (
            [1, 2, 3, 4, 5, 6],
            [
                "length",
                "density",
                1.000000000000000000e+00,
                1.115873735550648727e-01,
                2.000000000000000000e+00,
                1.509601468422716863e-01,
                3.000000000000000000e+00,
                1.641127053753461129e-01,
                4.000000000000000000e+00,
                1.641127053753461407e-01,
                5.000000000000000000e+00,
                1.509601468422717141e-01,
            ],
        ),
    ],
)
def test_calculate_kde(values, expected, tmp_path):
    """Check kde calculated as expected."""
    out_kde = tmp_path / "test_kde.tsv"
    calculate_kde(values, out_kde, step=1)
    result = []
    with open(out_kde) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            result.extend(row)
    # Empty TSV
    if len(result) == 2:
        assert result == expected
    # Has values
    else:
        # Matching headers
        assert result[:2] == expected[:2]
        # Pairup values and assert floats are similar to each other
        for exp, kde in zip(expected[2:], list(map(float, result[2:]))):
            assert np.isclose(exp, kde)


def test_too_close_too_end(short_subtelo_bam):
    """Test the too close too end for reads which have short sub telomere."""
    count_classifications = Counter()
    for record in short_subtelo_bam:
        b, c = find_telo_boundary(
            record,
            TELOMERE_MOTIF,
            filter_width=10,
            min_repeats=100,
            start_window=0.3,
            start_repeats=0.8,
            min_qual_non_telo=9,
            post_boundary_ccc_threshold=0.25,
        )
        count_classifications[c.value] += 1
    # These reads should be Good, were taken from Chr5,
    # where cutsite is ~60 bases from end of telomere boundary
    assert count_classifications["Good"] == 995
    assert count_classifications["TooCloseEnd"] == 23


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

    args.input_bam = str((SCRIPT_DIR / "static" / "test_main.bam").resolve())
    args.output_bam = "/dev/null"
    args.summary_tsv_name = summary_tsv
    args.sample = "TEST"
    args.filter_width = 10
    args.min_repeats = 100
    args.start_window = 0.3
    args.start_repeats = 0.5
    args.min_qual_non_telo = 9
    args.barcode = None
    args.skip_trimming = False
    args.max_errors = 5
    args.error_distance = 500
    args.post_boundary_ccc_threshold = 0.25

    # Run main function
    main(args)
    captured = capsys.readouterr()
    # Added a read that gets trimmed:
    assert (
        "Skipping read d89f9da1-cc54-414e-b014-e1fcb053cd4b, as it has no sequence left"
        in captured.err
    )
    assert (
        "Good: 26, TelomericOnly: 29, TooCloseEnd: 34, TooErrorful: 19, StartNotRepeats: 5, TooShort: 41, TooFewRepeats: 32, LowQuality: 11"  # noqa: E501
        in captured.err
    )
    assert summary_tsv.exists()
