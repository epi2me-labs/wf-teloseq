"""Tests for wf-glue `process_alignments.py`.

Hits every function, description of test data in tests/static.
"""

from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pysam
import pytest
from workflow_glue.process_alignments import (
    determine_haplotype,
    main,
    mean_quality
)


SCRIPT_DIR = Path(__file__).parent
TEST_BAM = SCRIPT_DIR / "static" / "test_ref_reads.bam"


# If this didn't work I would quit and open a bakery
@pytest.mark.parametrize(
    "ref, expected",
    [
        ("pat1", "Pat"),
        ("mat2", "Mat"),
        ("random", "UH"),
        ("PAT-XYZ", "Pat"),
        ("MAT-XYZ", "Mat"),
    ],
)
def test_assign_haplotype(ref, expected):
    """Test high tech haplotype assignment."""
    assert determine_haplotype(ref) == expected


# If this didn't would I would quit and open a patisserie
def test_mean_quality():
    """Test mean quality calculation."""
    # return tag value if present
    record = Mock()
    record.get_tag = Mock(return_value=5)
    assert mean_quality(record) == 5

    # redo calculation if tag not present
    record = Mock()
    record.get_tag = Mock(side_effect=KeyError('qs'))
    record.query_qualities = [10, 20, 30, 40, 50]
    assert np.isnan(mean_quality(record, trim=10))  # trimmed away
    assert np.isclose(mean_quality(record, trim=0), 16.5321, atol=0.001)  # all
    assert np.isclose(mean_quality(record, trim=3), 42.5963, atol=0.001)  # 40, 50


@pytest.fixture
def bam_file():
    """Yield BAM file for tests."""
    with pysam.AlignmentFile(TEST_BAM, "rb") as bam:
        yield bam


def test_main(tmp_path):
    """Integration test for the main function."""
    output_bam = tmp_path / "output.bam"
    boxplot_stats_tsv = tmp_path / "chr_box_plot_data.tsv"
    summary_stats_tsv = tmp_path / "summary_stats.tsv"
    qc_tsv = tmp_path / "qc_modes.tsv"
    contig_summary_tsv = tmp_path / "contig_summary.tsv"

    args = Mock()
    args.input_bam = TEST_BAM
    args.output_bam = output_bam
    args.sample = "TEST"
    args.summary_tsv_name = summary_stats_tsv
    args.boxplot_tsv_name = boxplot_stats_tsv
    args.qc_tsv_name = qc_tsv
    args.contig_summary_tsv_name = contig_summary_tsv
    args.identity_threshold = 0.8
    args.mapq_threshold = 20

    main(args)

    assert output_bam.exists()
    # Verify tagged reads
    with pysam.AlignmentFile(output_bam, "rb", check_sq=False) as bam:
        for record in bam:
            assert record.has_tag("HP")
            assert record.has_tag("tl")
        expected_record_count = 834
        assert (
            sum(x.mapped for x in bam.get_index_statistics())
            == expected_record_count
        )

    assert boxplot_stats_tsv.exists()
    assert qc_tsv.exists()
    assert summary_stats_tsv.exists()
