"""Tests for wf-glue `process_alignments.py`.

Hits every function, description of test data in tests/static.
"""

from pathlib import Path
from unittest.mock import Mock

import pysam
import pytest
from workflow_glue.process_alignments import (
    determine_haplotype,
    main,
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
