"""Tests for wf-glue `process_alignments.py`.

Hits every function, description of test data in tests/static.
"""

from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pysam
import pytest
from workflow_glue.process_alignments import determine_haplotype, main, mean_qscore


SCRIPT_DIR = Path(__file__).parent
TEST_MAIN_BAM = SCRIPT_DIR / "static" / "test_main_alignment.bam"


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
    assert mean_qscore(record) == 5

    # redo calculation if tag not present
    record = Mock()
    record.get_tag = Mock(side_effect=KeyError("qs"))
    record.query_qualities = [10, 20, 30, 40, 50]
    assert np.isnan(mean_qscore(record, trim=10))  # trimmed away
    assert np.isclose(mean_qscore(record, trim=0), 16.5321, atol=0.001)  # all
    assert np.isclose(mean_qscore(record, trim=3), 42.5963, atol=0.001)  # 40, 50


def test_main(tmp_path):
    """Integration test for the main function."""
    output_bam = tmp_path / "output.bam"
    boxplot_stats_tsv = tmp_path / "chr_box_plot_data.tsv"
    summary_stats_tsv = tmp_path / "summary_stats.tsv"
    qc_tsv = tmp_path / "qc_modes.tsv"
    contig_summary_tsv = tmp_path / "contig_summary.tsv"

    args = Mock()
    args.input_bam = TEST_MAIN_BAM
    args.output_bam = output_bam
    args.sample = "TEST"
    args.summary_tsv_name = summary_stats_tsv
    args.boxplot_tsv_name = boxplot_stats_tsv
    args.qc_tsv_name = qc_tsv
    args.contig_summary_tsv_name = contig_summary_tsv
    args.identity_threshold = 0.8
    args.mapq_threshold = 20
    args.base_stats_dir = "stats"

    main(args)

    assert output_bam.exists()
    # Verify tagged reads
    with pysam.AlignmentFile(output_bam, "rb", check_sq=False) as bam:
        for record in bam:
            assert record.has_tag("HP")
            assert record.has_tag("tl")
        expected_record_count = 992
        assert (
            sum(x.mapped for x in bam.get_index_statistics()) == expected_record_count
        )

    assert boxplot_stats_tsv.exists()
    assert qc_tsv.exists()
    # Check there are no NaNs in the table
    qc_table_text = qc_tsv.read_text()
    assert "nan" not in qc_table_text
    expected_qc = (
        "Status\tTotal reads\tMedian read length\tMedian Q score\tMedian identity\n"
        "TooFewRepeats\t7\t4864\t20.92\t0.99\n"
        "StartNotRepeats\t54\t6958\t15.76\t0.99\n"
        "TooCloseStart\t1\t864\t24.48\t0.00\n"
        "TooCloseEnd\t4\t2236\t20.73\t1.00\n"
        "LowSubTeloQual\t12\t7928\t9.90\t0.90\n"
        "TelomereOnly\t3\t2113\t20.35\t0.99\n"
        "TooErrorful\t69\t8769\t15.49\t0.97\n"
        "BadAlign\t14\t6512\t13.53\t0.95\n"
        "Good\t754\t8641\t20.47\t0.99\n"
    )
    assert qc_table_text == expected_qc
    assert summary_stats_tsv.exists()
