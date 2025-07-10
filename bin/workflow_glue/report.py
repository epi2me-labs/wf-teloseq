"""Create workflow report."""

import json
from pathlib import Path

from dominate import tags
from ezcharts import kdeplot
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.grid import IGridStyles
from ezcharts.layout.snippets.table import DataTable
from ezcharts.layout.util import css
from ezcharts.plots import Plot
import natsort
import pandas as pd
from si_prefix import si_format
from workflow_glue.util import get_named_logger, wf_parser  # noqa: ABS101


def _format_dataframes(
    dataframe,
    base_columns=None,
    thousand_sep_columns=["Read count", "Min length", "Q1", "Median length", "Q3", "Max length"],  # noqa: E501
):
    """Format columns in the dataframe before adding to report.

    `base_columns` are formatted to strings with a human readable SI prefix
    (to nearest thousand) applied.
    `thousand_sep_columns` are formatted to strings with ',' as the thousand separator.

    Both columns parameters are optional. If both are set to None, an unaltered
    dataframe is returned.
    """
    if base_columns is not None:
        dataframe[base_columns] = dataframe[base_columns].apply(
            lambda bases: f"{si_format(bases, precision=2)}B"
        )

    if thousand_sep_columns is not None:
        dataframe[thousand_sep_columns] = dataframe[thousand_sep_columns].apply(
            lambda x: x.apply(lambda y: f"{y:,.0f}")
        )
    return dataframe


def check_for_stats_file(stats_dir, sub_dir, stats_file):
    """Check for the presence of a stats data file in `stats_dir`.

    Returns True if there is at least one TSV file of that type.
    `stats_file` should be one of the following values:
    "kde", "boxplot", "contig", "qc", "unaligned"

    :param stats_dir: The top directory containing all staged stats directories
    :param sub_dir: Sub directory in the staged directory
    :param stats_file: The type of file to check for. For example
      `kde_stats`, `boxplot_stats`.
    """
    possible_types = {"kde", "boxplot", "contig", "qc", "unaligned"}
    if stats_file not in possible_types:
        raise ValueError(f"{stats_file} should be one of {possible_types}")
    # rglob doesn't follow symlinks, but the iterdir path does - so go
    # into the symlinked staged directory and check if the data file exists
    # Return true at the first one we find under the assumption that if one sample
    # has this data all samples must
    return any(
        (
            bool(
                tuple((p / sub_dir).rglob(f"*{stats_file}*.tsv"))
            )
            for p in stats_dir.iterdir()
        )
    )


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "wf-teloseq sequencing report",
        "wf-teloseq",
        args.params,
        args.versions,
        args.workflow_version,
    )

    with open(args.params, encoding="utf-8") as f:
        params = json.load(f)

    # This array is in the same order as the fastcat directories
    with open(args.metadata, encoding="utf-8") as f:
        meta = json.load(f)
        sample_input_order = [m["alias"] for m in meta]
        sample_human_order = natsort.natsorted(sample_input_order)
    root_dir = Path(".")
    stats_dir = args.report_stats_dir
    # Natsorting here is vital - previously this was not naturally sorted
    # which led to samples and fastcat stats directories being mismatched [CW-6342]
    # (Order of staged directories would no longer match order of meta array)
    sample_to_staged_dir = dict(
        zip(
            sample_input_order,
            natsort.natsorted(
                [
                    p for p in root_dir.iterdir()
                    if p.is_dir() and "staged" in p.name
                ]
            )
        )
    )
    with report.add_section("Read summary", "Read summary"):
        tags.p(
            f"""
            Read quality and length distributions, and the percentage yield for each
            read length, for reads above {params["min_length"]} bases and Q score
            {params["read_quality"]}. Use the dropdown menu to view the plots for
            the individual samples.
            """
        )

        # Filter both sample_ids and stats_dirs to ensure they match
        valid_samples = [
            (sample_alias, f"{stats_dir}/fastcat_stats")
            for sample_alias, stats_dir in sample_to_staged_dir.items()
        ]
        if valid_samples:
            sample_ids, stats_dirs = zip(*valid_samples)
            SeqSummary(seq_summary=stats_dirs, sample_names=sample_ids)

    if check_for_stats_file(root_dir, stats_dir, "qc"):
        with report.add_section("Filtering outcomes", "Filtering outcomes"):
            tags.p(
                """
                Reads must pass several filtering checks
                before being analysed. Reads which fail one check are excluded
                from further checks. A failing read is tagged with a valid SAM tag
                in the output BAM, in the format
                """,
                tags.code("qc:Z:<Status>"),
                ". Passing reads are tagged with the status: ",
                tags.code("Good"),
                ". Note that the median read length here is for the whole read, and "
                "not the telomeric repeats. The identity used here is gap-compressed"
                " identity, where a value of 1 represents a perfect alignment.",
            )
            tabs = Tabs()
            for sample in sample_human_order:
                with tabs.add_tab(sample):
                    # Set the column widths, 70% 30%
                    # Table (left) vs. definitions (right)
                    igrid = IGridStyles()
                    igrid.container = css(
                        "display: grid",
                        "grid-template-rows: repeat(1, 1fr)",
                        "grid-template-columns: 66% 34%",
                        "grid-column-gap: 10px",
                        "grid-row-gap: 20px",
                    )
                    with Grid(columns=None, styles=igrid):
                        file_path = sample_to_staged_dir[sample] / stats_dir / f"{sample}_qc_modes_metrics.tsv"  # noqa: E501
                        df = pd.read_csv(file_path, sep="\t")
                        df = _format_dataframes(
                            df, thousand_sep_columns=["Median read length"]
                        )
                        DataTable.from_pandas(
                            df, use_index=False, searchable=False, paging=False
                        )
                        with tags.div(style="padding-top: 2rem;"):
                            definitions = [
                                ('TooShort', 'The read was too short (less than 160 bases).'),  # noqa: E501
                                ('TooFewRepeats', 'There were fewer than 20 telomeric repeat motifs across the entire read.'),  # noqa: E501
                                ('TooCloseStart', 'The telomeric boundary is too close (within 60 bases) to the read start.'),  # noqa: E501
                                ('TooCloseEnd', 'The telomeric boundary is too close (within 60 bases) to the read end.'),  # noqa: E501
                                ('StartNotRepeats', 'The first 30% of the read is < 80% repeats.'),  # noqa: E501
                                ('LowSubTeloQual', 'The mean basecall Q score of the region after the boundary is below 9.'),  # noqa: E501
                                ('TelomereOnly', 'Sequence after boundary is rich in CCC kmers, and is likely extra telomere after a misidentified boundary.'),  # noqa: E501
                                ('TooErrorful', 'A large number of known basecalling error motifs have been observed in the subtelomere.'),  # noqa: E501
                                ('BadAlign', 'A low gap-compressed identity to the reference or low mapQ score. Only applicable to aligned reads.'),  # noqa: E501
                                ('Good', 'Passes all filtering, and is included in final analyses of telomere lengths.')  # noqa: E501
                            ]
                            tags.b("Status definitions (in order of filtering):")
                            status_list = tags.ul()

                            for label, description in definitions:
                                status_list += tags.li(
                                    tags.b(label), " - ", description
                                )

    else:
        tags.p("No QC statistics files were found.")
    # Check we have files to display
    if check_for_stats_file(root_dir, stats_dir, "unaligned"):
        with report.add_section(
            "Alignment-free telomere measurements",
            "Alignment-free telomere measurements"
        ):
            tags.p(
                """Bulk estimate(s) of telomere lengths from reads passing all QC
                   steps.
                   In the event no reads passed for a sample, all values will be 0.
                   Note: The min and max values are within the median +/- 1.5X
                   interquartile range.
                   The coefficient of variation (CV) is a statistical measure that
                   expresses the dispersion of data points around the mean, calculated
                   as the ratio of the standard deviation to the mean.
                   A lower value indicates less variation in measured telomere lengths.
                   """
            )
            dfs = []
            for sample in sample_human_order:
                unaligned_file_path = (
                    sample_to_staged_dir[sample] / stats_dir /
                    f"{sample}_telomere_unaligned_metrics.tsv"
                )
                dfs.append(pd.read_csv(unaligned_file_path, sep="\t"))
            df = pd.concat(dfs)
            if df.empty:
                tags.p("No reads passed filtering, no statistics were generated.")
            else:
                # Human readable yields, thousands sep (TODO - this might
                #  break sorting??)
                df = df.sort_values("Sample", key=natsort.natsort_key)
                df = _format_dataframes(df)
                DataTable.from_pandas(
                    df, use_index=False, searchable=False, paging=False
                )

    if check_for_stats_file(root_dir, stats_dir, "kde"):
        with report.add_section(
            "Single molecule telomere lengths",
            (
                "Kernel density estimate plot of alignment",
                "-free bulk telomere lengths"
            ),
        ):
            tags.p(
                """
                Kernel density plots of all QC passing reads for each sample.
                These plots show the distribution of lengths of all telomeric
                reads in the chosen sample.
                """
            )
            tabs = Tabs()
            for sample in sample_human_order:
                with tabs.add_tab(sample):
                    kde_file_path = sample_to_staged_dir[sample] / stats_dir / f"{sample}_kde_data.tsv"  # noqa: E501
                    df = pd.read_csv(kde_file_path, sep="\t")
                    if df.empty:
                        tags.p(
                            (
                                "KDE requires at least 2 reads to"
                                " pass telomere filtering steps."
                            )
                        )
                    else:
                        plt = kdeplot(
                            (df["density"].values, df["length"].values),
                            compute=False
                        )
                        plt._fig.x_range.range_padding = 0.0
                        EZChart(plt)

    # Check we have files to display
    if check_for_stats_file(root_dir, stats_dir, "contig"):
        with report.add_section(
            "Aligned telomere lengths per contig", "Per contig summary"
        ):
            tags.p(
                """Telomere length metrics grouped by contig.
                Only primary alignments from
                reads which passed the QC filtering stages above are considered.
                Supplementary and secondary alignments for reads are not included.
                Note: The min and max values are within the median +/-
                1.5X interquartile range.
                The coefficient of variation (CV) is a statistical measure that
                expresses the dispersion of data points around the mean, calculated
                as the ratio of the standard deviation to the mean.
                A lower value indicates less variation in measured telomere lengths.
                """
            )
            tabs = Tabs()
            for sample in sample_human_order:
                with tabs.add_tab(sample):
                    try:
                        # At this point there HAS to be only one file for this sample
                        # Fun little unpacking trick
                        contig_file_path = (
                            sample_to_staged_dir[sample] / stats_dir /
                            f"{sample}_contig_telomere_aligned_metrics.tsv"
                        )
                        df = pd.read_csv(contig_file_path, sep="\t")
                        df = _format_dataframes(df)
                        DataTable.from_pandas(df, use_index=False)
                    except pd.errors.EmptyDataError:
                        tags.p("""Sample had no reads which meet the criteria for stats
                                generation. Either all reads failed QC validation
                                (check Filtering outcomes table for absence of "Good"
                                reads), or there were no primary alignments.""")

    if check_for_stats_file(root_dir, stats_dir, "boxplot"):
        with report.add_section(
            "Aligned telomere length boxplots", "Aligned telomere length boxplots"
        ):
            tags.p(
                """Boxplot displaying aggregated lengths for telomeres
                grouped by alignment target contig. Only primary alignments from
                reads which passed the QC filtering stages above are considered.
                Supplementary and secondary alignments for reads are not included."""
            )
            tabs = Tabs()
            for sample in sample_human_order:
                with tabs.add_tab(sample):
                    boxplot_data_file_path = (
                        sample_to_staged_dir[sample] / stats_dir /
                        f"{sample}_boxplot_values.tsv"
                    )
                    df = pd.read_csv(boxplot_data_file_path, sep="\t")
                    if df.empty:
                        tags.p("""Sample had no reads which meet the criteria for stats
                            generation. Either all reads failed QC validation
                            (check Filtering outcomes table for absence of "Good"
                            reads), or there were no primary alignments.""")
                    else:
                        plt = Plot()
                        plt.tooltip = dict(trigger="axis")
                        plt.xAxis = dict(
                            name="Contig",
                            type="category",
                            axisLabel=dict(
                                rotate=90, interval=0, hideOverlap=False, fontSize=8
                            ),
                            axisTick=dict(show=False, interval=0),
                        )
                        plt.yAxis = dict(
                            name="Length",
                            type="value",
                        )
                        plt.dataset = [
                            dict(
                                dimensions=[
                                    "contig",
                                    "min",
                                    "Q1",
                                    "median",
                                    "Q3",
                                    "max",
                                ],
                                source=df.values.tolist(),
                            )
                        ]
                        plt.add_series(
                            dict(
                                type="boxplot",
                                name=sample,
                                encode=dict(
                                    itemName="contig",
                                    tooltip=["min", "Q1", "median", "Q3", "max"],
                                ),
                            )
                        )
                        plt.fix_axis_labels()
                        EZChart(plt)
    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument("--metadata", default="metadata.json", help="sample metadata")
    parser.add_argument(
        "--versions",
        required=True,
        help="directory containing CSVs containing name,version.",
    )
    parser.add_argument(
        "--params",
        default=None,
        required=True,
        help="A JSON file containing the workflow parameter key/values",
    )
    parser.add_argument(
        "--revision", default="unknown", help="git branch/tag of the executed workflow"
    )
    parser.add_argument(
        "--commit", default="unknown", help="git commit of the executed workflow"
    )
    parser.add_argument(
        "--workflow-version", default="unknown",
        help="The version of the workflow run to generate this report.",
    )
    parser.add_argument(
        "--report-stats-dir", default=Path("stats"), type=Path,
        help="Sub directory name containing stats TSV files.",
    )

    return parser
