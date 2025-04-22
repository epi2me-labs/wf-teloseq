"""Create workflow report."""

import json
from pathlib import Path

from dominate import tags
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
            (sample_alias, fastcat_dir)
            for sample_alias, fastcat_dir in zip(
                sample_input_order, sorted(args.fastcat_stats_dir.iterdir())
            )
        ]
        if valid_samples:
            sample_ids, stats_dirs = zip(*valid_samples)
        SeqSummary(stats_dirs, sample_names=sample_ids)

    if list(args.qc_stats_dir.rglob("*.tsv")):
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
                "not the telomeric repeats. The identity used here is Gap compressed"
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
                        (file_path,) = tuple(
                            args.qc_stats_dir.rglob(f"{sample}*.tsv")
                        )
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
                                ('StartNotRepeats', 'The first 30% of the read is < 80% repeats.'),  # noqa: E501
                                ('TooCloseEnd', 'The telomeric boundary is too close (within 80 bases) to the read end.'),  # noqa: E501
                                ('LowSubTeloQual', 'The mean basecall Q score of the region after the boundary is below 9.'),  # noqa: E501
                                ('TelomereOnly', 'Sequence after boundary is rich in CCC kmers, and is likely extra telomere after a misidentified boundary.'),  # noqa: E501
                                ('TooErrorful', 'A large number of known basecalling error motifs have been observed in the subtelomere.'),  # noqa: E501
                                ('BadAlign', 'A low gap-compressed identity to the reference. Only applicable to aligned reads.'),  # noqa: E501
                                ('Good', 'Passes all filtering, and is included in final analyses of telomere lengths.')  # noqa: E501
                            ]
                            tags.b("Status definitions (in order of filtering):")
                            status_list = tags.ul()

                            for label, description in definitions:
                                status_list += tags.li(
                                    tags.b(label), " - ", description
                                )  # noqa: E501

    else:
        tags.p("No QC statistics file was found.")
    # Check we have files to display
    if list(args.unaligned_stats_dir.rglob("*.tsv")):
        with report.add_section(
            "Alignment-free telomere measurements",
            "Alignment-free telomere measurements"
        ):
            tags.p(
                """Bulk estimate(s) of telomere lengths from reads passing all QC
                   steps.
                   In the event no reads passed for a sample, all values will be 0."""
            )
            dfs = []
            for sample in sample_human_order:
                (file_path,) = tuple(
                    args.unaligned_stats_dir.rglob(f"{sample}*.tsv")
                )
                dfs.append(pd.read_csv(file_path, sep="\t"))
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

    # Check we have files to display
    if list(args.contig_stats_dir.rglob("*.tsv")):
        with report.add_section(
            "Aligned telomere lengths per contig", "Per contig summary"
        ):
            tags.p("""Telomere length metrics grouped by alignment targets.
                Only primary alignments from
                reads which passed the QC filtering stages above are considered.
                Supplementary and secondary alignments for reads are not included.""")

            tabs = Tabs()
            for sample in sample_human_order:
                with tabs.add_tab(sample):
                    try:
                        # At this point there HAS to be only one file for this sample
                        # Fun little unpacking trick
                        (file_path,) = tuple(
                            args.contig_stats_dir.rglob(f"{sample}*.tsv")
                        )
                        df = pd.read_csv(file_path, sep="\t")
                        df = _format_dataframes(df)
                        DataTable.from_pandas(df, use_index=False)
                    except pd.errors.EmptyDataError:
                        tags.p("""Sample had no reads which meet the criteria for stats
                                generation. Either all reads failed QC validation
                                (check Filtering outcomes table for absence of "Good"
                                reads), or there were no primary alignments.""")

    if list(args.boxplot_stats_dir.rglob("*.tsv")):
        with report.add_section(
            "Aligned telomere length boxplots", "Aligned telomere length boxplots"
        ):
            tags.p(
                """Boxplot displaying aggregated lengths for telomeres
                grouped by alignment target contig. Only primary alignments from
                reads which passed the QC filtering stages above are considered.
                supplementary and secondary alignments for reads are not included."""
            )
            tabs = Tabs()
            for sample in sample_human_order:
                with tabs.add_tab(sample):
                    (file_path,) = tuple(
                        args.boxplot_stats_dir.rglob(f"{sample}*.tsv")
                    )
                    df = pd.read_csv(file_path, sep="\t")
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
        "--unaligned-stats-dir", type=Path, default=None,
        help="Directory containing the unaligned stats TSVS for this workflow.",
    )
    parser.add_argument(
        "--fastcat-stats-dir", type=Path, default=None,
        help="Directory containing the fastcat derived statistics for this workflow."
    )
    parser.add_argument(
        "--qc-stats-dir", type=Path, default=None,
        help="Directory containing the qc TSV files.",
    )
    parser.add_argument(
        "--contig-stats-dir", type=Path, default=None,
        help="Directory containing the qc TSV files.",
    )
    parser.add_argument(
        "--boxplot-stats-dir", type=Path, default=None,
        help="Directory containing the boxplot aggregated data TSV files.",
    )

    return parser
