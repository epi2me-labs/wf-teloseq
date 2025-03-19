"""Create workflow report."""

from collections import defaultdict
import json
import os
from pathlib import Path

from dominate import tags
from dominate.util import raw
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.grid import IGridStyles
from ezcharts.layout.snippets.table import DataTable
from ezcharts.layout.util import css
from ezcharts.plots import Plot
import pandas as pd
from si_prefix import si_format
from workflow_glue.util import get_named_logger, wf_parser  # noqa: ABS101


def _format_dataframes(dataframe, base_columns, thousand_sep_columns=None):
    """Format columns in the dataframe before adding to report.

    `base_columns` are formatted to strings with a human readable SI prefix
    (to nearest thousand) applied.
    `thousand_sep_columns` are formatted to strings with ',' as the thousand separator.

    Both columns parameters are optional. If both are set to None, an unaltered
    dataframe is returned.
    """
    if thousand_sep_columns is None:
        thousand_sep_columns = [
            'Read count', 'Min length', 'Mean length', 'Max length', 'N50']

    if base_columns:
        dataframe[base_columns] = dataframe[base_columns].apply(
            lambda bases: f"{si_format(bases, precision=2)}B"
        )

    if thousand_sep_columns:
        dataframe[thousand_sep_columns] = dataframe[thousand_sep_columns].apply(
            lambda x: x.apply(lambda y: f"{y:,.0f}")
        )
    return dataframe


def gather_sample_files(sample_details):
    """
    Create a dictionary of sample paths to data files in work directory.

    If alignment was skipped, alignment related files will be set to None.
    """
    expected_files = {
        "fastcat_stats": "data/{s}/fastcat_stats",
        "telomere_unaligned_stats": "data/{s}/{s}_telomere_unaligned_metrics.tsv",
        "telomere_aligned_stats": "data/{s}/{s}_telomere_aligned_metrics.tsv",
        "box_plot_stats": "data/{s}/{s}_boxplot_values.tsv",
        "qc_modes": "data/{s}/{s}_qc_modes_metrics.tsv",
        "contig_summary": "data/{s}/{s}_contig_telomere_aligned_metrics.tsv",
    }
    sample_files = defaultdict(dict)
    for sample in sample_details:
        for key, filename in expected_files.items():
            # If there is no {} returns unadulterated string
            file_path = Path(filename.format_map({"s": sample["alias"]}))
            if file_path.exists():
                sample_files[sample["alias"]][key] = file_path
            else:
                sample_files[sample["alias"]][key] = None
    return sample_files


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

    with open(args.metadata) as metadata:
        sample_details = sorted(
            [
                {"alias": d["alias"], "type": d["type"], "barcode": d["barcode"]}
                for d in json.load(metadata)
            ],
            key=lambda d: d["alias"],
        )
        # Choose sample files based on mapping report flag
        sample_files = gather_sample_files(sample_details)
    # this file only exists if alignment was performed
    alignment_performed = any(
        sample.get("telomere_aligned_stats") for sample in sample_files.values()
    )
    with open(args.params, encoding="utf-8") as f:
        params = json.load(f)

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
            (sample_id, sample_files[sample_id]["fastcat_stats"])
            for sample_id in sample_files
            if sample_files[sample_id]["fastcat_stats"]
            and os.path.exists(sample_files[sample_id]["fastcat_stats"])
        ]

        if valid_samples:
            sample_ids, stats_dirs = zip(*valid_samples)

        SeqSummary(stats_dirs, sample_names=sample_ids)

    with report.add_section(
        "Per sample telomere stats", "Per sample telomere stats"
    ):
        tags.p("""Statistics on all reads in a sample that have passed all filtering
                and have been identified as spanning beyond the telomere repeat
                boundary. Length derived metrics are calculated using the estimated
                length of the telomere repeats.
                Some key metrics are described in more detail below.""")
        tags.ul(
            tags.li("Yield - The total number of bases sequenced for passing reads in a sample."),  # noqa: E501
            tags.li("CV - Coefficient of Variation, the standard deviation divided by the mean. Useful for comparing variation in measured telomere lengths between samples."),  # noqa: E501
            tags.li("N50 - 50% of total sequenced bases are contained in reads of this length or longer.")  # noqa: E501
        )
        dfs = []
        for _sample, files in sample_files.items():
            if file_path := files.get("telomere_unaligned_stats"):
                dfs.append(pd.read_csv(file_path, sep="\t"))
        if dfs:
            df = pd.concat(dfs)
            # Human readable yields, thousands sep (TODO - this might break sorting??)
            df = _format_dataframes(df, base_columns="Yield")
            DataTable.from_pandas(df, use_index=False, searchable=False, paging=False)
        else:
            tags.p("No statistics were generated.")

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
            ". Note that the median read length here is for the whole read, and not "
            "the telomeric repeats. The identity used here is Gap compressed"
            " identity, where a value of 1 represents a perfect alignment."
        )
        tabs = Tabs()
        for sample, files in sample_files.items():
            with tabs.add_tab(sample):
                if file_path := files.get("qc_modes"):
                    # Set the column widths, 70% 30%
                    # Table (left) vs. definitions (right)
                    igrid = IGridStyles()
                    igrid.container = css(
                        "display: grid",
                        "grid-template-rows: repeat(1, 1fr)",
                        "grid-template-columns: 66% 34%",
                        "grid-column-gap: 10px",
                        "grid-row-gap: 20px")
                    with Grid(columns=None, styles=igrid):
                        df = pd.read_csv(file_path, sep="\t")
                        df = _format_dataframes(
                            df,
                            base_columns=None,
                            thousand_sep_columns=["Median read length"]
                        )
                        DataTable.from_pandas(
                            df, use_index=False, searchable=False, paging=False
                        )
                        with tags.div(style="padding-top: 2rem;"):
                            raw("""
                                <b>Status definitions (in order of filtering):</b>
                                <ul>
                                    <li><b>TooShort</b> - The read was too short
                                 (less than 160 bases).</li>
                                    <li><b>TooFewRepeats</b> - There were fewer
                                 than 20 telomeric repeat motifs across the entire
                                 read.</li>
                                    <li><b>StartNotRepeats</b> - The first 30% of
                                 the read is &lt; 80% repeats.</li>
                                    <li><b>TooCloseEnd</b> - The telomeric boundary
                                 is too close (within 80 bases) to the read end.</li>
                                    <li><b>LowSubTeloQual</b> - The mean basecall
                                 Q score of the region after the boundary
                                 is below 9.</li>
                                    <li><b>TooErrorful</b> - A large number of known
                                 basecalling error motifs have been observed in
                                 the subtelomere.</li>
                                    <li><b>BadAlign</b> - A low gap-compressed
                                 identity to the reference. Only applicable to
                                 aligned reads.</li>
                                    <li><b>Good</b> - Passes all filtering, and is
                                 included in final analyses of telomere lengths.</li>
                                </ul>
                                """)
                else:
                    tags.p("No QC statistics file was found.")

    if alignment_performed:
        with report.add_section(
            "Aligned telomere lengths per contig", "Per contig summary"
        ):
            tags.p("""Telomere length metrics grouped by alignment targets, as found
                in the provided reference.
                Only primary alignments from
                reads which passed the QC filtering stages above are considered.
                Supplementary and secondary alignments for reads are not included.""")

            tags.ul(
                tags.li("Yield - The total number of bases sequenced for passing reads in a sample."),  # noqa: E501
                tags.li("CV - Coefficient of Variation, the standard deviation divided by the mean. Useful for comparing variation in measured telomere lengths between aligned contigs."),  # noqa: E501
                tags.li("N50 - 50% of total sequenced bases are contained in reads of this length or longer.")  # noqa: E501
            )
            tabs = Tabs()
            for sample, files in sample_files.items():
                with tabs.add_tab(sample):
                    if contig_file_path := files.get("contig_summary"):
                        df = pd.read_csv(contig_file_path, sep="\t")
                        df = _format_dataframes(df, base_columns="Yield")
                        DataTable.from_pandas(df, use_index=False)
                    else:
                        tags.p("No alignment metrics file was found.")

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
            for sample, files in sample_files.items():
                with tabs.add_tab(sample):
                    if box_data_path := files.get("box_plot_stats"):
                        df = pd.read_csv(box_data_path, sep="\t")

                        plt = Plot()
                        plt.tooltip = dict(trigger="axis")
                        # plt.grid = dict(containLabel=True, left="left")
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
        "--workflow-version",
        default="unknown",
        help="The verison of the workflow run to generate this report.",
    )
    return parser
