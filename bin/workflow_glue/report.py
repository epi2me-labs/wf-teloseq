"""Create workflow report."""

from collections import defaultdict
import json
import os
from pathlib import Path
import re

from bokeh.io import output_file, save
from bokeh.models import ColumnDataSource, Range1d
from dominate import tags
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from workflow_glue.util import get_named_logger, wf_parser


def gather_sample_files(sample_details):
    """Create a dictionary of sample paths to data files in work directory.

    If mapping was skipped, mapping related files will be set to None.
    NB. These will be symlinks, just FYI
    """
    expected_files = {
        "fastcat_stats": "data/{s}/fastcat_stats",
        "rawcov": "data/{s}/telomere_length_metrics.csv",
        "rawlen": "data/{s}/read_telomere_lengths.csv",
        "telsubstats": "data/{s}/subtelomere_reads_stats.txt",
        "subtellen": "data/{s}/sample_raw_per_read_telomere_length.csv",
        "endmotif": "data/{s}/summary_telomere_motif.txt",
        "none": "data/{s}/{s}_results/{s}_NS_per_read_telomere_length.csv",
        "lenient": "data/{s}/{s}_results/{s}_LS_per_read_telomere_length.csv",
        "strict": "data/{s}/{s}_results/{s}_HS_per_read_telomere_length.csv",
        "strictcov": "data/{s}/{s}_results/{s}_HS_chr_arm_coverage.csv",
        "lenientcov": "data/{s}/{s}_results/{s}_LS_chr_arm_coverage.csv",
        "output": "data/{s}/{s}_results/{s}_combined_summary.csv",
        "rawstats": "raw.txt",
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


def load_csv_with_default(file_path, **kwargs):
    """Attempt to load a CSV file, or return an empty DataFrame if missing."""
    if file_path and os.path.exists(file_path):
        return pd.read_csv(file_path, **kwargs)
    else:
        # Return an empty DataFrame with the expected columns
        return pd.DataFrame(columns=["No data available"])


def create_boxplot(df, column_name, output_filename, chart_name):
    """Create a boxplot for the given column."""
    series = df[column_name]

    plt = BokehPlot(tools="save", x_range=[column_name], title=chart_name)
    p = plt._fig
    p.y_range = Range1d(start=series.min() - 10, end=series.max() + 10)

    kde = gaussian_kde(series)
    kde_vals = kde.evaluate(np.linspace(series.min(), series.max(), 100))

    kde_scale = max(kde_vals) / 0.4
    kde_x_pos = 0.5

    kde_vals_left = -kde_vals / kde_scale + kde_x_pos
    kde_vals_right = kde_vals[::-1] / kde_scale + kde_x_pos
    kde_vals = np.hstack([kde_vals_left, kde_vals_right])
    kde_support = np.hstack(
        [
            np.linspace(series.min(), series.max(), 100),
            np.linspace(series.max(), series.min(), 100),
        ]
    )

    source = ColumnDataSource(data={"x": kde_vals, "y": kde_support})
    p.patch("x", "y", source=source, alpha=0.3)

    padding_top = 10
    p.y_range = Range1d(
        start=series.min() - padding_top, end=series.max() + padding_top
    )

    q1, q3 = series.quantile([0.25, 0.75])
    iqr = q3 - q1
    qmin, q1, q2, q3, qmax = series.quantile([0, 0.25, 0.5, 0.75, 1])
    upper = min(qmax, q3 + 1.5 * iqr)
    lower = max(qmin, q1 - 1.5 * iqr)

    hbar_height = (qmax - qmin) / 500
    whisker_width = 0.1

    p.rect([column_name], lower, whisker_width, hbar_height, line_color="grey")
    p.rect([column_name], upper, whisker_width, hbar_height, line_color="grey")
    p.segment([column_name], upper, [column_name], q3, line_color="grey")
    p.segment([column_name], lower, [column_name], q1, line_color="grey")
    p.vbar([column_name], 0.2, q2, q3, line_color="black")
    p.vbar([column_name], 0.2, q1, q2, line_color="black")

    p.xaxis.major_label_orientation = "vertical"
    p.xaxis.ticker = []  # Removes the ticks
    p.xaxis.major_label_text_font_size = "0pt"  # Hide major tick labels
    p.xaxis.axis_label = column_name
    p.yaxis.axis_label = "Values"

    output_file(output_filename)
    save(p)
    EZChart(plt)


def create_chrboxplot(df, group_column, value_column, output_filename, chart_name):
    """Create a boxplot for chromosome data."""
    grouped = df.groupby(group_column)

    def sort_key(item):
        # Match chromosome like "chrN" where N is a number (prioritize chromosomes)
        chr_match = re.search(r"chr(\d+)", item, re.IGNORECASE)
        if chr_match:
            return (0, int(chr_match.group(1)))  # Chromosomes first, sorted numerically

        # Match contigs like "contig_N" (prioritize contigs after chromosomes)
        contig_match = re.search(r"contig_(\d+)", item, re.IGNORECASE)
        if contig_match:
            return (1, int(contig_match.group(1)))  # Contigs come after chromosomes

        # If neither chromosome nor contig, sex chromosomes at end
        return (2, item)

    # Sort the group keys accordingly
    sorted_group_keys = sorted(grouped.groups.keys(), key=sort_key)

    plt = BokehPlot(tools="save", x_range=sorted_group_keys, title=chart_name)
    p = plt._fig

    for i, group_key in enumerate(sorted_group_keys):
        group_data = grouped.get_group(group_key)
        series = group_data[value_column]

        if len(series) < 2:
            continue

        kde = gaussian_kde(series)
        kde_vals = kde.evaluate(np.linspace(series.min(), series.max(), 100))

        kde_scale = max(kde_vals) / 0.4
        offset = 0.5
        kde_x_pos = i + offset

        kde_vals_left = -kde_vals / kde_scale + kde_x_pos
        kde_vals_right = kde_vals[::-1] / kde_scale + kde_x_pos
        kde_vals = np.hstack([kde_vals_left, kde_vals_right])
        kde_support = np.hstack(
            [
                np.linspace(series.min(), series.max(), 100),
                np.linspace(series.max(), series.min(), 100),
            ]
        )

        source = ColumnDataSource(data={"x": kde_vals, "y": kde_support})
        p.patch("x", "y", source=source, alpha=0.3)

        q1, q3 = series.quantile([0.25, 0.75])
        iqr = q3 - q1
        qmin, q1, q2, q3, qmax = series.quantile([0, 0.25, 0.5, 0.75, 1])
        upper = min(qmax, q3 + 1.5 * iqr)
        lower = max(qmin, q1 - 1.5 * iqr)

        hbar_height = (qmax - qmin) / 500
        whisker_width = 0.1
        box_width = 0.2

        p.rect([kde_x_pos], lower, whisker_width, hbar_height, line_color="grey")
        p.rect([kde_x_pos], upper, whisker_width, hbar_height, line_color="grey")
        p.segment([kde_x_pos], upper, [kde_x_pos], q3, line_color="grey")
        p.segment([kde_x_pos], lower, [kde_x_pos], q1, line_color="grey")
        p.vbar([kde_x_pos], box_width, q2, q3, line_color="black")
        p.vbar([kde_x_pos], box_width, q1, q2, line_color="black")

    p.xaxis.axis_label = group_column
    p.yaxis.axis_label = value_column
    p.xaxis.major_label_orientation = "vertical"

    output_file(output_filename)
    save(p)
    EZChart(plt)


def process_text_files(folder_path):
    """Process text files in the given folder."""
    file_data = []

    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"):
            with open(os.path.join(folder_path, filename), "r") as file:
                content = file.read()
                file_data.append({"Filename": filename, "Content": content})

    if file_data:
        df = pd.DataFrame(file_data)
        return df
    else:
        return None


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "wf-teloseq sequencing report",
        "wf-teloseq",
        args.params,
        args.versions,
        args.workflow_version
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

    with report.add_section("Input summary", "Input summary"):
        tags.p(
            """
            The following plots show the read quality and length distributions as well
            as the base yield of the input data for each sample post filtering for
            read quality and length (use the dropdown menu
            to view the plots for the individual samples).
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

    with report.add_section("Metadata", "Metadata"):
        tabs = Tabs()
        for d in sample_details:
            with tabs.add_tab(d["alias"]):
                df = pd.DataFrame.from_dict(d, orient="index", columns=["Value"])
                df.index.name = "Key"
                DataTable.from_pandas(df)

    with report.add_section(
        "Telomere spanning read stats.", "Telomere spanning read stats."
    ):
        tags.p("""Statistics on reads that have been passed filtering and have been
                identified as spanning the telomere repeat boundary""")
        tabs = Tabs()
        for d, files in sample_files.items():
            with tabs.add_tab(d):
                df = load_csv_with_default(files["telsubstats"], sep="\t")
                DataTable.from_pandas(df, use_index=False)

    # raw and telomere filtered read statistics section
    if not args.mappingreport:
        with report.add_section("Telomere summary", "Telomere summary"):
            tags.p(
                "More stats on reads that have telomere and \
            non-telomere identified"
            )
            tabs = Tabs()
            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    df = load_csv_with_default(files["rawcov"], sep=",", header=0)
                    DataTable.from_pandas(df, use_index=False)

    # end motif section
    # with report.add_section("End motif", "End motif"):
    #    tags.p("End motif parsing is restricted to only reads in which the end of NB
    # barcode01 (default) is very confidently identified, so this will be
    # a subset of your reads")
    #    tabs = Tabs()
    #    for d, files in sample_files.items():
    #        with tabs.add_tab(d):
    #            df = pd.read_csv(files["endmotif"], sep="\t", header=0)
    #            DataTable.from_pandas(df, use_index=False)

    if not args.mappingreport:
        with report.add_section("Telomere length", "Telomere length"):
            tags.p("Unaligned = unmapped telomere reads")
            tabs = Tabs()
            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    df = load_csv_with_default(files["rawlen"], sep=",", header=0)
                    create_boxplot(
                        df,
                        "Telomere_length",
                        "unmapped_filters_boxplot.html",
                        "Unaligned: No Filters",
                    )

    if args.mappingreport:
        with report.add_section("Telomere summary", "Telomere summary"):
            tags.p(
                "Telomere statistics on reads that have telomere and \
            non-telomere identified (row 3 in Read no. section) before \
            and after mapped with different filters applied"
            )
            tags.p("FILTER")
            tags.p("Unaligned = Unaligned telomere reads")
            tags.p(
                "Mapped: LS filters = Low stringency: keep reads where the \
            end mapping position is at least 2000bp (default) beyond last telomere \
            motif or to cut site, which ever comes first. This is \
            to remove short telomere only reads that would not be \
            chr arm specific and also could be truncated."
            )
            tags.p(
                "Mapped: HS filters = High stringency: keep reads where the \
            start mapping position is before last telomere motif identification \
            and end mapping position is within 25 bp of cutsite (with exception of \
            cutsites beyond 45k). This is to ensure reads span subtelomere \
            and to limit mismapping and fragmented reads."
            )
            tags.p("REF_TYPE")
            tags.p(
                "UH = Unknown Haplotype, Mat = Maternal Haplotype, \
            Pat = Paternal Haplotype"
            )
            tags.p("TELOMERE_LENGTH CV")
            tags.p(
                "Coefficient Variation which is how much \
            variation there is relative to the average value. The \
            higher the CV, the greater the level of dispersion around the mean."
            )
            tabs = Tabs()
            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    df = load_csv_with_default(files["output"], sep=",", header=0)
                    # Check if the first column header is "Unnamed: 0" and rename
                    if df.columns[0] == "Unnamed: 0":
                        df.rename(columns={"Unnamed: 0": "Filter"}, inplace=True)
                    DataTable.from_pandas(df, use_index=False)

    if args.mappingreport:
        with report.add_section("Telomere len", "Telomere len"):
            tags.p(
                "Plotting telomere lengths for raw telomere identified \
            reads, and post mapping with different filters."
            )
            # tags.p("No filters = no additional filters")
            tags.p(
                "Mapped: LS filters = Low stringency: keep reads where the \
            end mapping position is at least 2000bp (default) beyond last telomere \
            motif or to cut site, which ever comes first. This is \
            to remove short telomere only reads that would not be \
            chr arm specific and also could be truncated."
            )
            tags.p(
                "Mapped: HS filters = High stringency: keep reads where the \
            start mapping position is before last telomere motif identification \
            and end mapping position is within 25 bp of cutsite (with exception of \
            cutsites beyond 45k). This is to ensure reads span subtelomere \
            and to limit mismapping and fragmented reads."
            )
            tags.p("Unaligned = Unaligned telomere reads")
            tabs = Tabs()

            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    with Grid():
                        df = load_csv_with_default(files["rawlen"], sep=",", header=0)
                        if "Telomere_length" in df.columns:
                            create_boxplot(
                                df,
                                "Telomere_length",
                                "unmapped_filters_boxplot.html",
                                "Unmapped: Raw reads",
                            )
                            # df1 = load_csv_with_default(
                            #     files["none"], sep=",", header=0)
                            # create_boxplot(
                            #     df1, 'Telomere_length',
                            #     "no_filters_boxplot.html", "Mapped: No Filters"
                            # )
                            df2 = load_csv_with_default(
                                files["lenient"], sep=",", header=0
                            )
                            if not df2.empty:
                                create_boxplot(
                                    df2,
                                    "Telomere_length",
                                    "Low_stringency_filters_boxplot.html",
                                    "Mapped: Low stringency",
                                )
                            df3 = load_csv_with_default(
                                files["strict"], sep=",", header=0
                            )
                            if not df3.empty:
                                create_boxplot(
                                    df3,
                                    "Telomere_length",
                                    "High_stringency_filters_boxplot.html",
                                    "Mapped: High stringency",
                                )
                        else:
                            tags.p("No data available")
                        # Raw plot that always have
                        df = load_csv_with_default(files["rawlen"], sep=",", header=0)

    if args.mappingreport:
        with report.add_section("Telomere chr len", "Telomere chr len"):
            tags.p(
                "Plotting telomere lengths after mapping to chromosome arms \
            with different filters."
            )
            tags.p(
                "Mapped: LS filters = Low stringency: keep reads where the \
            end mapping position is at least 2000bp (default) beyond last telomere \
            motif or to cut site, which ever comes first. This is \
            to remove short telomere only reads that would not be \
            chr arm specific and also could be truncated."
            )
            tags.p(
                "Mapped: HS filters = High stringency: keep reads where the \
            start mapping position is before last telomere motif identification \
            and end mapping position is within 25 bp of cutsite (with exception of \
            cutsites beyond 45k). This is to ensure reads span subtelomere \
            and to limit mismapping and fragmented reads."
            )
            tabs = Tabs()
            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    tabs2 = Tabs()
                    with tabs2.add_tab(
                        "Telomere length per \
                                        chromosome (High stringency Filters)"
                    ):
                        df3 = load_csv_with_default(files["strict"], sep=",", header=0)
                        if not df3.empty:
                            create_chrboxplot(
                                df3,
                                "Ref",
                                "Telomere_length",
                                "Chromosome_strict_boxplot.html",
                                "Chromosome: High stringency filters",
                            )
                        else:
                            tags.p("No data available")
                    with tabs2.add_tab(
                        "Telomere length per chromosome \
                                        (Low stringency Filters)"
                    ):
                        df3 = load_csv_with_default(files["lenient"], sep=",", header=0)
                        if not df3.empty:
                            create_chrboxplot(
                                df3,
                                "Ref",
                                "Telomere_length",
                                "Chromosome_lenient_boxplot.html",
                                "Chromosome: Lenient Filters",
                            )
                        else:
                            tags.p("No data available")

    if args.mappingreport:
        with report.add_section("Coverage (chr)", "Coverage (chr)"):
            tags.p(
                "Tables of telomere lengths after mapping to chromosome arms \
            with different filters. Mapq score filter is applied to each \
            of these conditions (default 4)"
            )
            tags.p(
                "Mapped: LS filters = Low stringency: keep reads where the \
            end mapping position is at least 2000bp (default) beyond last telomere \
            motif or to cut site, which ever comes first. This is \
            to remove short telomere only reads that would not be \
            chr arm specific and also could be truncated."
            )
            tags.p(
                "Mapped: HS filters = High stringency: keep reads where the \
            start mapping position is before last telomere motif identification \
            and end mapping position is within 25 bp of cutsite (with exception of \
            cutsites beyond 45k). This is to ensure reads span subtelomere \
            and to limit mismapping and fragmented reads."
            )
            tabs = Tabs()

            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    tabs2 = Tabs()
                    with tabs2.add_tab(
                        "Telomere chromosome coverage (High stringency Filter)"
                    ):
                        df = load_csv_with_default(
                            files["strictcov"], sep=",", header=0
                        )
                        DataTable.from_pandas(df, use_index=False)
                    with tabs2.add_tab(
                        "Telomere chromosome coverage (Low stringency Filter)"
                    ):
                        df = load_csv_with_default(
                            files["lenientcov"], sep=",", header=0
                        )
                        DataTable.from_pandas(df, use_index=False)

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
        help="The verison of the workflow run to generate this report.",
    )
    parser.add_argument(
        "--mappingreport",
        action="store_true",
        help="raw input reads telomere read lengths for boxplot (lenient)",
    )
    return parser
