"""Create workflow report."""
import json
import os
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


def gather_sample_files(sample_details, args):
    """Create dictionary for sample with paths to data files."""
    sample_files = {}
    for sample in sample_details:
        sample_dir = os.path.join(args.data[0], sample["alias"])
        sample_files[sample["alias"]] = {
            "fastcat_stats": os.path.join(sample_dir, "fastcat_stats"),
            "rawcov": os.path.join(sample_dir, "Sample_raw_Coverage.csv"),
            "rawlen": os.path.join(
                sample_dir, "Sample_raw_Per_Read_telomere_length.csv"),
            "rawstats": os.path.join(sample_dir, "raw.txt"),
            "telstats": os.path.join(sample_dir, "telomere.txt"),
            "telsubstats": os.path.join(sample_dir, "telomere_subtelomere.txt"),
            "subtellen": os.path.join(sample_dir, "subtelomere.txt")
        }
    return sample_files


def gather_sample_files_mapping(sample_details, args):
    """Create dictionary for sample with paths to data files for mapping."""
    sample_files = {}
    for sample in sample_details:
        sample_dir = os.path.join(args.data[0], sample["alias"])
        sample_files[sample["alias"]] = {
            "fastcat_stats": os.path.join(sample_dir, "fastcat_stats"),
            "endmotif": os.path.join(sample_dir, "summary_telomere_motif.txt"),
            "rawcov": os.path.join(sample_dir, "Sample_raw_Coverage.csv"),
            "rawlen": os.path.join(
                sample_dir, "Sample_raw_Per_Read_telomere_length.csv"),
            "none": os.path.join(
                sample_dir, "nofiltered_Per_Read_telomere_length.csv"),
            "lenient": os.path.join(
                sample_dir, "lowfiltered_Per_Read_telomere_length.csv"),
            "strict": os.path.join(
                sample_dir, "highfiltered_Per_Read_telomere_length.csv"),
            "strictcov": os.path.join(
                sample_dir, "highfiltered_chr_arm_Coverage.csv"),
            "lenientcov": os.path.join(
                sample_dir, "lowfiltered_chr_arm_Coverage.csv"),
            "rawstats": os.path.join(sample_dir, "raw.txt"),
            "telstats": os.path.join(sample_dir, "telomere.txt"),
            "telsubstats": os.path.join(
                sample_dir, "telomere_subtelomere.txt"),
            "output": os.path.join(sample_dir, "output.csv"),
            "subtellen": os.path.join(sample_dir, "subtelomere.txt")
        }
    return sample_files


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
    kde_support = np.hstack([
        np.linspace(series.min(), series.max(), 100),
        np.linspace(series.max(), series.min(), 100)
    ])

    source = ColumnDataSource(data={'x': kde_vals, 'y': kde_support})
    p.patch('x', 'y', source=source, alpha=0.3)

    padding_top = 10
    p.y_range = Range1d(
        start=series.min() - padding_top,
        end=series.max() + padding_top
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
    p.xaxis.axis_label = column_name
    p.yaxis.axis_label = 'Values'

    output_file(output_filename)
    save(p)
    EZChart(plt)


def create_chrboxplot(df, group_column, value_column, output_filename, chart_name):
    """Create a boxplot for chromosome data."""
    grouped = df.groupby(group_column)

    def sort_key(item):
        match = re.search(r'\w+_(\d+)', item)
        if not match:
            match = re.search(r'chr(\d+)', item, re.IGNORECASE)

        if match:
            return (1, int(match.group(1)))
        else:
            return (0, item)

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
        kde_support = np.hstack([
            np.linspace(series.min(), series.max(), 100),
            np.linspace(series.max(), series.min(), 100)
        ])

        source = ColumnDataSource(data={'x': kde_vals, 'y': kde_support})
        p.patch('x', 'y', source=source, alpha=0.3)

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
            with open(os.path.join(folder_path, filename), 'r') as file:
                content = file.read()
                file_data.append({'Filename': filename, 'Content': content})

    if file_data:
        df = pd.DataFrame(file_data)
        return df
    else:
        return None


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "wf-teloseq sequencing report", "wf-teloseq",
        args.params, args.versions, "0.0.3")

    with open(args.metadata) as metadata:
        sample_details = sorted([
            {
                'alias': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ], key=lambda d: d["alias"])

        # Choose sample files based on mapping report flag
        sample_files = (
            gather_sample_files_mapping(sample_details, args)
            if args.mappingreport
            else gather_sample_files(sample_details, args)
        )

    with report.add_section("Input summary", "Input sumamry"):
        tags.p(
            """
            The following plots show the read quality and length distributions as well
            as the base yield of the input data for each sample (use the dropdown menu
            to view the plots for the individual samples).
            """
            # TODO: say that this is post read length filtering once we use fastcat to
            # filter
        )
        sample_ids = tuple(sample_files.keys())
        stats_dirs = tuple(
            sample_files[sample_id]["fastcat_stats"] for sample_id in sample_ids
        )
        SeqSummary(stats_dirs, sample_names=sample_ids)

    with report.add_section("Metadata", "Metadata"):
        tabs = Tabs()
        for d in sample_details:
            with tabs.add_tab(d["alias"]):
                df = pd.DataFrame.from_dict(d, orient="index", columns=["Value"])
                df.index.name = "Key"
                DataTable.from_pandas(df)

    with report.add_section("Read no.", "Read no."):
        tags.p("Statistics on complete set of reads supplied to pipeline, the second \
        row are those reads that have identified telomere (x10 repeats) within the \
        first 60-500bp, and the third row is a further subset row 2 reads that \
        should not have telomere sequence (telomere repeats x4) in the last 70bp \
        as shortest cut site is beyond this so should have at least this amount \
        of sequence at end that is not telomere and don't want telomere only reads.")
        tabs = Tabs()
        for d, files in sample_files.items():
            with tabs.add_tab(d):
                df1 = pd.read_csv(files["rawstats"], sep="\t", header=0)
                df2 = pd.read_csv(files["telstats"], sep="\t", header=None, skiprows=1)
                df2.columns = df1.columns
                df3 = pd.read_csv(
                    files["telsubstats"], sep="\t", header=None, skiprows=1
                )
                df3.columns = df1.columns
                # Concatenate the DataFrames vertically (along rows) into final_df
                df1 = pd.concat([df1, df2, df3], ignore_index=True)
                DataTable.from_pandas(df1, use_index=False)

    # raw and telomere filtered read statistics section
    if not args.mappingreport:
        with report.add_section("Telomere summary", "Telomere summary"):
            tags.p("Telomere statistics on reads that have telomere and \
            non-telomere identified (row 3 in Read no. section)")
            tabs = Tabs()
            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    df = pd.read_csv(files["rawcov"], sep=",", header=0)
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
            tags.p("Raw reads = unmapped telomere reads")
            tabs = Tabs()
            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    df = pd.read_csv(files["rawlen"], sep=",", header=0)
                    create_boxplot(
                        df,
                        'Telomere_length2',
                        "unmapped_filters_boxplot.html",
                        "Unmapped: No Filters")

    if args.mappingreport:
        with report.add_section("Telomere summary", "Telomere summary"):
            tags.p("Telomere statistics on reads that have telomere and \
            non-telomere identified (row 3 in Read no. section) before \
            and after mapped with different filters applied")
            tags.p("Raw reads = unmapped telomere reads")
            tags.p("Mapped: No filters = no additional filters")
            tags.p("Mapped: Lenient Filters = keep reads where the end mapping \
            position is at least 80bp beyond last telomere motif. This is \
            to remove short telomere only reads that would not be chr arm \
            specific and also could be truncated.")
            tags.p("Mapped: Strict Filters = keep reads where the start mapping \
            position is before last telomere motif identification and end \
            mapping position is within 25 bp of cutsite (with exception of \
            cutsites beyond 45k). This is to ensure reads span subtelomere \
            and to limit mismapping and fragmented reads.")
            tabs = Tabs()
            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    df = pd.read_csv(files["output"], sep=",", header=0)
                    DataTable.from_pandas(df, use_index=False)

    with report.add_section("Subtelomere", "Subtelomere"):
        tags.p("Plotting read lengths after trimming the telomere off the reads \
        that have been subsetted as having telomere and non-telomere ends. \
        This can aid in diagnosis of prep issues if high proportion of \
        short subtelomeres which would not be typical of your sample.")
        tabs = Tabs()

        for d, files in sample_files.items():
            with tabs.add_tab(d):
                df = pd.read_csv(files["subtellen"], sep="\t")
                x_range = (1, 20000)
                plt = BokehPlot(
                    tools="save", x_range=x_range,
                    title="Subtelomere Length Histogram"
                )
                p = plt._fig
                # You can adjust the number of bins as needed but exclude outliers
                hist, edges = np.histogram(df, bins=1000)
                source = ColumnDataSource(
                    data={'hist': hist, 'left': edges[:-1], 'right': edges[1:]}
                )
                p.quad(
                    top='hist', bottom=0, left='left', right='right',
                    source=source, line_color="#007BA7",
                    line_width=2, fill_color="#007BA7"
                )
                p.xaxis.axis_label = 'Subtelomere Length'
                p.yaxis.axis_label = 'Read No.'
                EZChart(plt)

    if args.mappingreport:
        with report.add_section("Telomere len", "Telomere len"):
            tags.p("Plotting telomere lengths for raw telomere identified \
            reads, and post mapping with different filters.")
            tags.p("No filters = no additional filters")
            tags.p("Lenient = keep reads where the end mapping position \
            is at least 80bp beyond last telomere motif. This is \
            to remove short telomere only reads that would not be \
            chr arm specific and also could be truncated.")
            tags.p("Strict = keep reads where the start mapping position \
            is before last telomere motif identification and end mapping \
            position is within 25 bp of cutsite (with exception of \
            cutsites beyond 45k). This is to ensure reads span \
            subtelomere and to limit mismapping and fragmented reads.")
            tags.p("Raw reads = unmapped telomere reads")
            tabs = Tabs()

            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    with Grid():
                        df = pd.read_csv(files["rawlen"], sep=",", header=0)
                        create_boxplot(
                            df, 'Telomere_length2',
                            "unmapped_filters_boxplot.html",
                            "Unmapped: Raw reads"
                        )

                        df1 = pd.read_csv(files["none"], sep=",", header=0)
                        create_boxplot(
                            df1, 'Telomere_length2',
                            "no_filters_boxplot.html", "Mapped: No Filters"
                        )

                        df2 = pd.read_csv(files["lenient"], sep=",", header=0)
                        create_boxplot(
                            df2, 'Telomere_length2',
                            "Lenient_filters_boxplot.html",
                            "Mapped: Lenient"
                        )

                        df3 = pd.read_csv(files["strict"], sep=",", header=0)
                        create_boxplot(
                            df3, 'Telomere_length2',
                            "Strict_filters_boxplot.html", "Mapped: Strict"
                        )

    if args.mappingreport:
        with report.add_section("Telomere chr len", "Telomere chr len"):
            tags.p("Plotting telomere lengths after mapping to chromosome arms \
            with different filters.")
            tags.p("Lenient Filters = keep reads where the end mapping position \
            is at least 80bp beyond last telomere motif. This is to remove \
            short telomere only reads that would not be chr arm specific \
            and also could be truncated. (Recommended condition)")
            tags.p("Strict Filters = keep reads where the start mapping position \
            is before last telomere motif identification and end mapping \
            position is within 25 bp of cutsite (with exception of two \
            repeating chr arms). This is to ensure reads span subtelomere \
            and to limit mismapping and fragmented reads.")
            tabs = Tabs()
            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    tabs2 = Tabs()
                    with tabs2.add_tab("Telomere length per \
                                        chromosome (Strict Filters)"):
                        df3 = pd.read_csv(files["strict"], sep=",", header=0)
                        create_chrboxplot(
                            df3, 'Ref', 'Telomere_length2',
                            "Chromosome_strict_boxplot.html",
                            "Chromosome: Strict Filters"
                        )
                    with tabs2.add_tab("Telomere length per chromosome \
                                        (Lenient Filters)"):
                        df3 = pd.read_csv(files["lenient"], sep=",", header=0)
                        create_chrboxplot(
                            df3, 'Ref', 'Telomere_length2',
                            "Chromosome_lenient_boxplot.html",
                            "Chromosome: Lenient Filters"
                        )

    if args.mappingreport:
        with report.add_section("Coverage (chr)", "Coverage (chr)"):
            tags.p("Tables of telomere lengths after mapping to chromosome arms \
            with different filters. Mapq score filter is applied to each \
            of these conditions (default 4)")
            tags.p("No filters = no additional filters")
            tags.p("Lenient = keep reads where the end mapping position is at least \
            80bp beyond last telomere motif. This is to remove short telomere \
            only reads that would not be chr arm specific and also could be \
            truncated. (Recommended condition)")
            tags.p("Strict = keep reads where the start mapping position is before \
            last telomere motif identification and end mapping position is \
            within 25 bp of cutsite (with exception of cutsites beyond 45k). \
            This is to ensure reads span subtelomere and to limit mismapping \
            and fragmented reads.")
            tabs = Tabs()

            for d, files in sample_files.items():
                with tabs.add_tab(d):
                    tabs2 = Tabs()
                    with tabs2.add_tab("Telomere chromosome coverage (Strict Filter)"):
                        df = pd.read_csv(files["strictcov"], sep=",", header=0)
                        DataTable.from_pandas(df, use_index=False)
                    with tabs2.add_tab("Telomere chromosome coverage (Lenient Filter)"):
                        df = pd.read_csv(files["lenientcov"], sep=",", header=0)
                        DataTable.from_pandas(df, use_index=False)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--data", nargs="+", required=True,
        help="Collected outputs per sample")
    parser.add_argument(
        "--mappingreport", action="store_true",
        help="raw input reads telomere read lengths for boxplot (lenient)")
    return parser
