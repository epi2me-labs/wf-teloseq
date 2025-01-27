"""This script processes a FASTQ file to find telomere boundaries in sequences."""
import re
import sys

import numpy as np
import pysam

from .util import wf_parser  # noqa: ABS101


# MIT License
#
# Copyright (c) 2024 Ramin Kahidi
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# TODO: lots of magic numbers in here; which ones do we want to expose to the user?
# cjw: none
AREA_DIFFS_THRESHOLD = 0.2
COMPOSITIONCSTRAND = [["CCC", 3 / 6]]


def calc_graph_area(offsets, target_column, window_size):
    """Calculate the graph area for given offsets and window size."""
    data = np.array(offsets)
    transposed_data = data.T
    area_list = []
    row = transposed_data[target_column, :]
    for i in range(0, len(row) - window_size, 1):
        area = row[i: i + window_size].sum()
        area_list.append((area / window_size))
    return area_list


def find_telo_boundary(
    seq,
    composition_c_strand=[],
    telo_window=80,
    window_step=8,
    change_threshold=-20,
    plateau_detection_threshold=-60,
    target_pattern_index=-1,
    nucleotide_graph_area_window_size=700,
    return_last_discontinuity=False,
    secondary_search=True,
):
    """Find the boundary point of the telomere in the given sequence."""
    boundary_point = -1
    nt_offsets = []
    graph_area_window_size = int(nucleotide_graph_area_window_size / window_step)

    if not composition_c_strand:
        composition_c_strand = COMPOSITIONCSTRAND

    composition = composition_c_strand

    for i in range(0, len(seq) - telo_window, window_step):
        telo_seq = seq[i: i + telo_window]
        telo_len = len(telo_seq)
        telo_seq_upper = str(telo_seq.upper())

        current_offsets = []
        for nt_pattern_entry in composition:
            nt_pattern = nt_pattern_entry[0]
            pattern_composition = nt_pattern_entry[1]
            pattern_count = len(re.findall(nt_pattern, telo_seq_upper))
            if is_regex_pattern(nt_pattern):
                if len(nt_pattern_entry) != 3:
                    return -1
                regex_target_length = nt_pattern_entry[2]
                raw_offset_value = (
                    pattern_count * regex_target_length
                ) / telo_len - pattern_composition
                if pattern_composition != 0:
                    percent_offset_value = (
                        raw_offset_value / pattern_composition
                    ) * 100
                    current_offsets.append(percent_offset_value)
                else:
                    percent_offset_value = raw_offset_value * 100
                    current_offsets.append(percent_offset_value)
            else:
                raw_offset_value = (
                    pattern_count * len(nt_pattern)
                ) / telo_len - pattern_composition
                if pattern_composition != 0:
                    percent_offset_value = (
                        raw_offset_value / pattern_composition
                    ) * 100
                    current_offsets.append(percent_offset_value)
                else:
                    percent_offset_value = raw_offset_value * 100
                    current_offsets.append(percent_offset_value)

        nt_offsets.append(current_offsets)

    area_list = calc_graph_area(
        nt_offsets, target_pattern_index, graph_area_window_size
    )
    area_diffs = np.diff(area_list)
    index_at_threshold = -1

    if return_last_discontinuity:
        index_at_threshold = next(
            (
                y
                for y in range(len(area_list) - 2, 0, -1)
                if (
                    area_list[y] > change_threshold
                    and (0 > (area_list[y + 1] - area_list[y]))
                )
            ),
            index_at_threshold,
        )
        if index_at_threshold != -1:
            index_at_threshold = next(
                (
                    y
                    for y in range(index_at_threshold, len(area_list) - 2)
                    if (
                        area_list[y] < plateau_detection_threshold
                        and (0 > (area_list[y + 1] - area_list[y]))
                    )
                ),
                index_at_threshold,
            )
        else:
            index_at_threshold = next(
                (
                    y
                    for y in range(len(area_list) - 2)
                    if (
                        area_list[y] < plateau_detection_threshold
                        and (0 > (area_list[y + 1] - area_list[y]))
                    )
                ),
                index_at_threshold,
            )
    else:
        index_at_threshold = next(
            (
                y
                for y in range(len(area_list) - 2)
                if (
                    area_list[y] < plateau_detection_threshold
                    and (0 > (area_list[y + 1] - area_list[y]))
                )
            ),
            index_at_threshold,
        )

    if index_at_threshold == -1:
        return -1

    area_diffs_threshold = 0.1
    for x in range(index_at_threshold, len(area_diffs) - 1):
        if (
            abs(area_diffs[x]) < area_diffs_threshold
            or area_diffs[x] < area_diffs_threshold
            and area_diffs[x + 1] > area_diffs_threshold
        ):
            boundary_point = x * window_step
            break
    if boundary_point == -1:
        boundary_point = len(area_diffs) * window_step

    if secondary_search:
        telomere_offset = 500
        sub_telomere_offset = 1000
        telomere_offset_re = 30
        sub_telomere_offset_re = 30
        nt_pattern_entry = composition[target_pattern_index]
        nt_pattern = nt_pattern_entry[0]

        lower_index = boundary_point - telomere_offset
        upper_index = boundary_point + sub_telomere_offset
        if lower_index < 0:
            lower_index = 0
        if upper_index > len(seq):
            upper_index = len(seq)
        scan_seq = seq[lower_index:upper_index]
        if len(scan_seq) < 100:
            return boundary_point
        else:
            sec_boundary = find_telo_boundary(
                scan_seq,
                composition_c_strand,
                telo_window=80,
                window_step=8,
                change_threshold=change_threshold,
                plateau_detection_threshold=-60,
                target_pattern_index=-1,
                nucleotide_graph_area_window_size=100,
                return_last_discontinuity=return_last_discontinuity,
                secondary_search=False,
            )

            temp_boundary = boundary_point + sec_boundary - telomere_offset
            if lower_index == 0:
                temp_boundary = sec_boundary
            scan_seq = seq[
                temp_boundary
                - telomere_offset_re: temp_boundary
                + sub_telomere_offset_re
            ]

            pattern_composition = nt_pattern_entry[1]
            if pattern_composition != 1:
                return boundary_point
            else:
                scan_pattern = f"({nt_pattern})({nt_pattern})"
                matches = [match for match in re.finditer(scan_pattern, str(scan_seq))]
                if matches:
                    telo_end = matches[-1].end()
                    boundary_point = temp_boundary + telo_end - telomere_offset_re
                else:
                    boundary_point = temp_boundary

    return boundary_point


def is_regex_pattern(pattern):
    """Check if the given pattern is a regex pattern."""
    regex_chars = r"[\\^$.*+?()[\]{}|]"
    return any(char in pattern for char in regex_chars)


def main(args):
    """Run the entry point."""
    part1 = (
        "TTGAGGGAGTGCATTAGCATACAGGTGCTTGTTACATGTTGAGGGAGTGCATTAGCATACAGGTGCTTGTTACATGTT"
    )
    non_telo_seq = part1 * 10

    with pysam.FastxFile(args.fastq_file) as f:
        for record in f:
            # Append the known non-telomere sequence to both ends of the sequence
            modified_seq = record.sequence + non_telo_seq
            boundary = find_telo_boundary(modified_seq, secondary_search=True)
            sys.stdout.write(f"{record.name}\t{boundary}\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("TeloBoun")
    parser.add_argument("fastq_file", help="input FASTQ file")
    return parser
