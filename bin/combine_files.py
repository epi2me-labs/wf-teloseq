#!/usr/bin/env python
"""Script to combine statistics CSV."""

import csv

headers = [
    "Filter", "Read count", "Accuracy mean", "Accuracy std", "Telomere length mean",
    "Telomere length std", "Telomere length max", "Telomere length weighted mean",
    "Telomere length weighted std"
]

output_rows = []

file_names = {
    "raw_coverage.csv": "Raw Reads",
    "nofiltered.csv": "Mapped: No Filter",
    "lowfiltered.csv": "Mapped: Lenient Filter",
    "highfiltered.csv": "Mapped: Strict Filter"
}

for f, header in file_names.items():
    rows = []
    with open(f, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        if f != "raw_coverage.csv":
            next(reader)
        row = next(reader)

        if f == "raw_coverage.csv":
            row = [row[0], "N/A", "N/A"] + [row[1]] + [row[2]] + \
                    [row[3]] + ["N/A", "N/A"]
        rows.append(row)

    output_rows.append([header] + rows[0])

with open('output.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(headers)
    writer.writerows(output_rows)
