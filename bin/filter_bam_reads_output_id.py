#!/usr/bin/env python
"""Filter and extract data from input files."""

import argparse

import pandas as pd

# Create argument parser
parser = argparse.ArgumentParser(description='Filter and extract data from input')
parser.add_argument('text_file', help='Input text file name')
parser.add_argument('bed_file_cutsites', help='Input bed file cutsites name')
parser.add_argument('bed_file_telomereend', help='Input bed file telomereend name')
parser.add_argument('output_file', help='Output file name')

# Add arguments for filter conditions
parser.add_argument('--lenient', action='store_true', help='Apply filter lenient')
parser.add_argument('--strict', action='store_true', help='Apply filter strict')
parser.add_argument('--none', action='store_true', help='Apply filter none')

args = parser.parse_args()

# Read input files
text_file = pd.read_csv(args.text_file, sep='\t')
bed_file_cutsites = pd.read_csv(
    args.bed_file_cutsites,
    sep='\t',
    header=None,
    names=['Ref', 'StartPos', 'junk'],
    dtype={'Ref': str}
)
bed_file_telomereend = pd.read_csv(
    args.bed_file_telomereend,
    sep='\t',
    header=None,
    names=['Ref', 'telPos1', 'junk2'],
    dtype={'Ref': str}
)

# Merge data
merged_data = pd.merge(text_file, bed_file_cutsites, on='Ref', how='inner')
merged_data2 = pd.merge(merged_data, bed_file_telomereend, on='Ref', how='inner')

# Apply filters based on args
if args.lenient:
    filtered_data = pd.DataFrame()  # Initialize an empty DataFrame for filtered data

    for _, row in merged_data2.iterrows():
        if row['telPos1'] + 80 > row['StartPos']:
            if (row['IsSec'] == 0 and row['IsSup'] == 0 and
                    row['StartPos'] - 25 <= row['EndPos']):
                filtered_data = pd.concat([filtered_data, pd.DataFrame([row])])
        else:
            if ((row['telPos1'] + 80) < row['EndPos'] and
                    row['IsSec'] == 0 and row['IsSup'] == 0):
                filtered_data = pd.concat([filtered_data, pd.DataFrame([row])])

elif args.strict:
    # Create an empty DataFrame to store the filtered data
    filtered_data = pd.DataFrame()

    for _, row in merged_data2.iterrows():
        if row['Ref'] in {"chr21_PATERNAL_P", "chr21_MATERNAL_P"}:
            if (row['Pos'] < row['telPos1'] and
                    row['IsSec'] == 0 and
                    row['IsSup'] == 0):
                filtered_data = pd.concat([filtered_data, pd.DataFrame([row])])
        else:
            if row['StartPos'] - 25 <= row['EndPos']:
                filtered_data = pd.concat([filtered_data, pd.DataFrame([row])])

elif args.none:
    filtered_data = merged_data2[
        (merged_data2['IsSec'] == 0) &
        (merged_data2['IsSup'] == 0)
    ]

# Extract and save output
try:
    output_data = filtered_data[['Read']]
    output_data.to_csv(args.output_file, sep='\t', index=False, header=False)
except KeyError:
    with open(args.output_file, 'w') as f:
        f.write("KeyError: One of the expected columns is missing.")
