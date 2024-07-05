#!/usr/bin/env python
"""Script to process error reads data and print ID results to a file."""

import sys

import numpy as np
import pandas as pd

input_file = sys.argv[1]
input_file2 = sys.argv[2]

# Read the data from your file into a pandas DataFrame
df = pd.read_csv(
    input_file,
    delimiter='\t',
    header=None,
    names=['seqID', 'patternName', 'pattern', 'strand', 'start', 'end', 'matched'],
    skiprows=1
)

# Import reference fai as dataframe for list of chr that should be in list
rawreadset = pd.read_csv(input_file2, sep="\t", header=None, names=['seqID', 'start'])

# Merge rawreadset by seqID
df = df.merge(rawreadset, on='seqID')

# Convert the 'start' columns to numeric
df['start_x'] = pd.to_numeric(df['start_x'], errors='coerce')
df['start_y'] = pd.to_numeric(df['start_y'], errors='coerce')

# Remove rows where start_x is greater than start_y as errors only occur within telomere
df = df[df['start_x'] <= df['start_y']]

# Group the DataFrame by 'seqID'
grouped = df.groupby('seqID')

# Define the sliding window size (500 base pairs)
window_size = 500

# Prepare to write output to a file
output_file = 'errors.txt'
with open(output_file, 'w') as f:
    # Iterate through groups and create a sliding window to count rows
    for name, group in grouped:
        start_positions = group['start_x'].to_numpy()
        for start in start_positions:
            count = (np.abs(start_positions - start) <= window_size).sum()
            if count >= 5:
                f.write(f'{name}\n')
