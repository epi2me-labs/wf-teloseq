#!/usr/bin/env python
"""
Script to process telomere length coverage raw data.

Then generate summaries and plots.
"""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

input_file2 = sys.argv[1]  # Assuming this is the first dataset
telomere_length2_file = sys.argv[2]  # Second dataset with Telomere_length2
input_basename = os.path.basename(input_file2)


def parse_motif(fname):
    """Parse motif and return DataFrame with columns 'Read' and 'Telomere_length'."""
    cols = {'seqID': str, 'end': int}
    df = pd.read_csv(fname, sep="\t", dtype=cols, usecols=cols.keys())
    df.columns = ["Read", "Telomere_length"]
    return df


def parse_telomere_length2(fname):
    """Parse telomere length, return DataFrame with 'Read' and 'Telomere_length2'."""
    df = pd.read_csv(fname, sep="\t", header=None, names=["Read", "Telomere_length2"])
    return df


# Load the datasets
dfmotif = parse_motif(input_file2)
dftelomere2 = parse_telomere_length2(telomere_length2_file)

# Merge the DataFrames on the Read column
merged_df = pd.merge(dfmotif, dftelomere2, on="Read")

# Define aggregation functions
agg_functions = {
    'Telomere_length2': ['mean', 'std', 'max']
}

# Aggregate statistics
summary_stats = merged_df.agg(agg_functions).round(0)
summary_stats['Read_count'] = len(merged_df)
read_count = len(merged_df)

# Create a new summary DataFrame with the desired format
summary_formatted = pd.DataFrame({
    'Read No': [read_count],
    'Telomere Mean': [summary_stats['Telomere_length2']['mean']],
    'Telomere SD': [summary_stats['Telomere_length2']['std']],
    'Telomere Max': [summary_stats['Telomere_length2']['max']]
})

# Save the formatted summary DataFrame to a CSV
summary_formatted.to_csv(f'{input_basename}_raw_Coverage.csv', index=False)

# Plotting remains the same, adapted to use filtered DataFrame
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
merged_df.boxplot(column='Telomere_length2', fontsize=12, ax=ax2)
ax2.set_ylabel("Telomere Length")
ax2.xaxis.set_label_text("Filtered")
ax2.set_title("Filtered Telomere Length Distribution")
fig2.savefig(f'{input_basename}_raw_Boxplot_of_Telomere_length.pdf')

# Save per-Read information including telomere length
merged_df.to_csv(f'{input_basename}_raw_Per_Read_telomere_length.csv', index=False)
