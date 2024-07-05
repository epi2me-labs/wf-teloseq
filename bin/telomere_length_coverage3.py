#!/usr/bin/env python
"""Script to analyze telomere lengths and coverage data."""

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

input_file = sys.argv[1]
input_file2 = sys.argv[2]
input_file3 = sys.argv[3]
coverage_file_path = sys.argv[4]
tel_length = sys.argv[5]

# Read the first line and strip newline characters
with open(coverage_file_path, 'r') as file:
    coverage_value_str = file.readline().strip()

# Convert the string value to float
min_coverage = float(coverage_value_str)


def parse_seqkit(fname):
    """Parse seqkit output file and return a DataFrame."""
    cols = {
        'Read': str, 'Ref': str, 'Acc': float, 'ReadLen': int, 'RefCov': float,
        'EndPos': int, 'RightClip': int, 'LeftClip': int, 'Strand': str,
        'ReadAln': int, 'ReadCov': float, 'IsSec': bool, 'IsSup': bool
    }
    df = pd.read_csv(fname, sep="\t", dtype=cols, usecols=cols.keys())
    df['Clipped'] = df['ReadLen'] - df['ReadAln']
    df['Type'] = 'Primary'
    df.loc[df['IsSec'], 'Type'] = 'Secondary'
    df.loc[df['IsSup'], 'Type'] = 'Supplementary'
    return df


def parse_motif(fname):
    """Parse telomere motif file and return a DataFrame."""
    cols = {'seqID': str, 'end': int}
    df = pd.read_csv(fname, sep="\t", dtype=cols, usecols=cols.keys())
    df.columns = ["Read", "Telomere_length"]
    return df


def parse_telomere_length2(fname):
    """Parse telomere length file and return a DataFrame."""
    df = pd.read_csv(fname, sep="\t", header=None, names=["Read", "Telomere_length2"])
    return df


# Parse input files
dfreads = parse_seqkit(input_file)
dfmotif = parse_motif(input_file2)
dftel = parse_telomere_length2(tel_length)

# Filter primary alignments and specific strand
rslt_df2 = dfreads[dfreads['Type'] == 'Primary']
rslt_df = rslt_df2[rslt_df2['Strand'] == '1']

# Merge dataframes
dfmerged = rslt_df.merge(dfmotif, on='Read')
dfmerged2 = dfmerged.merge(dftel, on='Read')

# Group by reference and aggregate statistics
Summary = dfmerged2.groupby('Ref').agg(**{
    'Coverage': ('Read', 'count'),
    'avg_accuracy': ('Acc', 'mean'),
    'SD_accuracy': ('Acc', 'std'),
    'Telomere_mean': ('Telomere_length2', 'mean'),
    'Telomere_sd': ('Telomere_length2', 'std'),
    'Telomere_max': ('Telomere_length2', 'max')
})

# Format specific columns to two decimal places
Summary['Coverage'] = Summary['Coverage'].round(0).astype(int)
Summary['avg_accuracy'] = Summary['avg_accuracy'].round(2)
Summary['SD_accuracy'] = Summary['SD_accuracy'].round(2)
Summary['Telomere_mean'] = Summary['Telomere_mean'].round(2)
Summary['Telomere_sd'] = Summary['Telomere_sd'].round(2)
Summary['Telomere_max'] = Summary['Telomere_max'].round(0).astype(int)

# Import reference fai as dataframe for list of chr that should be in list
cols = {'Ref': str, 'length': int, 'length2': int, 'length3': int, 'length4': int}
referenceset = pd.read_csv(
    input_file3, sep="\t", header=None,
    usecols=cols.keys(), names=cols.keys(), dtype=cols
)

# Merge with Summary dataframe
Summary2 = referenceset['Ref'].to_frame().merge(Summary, on='Ref', how='left').fillna(0)

# Filter Summary dataframe by minimum coverage
Summary2['Coverage'] = Summary2['Coverage'].astype(int)
filtered_summary = Summary2[Summary2['Coverage'] >= min_coverage]

# Save filtered summary to CSV
input_basename = os.path.basename(input_file)
filtered_summary.to_csv(f'{input_basename}_chr_arm_Coverage.csv', index=False)

# Remove low coverage contigs from merged dataframe
dfmerged_corrected = dfmerged2[dfmerged2['Ref'].isin(filtered_summary['Ref'])]

# Save per-read information including telomere length to CSV
dfmerged_corrected.to_csv(
    f'{input_basename}_Per_Read_telomere_length.csv',
    index=False
)

# Summarize read coverage for each chr arm with boxplot
# FIGURES
fig = plt.figure()
fig.set_figheight(15)
fig.set_figwidth(20)
ax1 = fig.add_subplot(111)
dfmerged_corrected.boxplot(
    column='Telomere_length', by='Ref', fontsize=12, ax=ax1, rot=90
)
ax1.xaxis.set_label_text("Chr arm")
fig.subplots_adjust(hspace=0.5)
ax1.set_ylabel("Telomere Length")
ax1.set_title("Telomere Length Distribution")
fig.savefig(f'{input_basename}_chr_arm_Boxplot_of_Telomere_length.pdf')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
dfmerged_corrected.boxplot(column='Telomere_length', fontsize=12, ax=ax2)
ax2.set_ylabel("Telomere Length")
ax2.xaxis.set_label_text(input_basename)
ax2.set_title("Telomere Length Distribution")
fig2.savefig(f'{input_basename}_Boxplot_of_Telomere_length.pdf')

# Aggregate statistics for overall summary
agg_count = dfmerged_corrected['Read'].count()
agg_mean_acc = dfmerged_corrected['Acc'].mean()
agg_std_acc = dfmerged_corrected['Acc'].std()
agg_mean_telomere = dfmerged_corrected['Telomere_length'].mean()
agg_std_telomere = dfmerged_corrected['Telomere_length'].std()
agg_max_telomere = dfmerged_corrected['Telomere_length'].max()

# Create summary DataFrame
Summary4 = pd.DataFrame({
    ('Read', 'count'): agg_count,
    ('Acc', 'mean'): agg_mean_acc,
    ('Acc', 'std'): agg_std_acc,
    ('Telomere_length', 'mean'): agg_mean_telomere,
    ('Telomere_length', 'std'): agg_std_telomere,
    ('Telomere_length', 'max'): agg_max_telomere
}, index=[0])

# Format specific columns to two decimal places
columns_to_round = [
    ('Acc', 'mean'),
    ('Acc', 'std'),
    ('Telomere_length', 'mean'),
    ('Telomere_length', 'std'),
    ('Telomere_length', 'max')
]

for column in columns_to_round:
    Summary4[column] = Summary4[column].round(2)

# Calculate overall weighted average
overall_weighted_mean = Summary2['Telomere_mean'].mean().round(2)

# Add weighted mean to the summary DataFrame
new_column = pd.MultiIndex.from_tuples([('Telomere length', 'weighted mean')])
Summary4[new_column] = overall_weighted_mean

# Calculate standard deviation of weighted means
weighted_std = Summary2['Telomere_mean'].std().round(2)

# Add std of weighted mean to the summary DataFrame
new_column = pd.MultiIndex.from_tuples([('Telomere length', 'std weighted mean')])
Summary4[new_column] = weighted_std

# Save the summary DataFrame to CSV
Summary4.to_csv(f'{input_basename}.csv', index=False)
