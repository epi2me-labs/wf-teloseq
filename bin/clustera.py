#!/usr/bin/env python
"""
Script for clustering reads based on sequence and motif information.

Usage:
python clustera.py <input_seqkitbam> <input_motif> <clustersize>
"""

import os
import sys

import numpy as np
import pandas as pd


# Get input and output files from command line
input_seqkitbam = sys.argv[1]
input_motif = sys.argv[2]
clustersize = round(float(sys.argv[3]))


# Import dataframes
df1 = pd.read_csv(input_seqkitbam, delimiter='\t')
cols = ['Read', 'motif']
df2 = pd.read_csv(input_motif, names=cols, header=None, delimiter='\t')

# Merge dataframes on read id
merged_df = pd.merge(df1, df2, on='Read', how='left')


def dynamic_group_values2(df, column, window=7):
    """
    Apply dynamic grouping based on moving window mean and condition.

    Args:
    df (pd.DataFrame): DataFrame containing data.
    column (str): Column name for which grouping is applied.
    window (int): Size of the moving window.

    Returns:
    pd.DataFrame: DataFrame with added group labels.
    """
    df = df.sort_values(by=column)
    group = 0
    df[f'{column}_group'] = 0
    previous_condition = None

    for i in range(window - 1, len(df)):
        window_vals = df[column].iloc[i-window+1:i+1]
        mean = np.mean(window_vals)
        diff = abs(window_vals.iloc[-1] - window_vals.iloc[0])
        condition = diff < (mean * 0.1)

        if previous_condition is None or condition != previous_condition:
            group += 1

        df[f'{column}_group'].iloc[i] = group
        previous_condition = condition

    return df


def process_group(group_df):
    """
    Process each group to filter based on group size.

    Args:
    group_df (pd.DataFrame): DataFrame of grouped data.

    Returns:
    pd.DataFrame: Filtered DataFrame based on group size.
    """
    # Get group sizes
    group_sizes = group_df.groupby('final_grouping')['Read'].count()

    # Filter for size >= 60
    min_size = int(clustersize)
    valid_groups = group_sizes[group_sizes >= min_size].index

    # Filter original df
    filtered_df = group_df[group_df['final_grouping'].isin(valid_groups)]

    if len(filtered_df) > 0:
        return filtered_df[['Read', 'Ref', 'final_grouping']]

    return pd.DataFrame()


df_filtered3 = merged_df[merged_df['RefCov'] > 0.30]
df_filtered2 = df_filtered3[df_filtered3['IsSup'] == 0]
df_filtered = df_filtered2[df_filtered2['IsSec'] == 0]

# Group df by Ref
ref_groups = df_filtered.groupby('Ref')

# List to hold the processed DataFrames
grouped_dfs = []

# Iterate over each group
for ref, group in ref_groups:
    # Apply dynamic_group_values2 to each column of interest
    group = dynamic_group_values2(group, 'Pos')
    group = dynamic_group_values2(group, 'motif')

    # Append the processed group to the list
    grouped_dfs.append(group)

# Concatenate all the processed groups
df_final = pd.concat(grouped_dfs)

# Create final grouping based on columns
df_final['final_grouping'] = df_final['Ref'].astype(str) + '_' \
    + df_final['Pos_group'].astype(str) + '_' + df_final['motif_group'].astype(str)

# Group by final grouping and apply processing
result = df_final.groupby('Ref').apply(process_group).reset_index(drop=True)

# Get base name of input
base_name = os.path.basename(input_seqkitbam)
filename = os.path.splitext(base_name)[0]

# Output results to CSV files
result.to_csv(f'{filename}_final_grouping_of_Reads.csv', index=False)
df_final.to_csv(f'{filename}_all_information.csv', index=False)

# Group result df by final_grouping
grouped = result.groupby('final_grouping')

# Iterate over groups
for name, group in grouped:
    # File path
    out_path = os.path.join('', f'{name}.txt')

    # Write IDs to file
    group['Read'].to_csv(out_path, index=False, header=False)
