"""Script for clustering and processing sequence data based on various criteria."""
# TODO: the magic numbers in here should be exposed to the user
import os

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def dynamic_group_values2(df, column, window=11):
    """Perform dynamic grouping on the dataframe based on the specified column."""
    df = df.sort_values(by=column)
    group = 0
    df[f'{column}_group'] = 0
    previous_condition = None

    for i in range(window - 1, len(df)):
        window_vals = df[column].iloc[i - window + 1:i + 1]
        mean = np.mean(window_vals)
        diff = abs(window_vals.iloc[-1] - window_vals.iloc[0])
        condition = diff < (mean * 0.1)

        if previous_condition is None or condition != previous_condition:
            group += 1

        df[f'{column}_group'].iloc[i] = group
        previous_condition = condition
    return df


def indel_group_values2(df, column):
    """Create a new column 'indel_group' based on presence of indels."""
    df['indel_group'] = np.where(df[column] > 0, 1, 0)
    return df


def process_group(group_df, min_cluster_size):
    """Process each group of data to filter based on size and criteria."""
    group_sizes = group_df.groupby('final_grouping')['Read'].count()
    valid_groups = group_sizes[group_sizes >= min_cluster_size].index
    filtered_df = group_df[group_df['final_grouping'].isin(valid_groups)]

    if len(filtered_df) > 0:
        return filtered_df[['Read', 'Ref', 'final_grouping']]

    return pd.DataFrame()


def main(args):
    """Run the entry point."""
    # Load dataframes
    df1 = pd.read_csv(args.input_seqkitbam, delimiter='\t')
    cols = ['Read', 'motif']
    df2 = pd.read_csv(args.input_motif, names=cols, header=None, delimiter='\t')
    cols2 = ['Read', 'indelnumber']
    df_indels = pd.read_csv(args.indel_counts, names=cols2, header=None, delimiter='\t')

    # Merge dataframes on read id
    intermediate_df = pd.merge(df1, df2, on='Read', how='left')
    merged_df = pd.merge(intermediate_df, df_indels, on='Read', how='left')

    # Exclude rows where 'Acc' is above 99.8
    merged_df = merged_df[merged_df['Acc'] <= 99.8]

    # Filter merged_df
    df_filtered3 = merged_df[merged_df['RefCov'] > 0.30]
    df_filtered2 = df_filtered3[df_filtered3['IsSup'] == 0]
    df_filtered = df_filtered2[df_filtered2['IsSec'] == 0]

    # Group by Ref
    ref_groups = df_filtered.groupby('Ref')

    # List to hold processed DataFrames
    grouped_dfs = []

    # Iterate over each group
    for ref, group in ref_groups:
        group = indel_group_values2(group, 'indelnumber')
        grouped_dfs.append(group)

    # Concatenate all processed groups
    df_final = pd.concat(grouped_dfs)

    # Generate final_grouping column
    df_final['final_grouping'] = (
        df_final['Ref'].astype(str) + '_round2_' + df_final['indel_group'].astype(str)
    )

    # Process groups and generate results
    result = (
        df_final.groupby("Ref")
        .apply(process_group, args.min_cluster_size)
        .reset_index(drop=True)
    )

    # Get base name of input file
    base_name = os.path.basename(args.input_seqkitbam)
    filename = os.path.splitext(base_name)[0]

    # Output results to CSV files
    result.to_csv(f'{filename}_final_grouping_of_Reads.csv', index=False)
    df_final.to_csv(f'{filename}_all_information.csv', index=False)

    # Group result df by final_grouping and write IDs to text files
    grouped = result.groupby('final_grouping')
    for name, group in grouped:
        out_path = os.path.join('', f'{name}.txt')
        group['Read'].to_csv(out_path, index=False, header=False)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("cluster2b")
    parser.add_argument(
        "input_seqkitbam",
        help="`seqkit bam` TSV stats file",
    )
    parser.add_argument(
        "input_motif",
        help="file with tab-separated alignment motifs and counts",
    )
    parser.add_argument(
        "indel_counts",
        help="headerless TSV file with columns: contig, read_id, indel_count",
    )
    parser.add_argument(
        "min_cluster_size",
        type=int,
        help="minimum cluster size",
    )
    return parser
