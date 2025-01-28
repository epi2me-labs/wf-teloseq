"""Script: Name de novo contigs from BLAST output to pangenome."""

from collections import Counter
import re

from Bio import SeqIO
import pandas as pd
import pysam

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entry point."""
    parser = wf_parser("RenameFa")
    parser.add_argument(
        "blast",
        help="Input BLAST file",
    )
    parser.add_argument(
        "reference_file",
        help="Reference FASTA file",
    )
    parser.add_argument(
        "output_fasta",
        help="Output FASTA file with renamed contigs",
    )
    parser.add_argument(
        "output_summary",
        help="Output summary CSV file",
    )
    parser.add_argument(
        "motif_threshold",
        type=int,
        default=40,
        help="Motif count threshold for chr21_p naming (default: 40)",
    )
    parser.add_argument(
        "max_sstart",
        type=int,
        default=20,
        help="Maximum sstart value for filtering (default: 20)"
    )
    parser.add_argument(
        "min_pident",
        type=float,
        default=95.0,
        help="Minimum percent identity for filtering (default: 95.0)"
    )
    parser.add_argument(
        "special_chr_arms",
        nargs='+',
        default=[
            'chr9q', 'chr11p', 'chr13p', 'chr15p', 'chr19p',
            'chr6p', 'chr20q', 'chr3q', 'chr6q', 'chr7p', 'chr14p'
        ],
        help="List of special chromosome arms with specific naming rules",
    )
    return parser


def extract_chr_arm(sseqid):
    """Extract chromosome and arm information from sseqid."""
    parts = sseqid.split('_')
    for i, part in enumerate(parts):
        if part.startswith('chr'):
            chr_num = part  # e.g., 'chr16'
            # Check if the preceding part is 'p' or 'q'
            if i > 0 and parts[i - 1] in ['p', 'q']:
                arm = parts[i - 1]
                chr_arm = f"{chr_num}{arm}"
                return chr_arm
            else:
                # If no 'p' or 'q' preceding, default to chromosome without arm
                return chr_num
    return None


def assign_contig_names(df_filtered, special_chr_arms):
    """Assign names to contigs based on the BLAST hit information."""
    # Create a dictionary to store assigned names and counts
    assigned_names = {}
    naming_reasons = {}
    generic_counter = 1

    # This dictionary will store a list of (qseqid, assigned_name) for each chr_arm.
    chr_arm_assignments = {}

    df_filtered['chr_arm'] = df_filtered['sseqid'].apply(extract_chr_arm)

    for qseqid, group in df_filtered.groupby('qseqid'):
        top_hits = group.head(5)

        # Only count valid chr_arm values and remove any None values
        valid_chr_arms = [
            arm for arm in top_hits['chr_arm'].tolist()
            if arm is not None
        ]
        chr_arm_counts = Counter(valid_chr_arms)

        new_name = None
        reason = None

        if chr_arm_counts:
            most_common_arm, most_common_count = chr_arm_counts.most_common(1)[0]

            # Check for equal pident in top 2 hits first
            if len(top_hits) >= 2:
                top_pidents = top_hits['pident'].head(2).values
                top_chr_arms = top_hits['chr_arm'].head(2).values
                if (
                    round(top_pidents[0], 3) == round(top_pidents[1], 3)
                    and top_chr_arms[0] != top_chr_arms[1]
                ):
                    new_name = f"contig_{generic_counter}"
                    reason = (
                        f"Top 2 hits have equal percent identity and different chr "
                        f"({top_pidents[0]:.3f}%)"
                    )
                    generic_counter += 1

            # If no equal pident match, check chromosome arm rules
            if new_name is None:
                if most_common_count >= 2:
                    if most_common_arm in special_chr_arms:
                        new_name = f"contig_{generic_counter}"
                        reason = f"{most_common_arm} in {special_chr_arms}"
                        generic_counter += 1
                    else:
                        assigned_names[most_common_arm] = (
                            assigned_names.get(most_common_arm, 0)
                            + 1
                        )
                        suffix_count = assigned_names[most_common_arm]
                        suffix = chr(96 + suffix_count)  # 'a' = 97 in ASCII
                        new_name = f"{most_common_arm}_{suffix}"
                        reason = (
                            f"At least 2 out of 5 top hits match {most_common_arm}"
                            f" ({most_common_count} matches)"
                        )

                        # Track this assignment
                        if most_common_arm not in chr_arm_assignments:
                            chr_arm_assignments[most_common_arm] = []
                        chr_arm_assignments[most_common_arm].append((qseqid, new_name))

                        # If we have assigned a third contig to the same arm,
                        # Revert all three to generic
                        if len(chr_arm_assignments[most_common_arm]) == 3:
                            # Revert all three assignments
                            for assigned_qseqid, old_assigned_name in (
                                chr_arm_assignments[most_common_arm]
                            ):
                                # Find those rows in df_filtered and update
                                mask = (df_filtered['qseqid'] == assigned_qseqid)
                                df_filtered.loc[mask, 'new_name'] = (
                                    f"contig_{generic_counter}"
                                )
                                naming_reasons[assigned_qseqid] = (
                                    "More than two contigs matched the same chr_arm, \
                                    reverted to generic naming"
                                )
                                generic_counter += 1

                            # Reset assignments for this arm
                            chr_arm_assignments[most_common_arm] = []
                            assigned_names[most_common_arm] = 0

                            # The current qseqid was also part of the reverted set,
                            # So don't assign again.
                            new_name = df_filtered.loc[group.index, 'new_name'].iloc[0]
                            reason = naming_reasons[qseqid]

        # Default generic naming if no rules matched
        if new_name is None:
            new_name = f"contig_{generic_counter}"
            reason = "No specific naming rule matched - assigned generic name"
            generic_counter += 1

        # Store the name and reason if not already stored by the revert code
        if qseqid not in naming_reasons:
            naming_reasons[qseqid] = reason
        df_filtered.loc[group.index, 'new_name'] = new_name

    return df_filtered, naming_reasons, generic_counter


def filter_dataframe(df, max_sstart, min_pident):
    """Filter the DataFrame based on given criteria."""
    df_filtered = df[
        (df['sstart'] < max_sstart) &
        (df['pident'] > min_pident) &
        (~df['sseqid'].str.contains('p_chr21'))
    ]

    return df_filtered


def update_names_based_on_sequence(
    df_top_hits,
    fasta_file,
    naming_reasons,
    motif_threshold,
):
    """Update contig names based on sequence content."""
    contig_sequences = {
        record.id: str(record.seq)
        for record in SeqIO.parse(fasta_file, "fasta")
    }

    chr21_p_count = 0
    for idx, row in df_top_hits.iterrows():
        qseqid = row['qseqid']
        sequence = contig_sequences.get(qseqid, '')
        motif_count = sequence.count('GGAATGGAAT')

        if motif_count > motif_threshold:
            chr21_p_count += 1
            new_name = f"chr21p_{chr(96 + chr21_p_count)}"
            df_top_hits.loc[idx, 'new_name'] = new_name
            naming_reasons[qseqid] = f"Contains >{motif_threshold} GGAATGGAAT motifs"

    return df_top_hits, naming_reasons


def create_summary_file(df_top_hits, df, summary, naming_reasons, generic_counter):
    """Create a summary CSV file detailing the renaming process."""
    summary_data = []

    # Get list of actual contigs from the FASTA file
    processed_qseqids = set()

    # Process all qseqids from the actual data (excluding header)
    for qseqid in df['qseqid'].unique():
        # Skip if this is the header row
        if qseqid == 'qseqid':
            continue

        processed_qseqids.add(qseqid)

        # Get the new name either from df_top_hits or assign a generic name
        if qseqid in df_top_hits['qseqid'].values:
            new_name = df_top_hits[df_top_hits['qseqid'] == qseqid]['new_name'].iloc[0]
            reason = naming_reasons.get(qseqid, "Unknown reason")
        else:
            new_name = f"contig_{generic_counter}"
            reason = "Not present in filtered BLAST results - assigned generic name"
            generic_counter += 1

        summary_data.append({
            'qseqid': qseqid,
            'new_name': new_name,
            'reason': reason
        })

    # Create DataFrame from summary data
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary, index=False)

    return summary_df


def sort_fasta(fasta_file, summary_df):
    """Sort the FASTA file based on the new contig names in summary."""
    # Create the rename dictionary from summary DataFrame
    rename_dict = dict(zip(summary_df['qseqid'], summary_df['new_name']))

    # Load and rename sequences
    sequences = []
    with pysam.FastxFile(fasta_file) as fh:
        for entry in fh:
            if entry.name in rename_dict:
                entry.name = rename_dict[entry.name]
                sequences.append(entry)

    # Sort sequences
    def sort_key(name):
        chr_match = re.match(r'chr(\d+|X|Y)', name)
        contig_match = re.match(r'contig_(\d+)', name)

        if chr_match:
            chr_val = chr_match.group(1)
            if chr_val.isdigit():
                return (0, int(chr_val))
            elif chr_val == 'X':
                return (0, float('inf') - 1)
            elif chr_val == 'Y':
                return (0, float('inf'))
        elif contig_match:
            return (1, int(contig_match.group(1)))
        return (2, float('inf'))

    sequences_sorted = sorted(sequences, key=lambda record: sort_key(record.name))

    return sequences_sorted


def write_sorted_fasta(sequences, output_fasta_file):
    """Write sorted sequences to the output FASTA file."""
    with open(output_fasta_file, 'w') as out_handle:
        for seq in sequences:
            out_handle.write(f">{seq.name}\n{seq.sequence}\n")


def main(args):
    """Identify useful blast hits for renaming."""
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]

    # Read BLAST data, skipping the header row since we provide our own column names
    df = pd.read_csv(args.blast, sep=r'\s+', names=columns, engine='python', skiprows=1)

    df['sseqid'] = df['sseqid'].astype(str)
    df['sstart'] = pd.to_numeric(df['sstart'], errors='coerce')
    df['pident'] = pd.to_numeric(df['pident'], errors='coerce')

    # Assign contig names and get naming reasons
    df_filtered = filter_dataframe(df, args.max_sstart, args.min_pident)

    # Assign contig names and get naming reasons
    df_top_hits, naming_reasons, generic_counter = assign_contig_names(
        df_filtered, args.special_chr_arms)

    # Update names based on sequence content
    df_top_hits, naming_reasons = update_names_based_on_sequence(
        df_top_hits, args.reference_file, naming_reasons, args.motif_threshold)

    # Create summary file
    summary_df = create_summary_file(
        df_top_hits, df, args.output_summary, naming_reasons, generic_counter)

    # Sort and write FASTA file
    sorted_sequences = sort_fasta(args.reference_file, summary_df)

    write_sorted_fasta(sorted_sequences, args.output_fasta)
