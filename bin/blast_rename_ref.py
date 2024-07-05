#!/usr/bin/env python
"""Script to perform BLAST-based renaming of reference sequences."""

import argparse
from collections import defaultdict
import os
import tempfile

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd

# Collect inputs
parser = argparse.ArgumentParser()
parser.add_argument('-db', required=True, help="DB FASTA")
parser.add_argument('-q', required=True, help="Query multi FASTA")
parser.add_argument('-i', required=True, help="Input CSV with BLAST DB info")
parser.add_argument('-o', default='renamed.fa', help="Output FASTA")
args = parser.parse_args()

db_fasta = args.db
query_fasta = args.q
info_csv = args.i
out_fasta = args.o

# Read CSV and create mapping
db_info = pd.read_csv(info_csv)
name_to_group = dict(zip(db_info.iloc[:, 0], db_info.iloc[:, 6]))

# Define thresholds
min_ident = 95
min_alignment = 80

# BLAST params
db = "mydb"
outfmt = 6

# New ID tracker
new_id_tracker = defaultdict(list)
no_hit_counter = 26  # Starting group number for no hits

# 1. Make BLAST DB
cmd = f"makeblastdb -in {db_fasta} -dbtype nucl -out {db}"
os.system(cmd)

# 2. BLAST each query
for seq in SeqIO.parse(query_fasta, "fasta"):
    with tempfile.NamedTemporaryFile(suffix=".fasta") as temp:
        SeqIO.write([seq], temp.name, "fasta")
        blast_out = f"{seq.id}.tsv"

        # BLAST
        cline = NcbiblastnCommandline(
            query=temp.name,
            db=db, outfmt=outfmt, out=blast_out
        )
        stdout, stderr = cline()

    if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
        # 3. Load BLAST table
        table = pd.read_csv(blast_out, sep="\t", header=None)

        # 4. Process BLAST hits and get top hit
        top_hits = table.sort_values([3, 2], ascending=[0, 0]).groupby(0).head(1)
        top = top_hits[top_hits[0] == seq.id]

        if len(top_hits) > 0:
            identity = top_hits[2].values[0]
            aln_len = top_hits[3].values[0]

            if identity >= min_ident and aln_len / len(seq.seq) >= min_alignment / 100:

                top_hit_name = top[1].values[0]
                group_number = name_to_group.get(top_hit_name, no_hit_counter)
                new_id = f'Group_{group_number}_{top_hit_name}'
                new_id_tracker[new_id].append(seq)
            else:
                if not seq.id.startswith("Group"):
                    new_id = f'Group_{no_hit_counter}_contig_{len(new_id_tracker)+1}'
                    no_hit_counter += 1
                    new_id_tracker[new_id].append(seq)
                else:
                    # If already has Group naming, keep the original ID
                    new_id_tracker[seq.id].append(seq)
    else:
        # 5. Handle the case where there are no BLAST hits
        if not seq.id.startswith("Group"):
            new_id = f'Group_{no_hit_counter}_contig_{len(new_id_tracker)+1}'
            no_hit_counter += 1
            new_id_tracker[new_id].append(seq)
        else:
            # If already has Group naming, keep the original ID
            new_id_tracker[seq.id].append(seq)

# 6. Rename if multiple occurrences
final_seqs = []
for new_id, seqs in new_id_tracker.items():
    if len(seqs) > 1:
        # If multiple sequences have the same new ID, add suffixes
        for i, seq in enumerate(seqs, start=1):
            seq.id = f'{new_id}_{chr(64 + i)}'  # A, B, etc.
            seq.seq = seq.seq + "GATATC"
    else:
        # If only one sequence has this new ID, keep it as is
        seqs[0].id = new_id
        seqs[0].seq = seqs[0].seq + "GATATC"
    final_seqs.extend(seqs)

# 7. Write output
SeqIO.write(final_seqs, out_fasta, "fasta")
