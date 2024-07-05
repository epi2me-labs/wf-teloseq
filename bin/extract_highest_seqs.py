#!/usr/bin/env python
"""
This script extracts the sequence with the highest count.

Using a given input file and writes it to an output file.
"""

from collections import defaultdict
import os
import sys

input_filename = sys.argv[1]
min_seqs = round(float(sys.argv[2]))

sequences = defaultdict(str)

with open(input_filename) as f:
    header = None
    curr_seq = []

    for line in f:
        line = line.strip()

        if line.startswith('>'):
            if header:
                sequences[header] = ''.join(curr_seq)
            header = line
            curr_seq = []
        else:
            curr_seq.append(line)

    sequences[header] = ''.join(curr_seq)

if sequences:
    highest_header = max(sequences, key=lambda h: int(h.split(';')[1].split('=')[1]))
    highest_seqs = int(highest_header.split(';')[1].split('=')[1])

    if highest_seqs >= min_seqs:
        base_filename, ext = os.path.splitext(input_filename)
        output_filename = f'{base_filename}_highestseqs{ext}'
        name = base_filename.split('.')[0]

        with open(output_filename, 'w') as f:
            f.write('>' + name + '\n')
            f.write(sequences[highest_header] + '\n')
