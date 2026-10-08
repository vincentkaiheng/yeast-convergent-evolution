#!/bin/env python3
# -*- coding: UTF-8 -*-

import os
import pandas as pd
from Bio import SeqIO

def load_mapping(path):
    df = pd.read_csv(path, sep='\t', header=None, names=["accession_number", "names", "numbering", "binary", "file_names"], dtype=str)
    df.set_index('accession_number', inplace=True)
    return df

genome_stat_path = './genome_stat.tsv'
gca_map_path    = '/home/yankh/yeast_converge/genomes/genomic_data/GCA.txt'
extra_map_path  = '/home/yankh/yeast_converge/genomes/genomic_data/extra.txt'
# Example: ./76_faa/
faa_dir         = '/home/yankh/76_converge/76_faa'

genome_df = pd.read_csv(genome_stat_path, sep='\t', dtype=str)
gca_map   = load_mapping(gca_map_path)
extra_map = load_mapping(extra_map_path)

new_names   = []
gene_counts = []

for orig in genome_df['Genome']:
    if orig.startswith('GCA_'):
        # Extract accession part e.g. GCA_000002945.2
        import re
        m = re.match(r'^(GCA_\d+\.\d+)', orig)
        if not m:
            raise ValueError(f"Could not parse accession from {orig}")
        acc = m.group(1)
        if acc not in gca_map.index:
            raise KeyError(f"Accession {acc} not found in GCA mapping file")
        mapped = gca_map.at[acc, 'names']
        file_base = gca_map.at[acc, 'file_names']
    else:
        # Remove .fna suffix
        key = os.path.splitext(orig)[0]
        if key not in extra_map.index:
            raise KeyError(f"Key {key} not found in extra mapping file")
        mapped = extra_map.at[key, 'names']
        file_base = extra_map.at[key, 'file_names']

    new_names.append(mapped)

    fasta_file = os.path.join(faa_dir, file_base + '.faa')
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(f"Fasta file not found: {fasta_file}")

    # Count sequences
    count = sum(1 for _ in SeqIO.parse(fasta_file, 'fasta'))
    # Convert to K unit
    gene_counts.append(count / 1000.0)

genome_df['Genome'] = new_names
genome_df['Gene(K)'] = [f"{x:.3f}" for x in gene_counts]

output_path = './genome_stat_updated.tsv'
genome_df.to_csv(output_path, sep='\t', index=False)
print(f"Updated genome_stat saved to: {output_path}")
