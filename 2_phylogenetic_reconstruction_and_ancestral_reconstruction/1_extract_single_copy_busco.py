#!/bin/env python3
# -*- coding: UTF-8 -*-

import os
from Bio import SeqIO
from collections import defaultdict

base_dir = "./busco_result"
output_dir = "./combined_single_copy_sequences"
os.makedirs(output_dir, exist_ok=True)

# Keep the original sequence IDs for species labels in step 5.
gene_sequences = defaultdict(list)

for species_dir in os.listdir(base_dir):
    species_path = os.path.join(base_dir, species_dir)
    single_copy_dir = os.path.join(species_path, "run_fungi_odb12.2", "busco_sequences", "single_copy_busco_sequences")
    
    if not os.path.exists(single_copy_dir):
        print(f"Warning: {single_copy_dir} not found; skipping this species.")
        continue
    
    for file in os.listdir(single_copy_dir):
        if file.endswith('.faa'):
            gene_id = file.split('at')[0]
            file_path = os.path.join(single_copy_dir, file)
            with open(file_path, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    gene_sequences[gene_id].append(record)

for gene_id, sequences in gene_sequences.items():
    output_file = os.path.join(output_dir, f"{gene_id}.faa")
    with open(output_file, 'w') as out_handle:
        SeqIO.write(sequences, out_handle, 'fasta')

print(f"Single-copy BUSCO sequences saved to {output_dir}")
