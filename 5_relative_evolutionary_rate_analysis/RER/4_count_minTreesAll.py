#!/usr/bin/env python3

import os

input_dir = './shared_orthologs_sp_only'

expected_species_count = 76

qualified_files = []

for filename in os.listdir(input_dir):
    if not filename.endswith('.fa'):
        continue
    
    file_path = os.path.join(input_dir, filename)
    

    with open(file_path, 'r', encoding='utf-8') as f:
        seq_count = sum(1 for line in f if line.startswith('>'))
    
    if seq_count == expected_species_count:
        qualified_files.append(filename)

print(f" A total of {len(qualified_files)} fasta files containing all {expected_species_count} species")
