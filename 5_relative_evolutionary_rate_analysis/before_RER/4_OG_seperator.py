#!/usr/bin/env python3

import os
import shutil

# Headers must use species_id|sequence_id.
pruned_alignment_dir = './trimmed_OGs/phylopypruner_output/output_alignments'
yeast_only_dir = './yeast_only_orthologs'
nonyeast_only_dir = './nonyeast_only_orthologs'
shared_ortholog_dir = './shared_orthologs'
foreground_list_file = './foreground_list.txt'

with open(foreground_list_file, 'r', encoding='utf-8') as f:
    foreground_species = set(line.strip() for line in f if line.strip())

os.makedirs(yeast_only_dir, exist_ok=True)
os.makedirs(nonyeast_only_dir, exist_ok=True)
os.makedirs(shared_ortholog_dir, exist_ok=True)

def parse_species_id_from_header(header_line: str) -> str:
    return header_line[1:].split('|', 1)[0].strip()

for filename in os.listdir(pruned_alignment_dir):
    if not filename.endswith('.fa'):
        continue

    file_path = os.path.join(pruned_alignment_dir, filename)

    contains_foreground = False
    contains_background = False

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if not line.startswith('>'):
                continue

            species_id = parse_species_id_from_header(line)
            if species_id in foreground_species:
                contains_foreground = True
            else:
                contains_background = True

            if contains_foreground and contains_background:
                break

    if contains_foreground and contains_background:
        shutil.copy(file_path, os.path.join(shared_ortholog_dir, filename))
    elif contains_foreground:
        shutil.copy(file_path, os.path.join(yeast_only_dir, filename))
    elif contains_background:
        shutil.copy(file_path, os.path.join(nonyeast_only_dir, filename))

print("Ortholog alignments classified.")
