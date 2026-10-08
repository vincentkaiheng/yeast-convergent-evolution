#!/usr/bin/env python3

import pandas as pd
import re

fungi_name_df = pd.read_csv('../../1_data_acquisition_and_functional_annotation/related_files/76_fungi_trait.txt', sep='\t')
fungi_name_dict = dict(zip(fungi_name_df['numbering'].astype(str), fungi_name_df['names']))

with open('related_files/gene_trees.tre', 'r') as file:
    lines = file.readlines()

new_lines = []
for line in lines:
    parts = line.strip().split('\t')
    tree = parts[1]

    tree = re.sub(r'[,\(](\d+):', lambda m: f'{m.group(0)[0]}{fungi_name_dict.get(m.group(1), m.group(1))}:', tree)

    new_lines.append(f"{parts[0]}\t{tree}\n")

with open('related_files/new_gene_trees.tre', 'w') as file:
    file.writelines(new_lines)
