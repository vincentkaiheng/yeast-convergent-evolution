#!/bin/env python3
# -*- coding: UTF-8 -*-

from pathlib import Path

import pandas as pd

def process_table(orthogroups_file, fungi_names_file, output_file):
    fungi_df = pd.read_csv(fungi_names_file, sep='\t', usecols=['names', 'file_names'])
    name_mapping = dict(zip(fungi_df['file_names'], fungi_df['names']))
    
    ortho_df = pd.read_csv(orthogroups_file, sep='\t')
    
    ortho_df.insert(0, 'desc', ortho_df['Orthogroup'])
    
    ortho_df.rename(columns=name_mapping, inplace=True)
    
    ortho_df = ortho_df.drop(columns=['Total'])
    
    ortho_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Table saved to: {output_file}")

def filter_gene_families(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    
    filtered_df = df[(df.iloc[:, 2:] < 100).all(axis=1)]
    
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Filtered families saved to: {output_file}")

base_dir = Path(__file__).resolve().parent
related_dir = base_dir / 'cafe_related_files'
trait_file = (base_dir.parents[1] / '1_data_acquisition_and_functional_annotation'
              / 'related_files' / '76_fungi_trait.txt')
process_table(related_dir / 'Orthogroups.GeneCount.tsv', trait_file,
              related_dir / 'gene_families.txt')
filter_gene_families(related_dir / 'gene_families.txt',
                     related_dir / 'gene_families_filtered.txt')
