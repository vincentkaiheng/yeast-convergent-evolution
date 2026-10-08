#!/usr/bin/env python3

import pandas as pd
import os

orthogroups_path = '../3_orthogroup_based_analysis_of_gene_family_evolution_across_independent_yeast_lineages/orthofinder_result/Orthogroups.GeneCount.tsv'
ml_label_path   = './related_files/ml_76_matrix_sens_label.txt'
output_path    = './related_files/ml_76_matrix_sens_labeled.tsv'

og = pd.read_csv(orthogroups_path, sep='\t', header=0)
og_t = og.set_index('Orthogroup').T.reset_index()
og_t = og_t.rename(columns={'index': 'species_full'})

ml = pd.read_csv(ml_label_path, sep='\t', header=0)
ml['species_name'] = ml['species_name'].str.strip()
if len(ml) != 76 or ml['species_name'].duplicated().any():
    raise ValueError('Expected 76 unique species in the sensitivity label table.')

matched_rows = []
unmatched = []

for _, row in ml.iterrows():
    sp_short = str(row['species_name']).strip()

    mask = og_t['species_full'].str.lower().str.contains(sp_short.lower(), regex=False, na=False)
    if mask.sum() == 1:
        counts = og_t.loc[mask, og_t.columns != 'species_full'].iloc[0]

        merged = pd.concat([row, counts], axis=0)
        matched_rows.append(merged)
    else:
        unmatched.append(sp_short)

if unmatched:
    print('Species without a unique matrix match:')
    for sp in unmatched:
        print(sp)
    raise ValueError("Species matching failed; no output was written.")

if matched_rows:
    result = pd.DataFrame(matched_rows)

    cols = list(ml.columns) + [c for c in result.columns if c not in ml.columns]
    result[cols].to_csv(output_path, sep='\t', index=False)
    print(f'Labeled matrix saved to: {os.path.abspath(output_path)}')
else:
    print('No species matched; no output was written.')
