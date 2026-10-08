#!/bin/env python3
# -*- coding: UTF-8 -*-

import os
import glob
import pandas as pd

df = pd.read_csv('genome_stat.tsv', sep='\t', dtype=str)
te_contents = []

for genome in df['Genome']:
    tbl_path = os.path.join('76_fna', f'{genome}.tbl')
    if not os.path.isfile(tbl_path):
        te_contents.append('NA')
        continue
    with open(tbl_path, 'r') as f:
        te_pct = 'NA'
        for line in f:
            line = line.strip()
            if line.startswith('bases masked'):
                start = line.find('(')
                end = line.find(')')
                if start != -1 and end != -1:
                    content = line[start+1:end].strip()
                    te_pct = content.replace('%', '').strip()
                break
        te_contents.append(te_pct)

df['TE content(%)'] = te_contents
df.to_csv('genome_stat.tsv', sep='\t', index=False)

print('Added TE content (%) to genome_stat.tsv.')
