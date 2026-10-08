#!/bin/env python3
# -*- coding: UTF-8 -*-

from pathlib import Path

import pandas as pd

# Example: ../../3_orthogroup_based_analysis_of_gene_family_evolution_across_independent_yeast_lineages/orthofinder_result/Orthogroups.GeneCount.tsv
orthogroup_file = "/home/yankh/76_converge/ml/Orthogroups.GeneCount.tsv"
mapping_file = Path(__file__).resolve().parents[1] / "related_files" / "76_fungi_trait.txt"
output_file = "./Orthogroups.GeneCount_renamed.tsv"

mapping_df = pd.read_csv(mapping_file, sep="\t", header=0)
mapping_dict = dict(zip(mapping_df["file_names"], mapping_df["names"]))

df = pd.read_csv(orthogroup_file, sep="\t", header=0)

# Remove last column if its header is "Total"
if df.columns[-1] == "Total":
    df = df.iloc[:, :-1]

new_columns = [mapping_dict.get(col, col) for col in df.columns]
df.columns = new_columns

df.to_csv(output_file, sep="\t", index=False)
