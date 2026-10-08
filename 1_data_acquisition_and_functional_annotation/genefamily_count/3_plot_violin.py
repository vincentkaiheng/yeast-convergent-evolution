#!/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import numpy as np

input_file = "./Orthogroups.presence_counts.tsv"
output_file = "./gene_family_violin_plot.png"

df = pd.read_csv(input_file, sep="\t")

df['group'] = df['trait'].map({1: 'General Yeast', 0: 'Multicellular Fungi'})

uni_counts = df.loc[df['trait'] == 1, 'present_OG_count'].values
multi_counts = df.loc[df['trait'] == 0, 'present_OG_count'].values

# Mann-Whitney U test (one-sided: uni < multi)
stat, p_value = mannwhitneyu(uni_counts, multi_counts, alternative='less')

if p_value < 0.001:
    sig_label = '***'
elif p_value < 0.01:
    sig_label = '**'
elif p_value < 0.05:
    sig_label = '*'
else:
    sig_label = 'ns'

sns.set(style="whitegrid")
plt.figure(figsize=(6, 6))

ax = sns.violinplot(
    x='group',
    y='present_OG_count',
    data=df,
    inner='box',
    palette=["#ddb7aa", "#d0e1e3"],
    cut=0
)

sns.stripplot(
    x='group',
    y='present_OG_count',
    data=df,
    color='k',
    size=4,
    jitter=True,
    alpha=0.6
)

ax.set_xlabel("")
ax.set_ylabel("Number of Gene Families", fontsize=12)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=12)

y_ticks = [3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500]
ax.set_yticks(y_ticks)
ax.set_ylim(min(y_ticks), max(y_ticks))

y_max_data = df['present_OG_count'].max()
y_sig = min(max(y_ticks), y_max_data + (max(y_ticks)-min(y_ticks))*0.05)  # 5% above max or top tick

x1, x2 = 0, 1  # positions of the two groups
ax.plot([x1, x1, x2, x2], [y_sig-50, y_sig, y_sig, y_sig-50], lw=1.5, c='k')
ax.text((x1+x2)/2, y_sig + 30, sig_label, ha='center', va='bottom', fontsize=14)

plt.tight_layout()
plt.savefig(output_file, dpi=300)
plt.close()

print(f"Violin plot saved as: {output_file}, Mann-Whitney p-value = {p_value:.3g}, significance = {sig_label}")
