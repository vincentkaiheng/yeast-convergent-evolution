#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet
import os

tsv_path = "./PFAM_upset_table.tsv"
figsize = (18, 9)
dpi = 600
output_pdf = "./UpSet_shared_novel_PFAMs.pdf"
output_png = "./UpSet_shared_novel_PFAMs.png"

class_order = [
    "Pucciniomycotina-yeast",
    "Ustilaginomycotina-yeast",
    "Agaricomycotina-yeast",
    "Taphrinomycotina-yeast",
    "Saccharomycotina-yeast",
    "Pezizomycotina-yeast",
]

df = pd.read_csv(tsv_path, sep="\t")

df = df[['PFAM'] + class_order]

indicator_df = df[class_order].astype(bool)

metadata_df = df[['PFAM']]

upset_data = metadata_df.copy()
upset_data.index = pd.MultiIndex.from_frame(indicator_df)

upset_data = upset_data[upset_data.index.to_frame().sum(axis=1) >= 2]

plt.rcParams.update({"font.size": 11})
fig = plt.figure(figsize=figsize, dpi=dpi)

upset = UpSet(
    upset_data,
    subset_size='count',
    show_counts=True,
    show_percentages=False,
    sort_by='cardinality',
    sort_categories_by='cardinality',
    orientation='horizontal',
    facecolor="#1f77b4",
    totals_plot_elements=6,
)

upset.plot(fig=fig)
fig.suptitle("Number of shared novel PFAMs", fontsize=18, y=0.98, weight='bold')

os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(output_pdf, dpi=dpi, bbox_inches='tight')
plt.savefig(output_png, dpi=dpi, bbox_inches='tight')
print(f"PDF: {output_pdf}")
print(f"PNG: {output_png}")
plt.show()
