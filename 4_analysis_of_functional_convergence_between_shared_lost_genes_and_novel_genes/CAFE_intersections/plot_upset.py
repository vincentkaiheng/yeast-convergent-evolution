#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet
import os
import argparse

parser = argparse.ArgumentParser(description="Plot shared CAFE expansions or contractions.")
parser.add_argument("--change", choices=["expanded", "contracted"], default="expanded")
change = parser.parse_args().change

tsv_path = f"./OGs_{change}_upset_outputs/og_upset_input.tsv"
figsize = (18, 9)
dpi = 600
output_pdf = f"./UpSet_shared_{change}_OGs.pdf"
output_png = f"./UpSet_shared_{change}_OGs.png"

name_mapping = {
    "Saccharomycotina_yeast": "Node 6",
    "Pezizomycotina_yeast": "Node 5",
    "Taphrinomycotina_yeast": "Node 4",
    "Agaricomycotina_yeast": "Node 3",
    "Ustilaginomycotina_yeast": "Node 2",
    "Pucciniomycotina_yeast": "Node 1"
}
target_order = list(name_mapping.values())

df = pd.read_csv(tsv_path, sep="\t")
df = df.rename(columns=name_mapping)

indicator_df = df[target_order].astype(bool)
metadata_df = df[["OG"]]

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
    sort_categories_by=None,
    orientation='horizontal',
    facecolor="#7c7c7c",
    totals_plot_elements=6,
)

upset.style_subsets(min_degree=3, facecolor="#b89352")

ax_dict = upset.plot(fig=fig)

if "intersections" in ax_dict:
    ax_dict["intersections"].set_ylabel(f"Number of shared {change} OGs",
                                        fontsize=8
                                        )

if "totals" in ax_dict:
    for text in ax_dict["totals"].texts:
        text.set_visible(False)

for ax in fig.axes:
    if "Set sizes" in ax.get_xlabel() or "Set sizes" in ax.get_ylabel():
        ax.invert_xaxis()
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

os.makedirs(os.path.dirname(output_pdf), exist_ok=True)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig(output_pdf, dpi=dpi, bbox_inches='tight', transparent=True)
plt.savefig(output_png, dpi=dpi, bbox_inches='tight', transparent=True)
plt.show()
