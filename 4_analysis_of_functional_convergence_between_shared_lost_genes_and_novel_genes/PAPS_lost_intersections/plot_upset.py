#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet
import os

tsv_path = "./OGs_upset_outputs/upset_input_OG.tsv"
figsize = (18, 9)
dpi = 600
output_pdf = "./UpSet_shared_lost_OGs.pdf"
output_png = "./UpSet_shared_lost_OGs.png"

name_mapping = {
    "Saccharomycotina": "Node 6",
    "Pezizomycotina": "Node 5",
    "Taphrinomycotina": "Node 4",
    "Agaricomycotina": "Node 3",
    "Ustilaginomycotina": "Node 2",
    "Pucciniomycotina": "Node 1"
}
target_order = list(name_mapping.values())

df = pd.read_csv(tsv_path, sep="\t")
df = df.rename(columns=name_mapping)

indicator_df = df[target_order].astype(bool)
metadata_df = df[["OG", "HG"]]

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

upset.style_subsets(min_degree=3, facecolor="#b26d68")

ax_dict = upset.plot(fig=fig)

if "intersections" in ax_dict:
    ax_dict["intersections"].set_ylabel("Number of shared lost OGs",
                                        fontsize=10
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
