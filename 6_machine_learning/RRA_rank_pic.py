#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

OUTPUT_DIR = Path(
    "results/ml_output_rra"
)

CORE_RESULT_PATH = (
    OUTPUT_DIR /
    "core_candidates_rra.tsv"
)

PNG_PATH = (
    OUTPUT_DIR /
    "core_candidates_rra_ranking.png"
)

PDF_PATH = (
    OUTPUT_DIR /
    "core_candidates_rra_ranking.pdf"
)

SVG_PATH = (
    OUTPUT_DIR /
    "core_candidates_rra_ranking.svg"
)

TARGET_OG = ""
MAX_CANDIDATES_TO_PLOT = 24

core = pd.read_csv(
    CORE_RESULT_PATH,
    sep="\t"
)

if core.empty:
    raise RuntimeError(
        "No core candidates are available to plot."
    )

if "core_rank" not in core.columns:
    raise ValueError(
        "Input file is missing core_rank."
    )

core = (
    core.sort_values("core_rank")
    .head(MAX_CANDIDATES_TO_PLOT)
    .copy()
)

core["plot_score"] = -np.log10(
    np.clip(
        core["rra_q_bh"].astype(float),
        1e-300,
        1.0
    )
)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none"
})

DEFAULT_COLOR = "#2878B5"
TARGET_COLOR = "#D55E00"
GRID_COLOR = "#D9E1E8"
TEXT_COLOR = "#263238"

n_candidates = len(core)

figure_height = max(
    4.5,
    0.48 * n_candidates + 1.8
)

fig = plt.figure(
    figsize=(12.5, figure_height),
    constrained_layout=True
)

grid = fig.add_gridspec(
    nrows=1,
    ncols=2,
    width_ratios=[1.05, 1.25]
)

ax_rank = fig.add_subplot(grid[0, 0])
ax_heat = fig.add_subplot(grid[0, 1])

y_positions = np.arange(n_candidates)

bar_colors = [
    TARGET_COLOR
    if og == TARGET_OG
    else DEFAULT_COLOR
    for og in core["orthogroup"]
]

bars = ax_rank.barh(
    y_positions,
    core["plot_score"],
    color=bar_colors,
    edgecolor="white",
    linewidth=0.7,
    height=0.72
)

ax_rank.set_yticks(y_positions)
ax_rank.set_yticklabels(
    core["orthogroup"]
)

ax_rank.invert_yaxis()

ax_rank.set_xlabel(
    r"$-\log_{10}$(BH-adjusted RRA P value)"
)

ax_rank.set_title(
    "A  Robust rank aggregation",
    loc="left",
    fontweight="bold",
    color=TEXT_COLOR
)

ax_rank.grid(
    axis="x",
    linestyle="--",
    linewidth=0.6,
    alpha=0.7,
    color=GRID_COLOR
)

ax_rank.set_axisbelow(True)

for spine in ["top", "right", "left"]:
    ax_rank.spines[spine].set_visible(False)

ax_rank.spines["bottom"].set_color("#80919D")

x_offset = max(core["plot_score"].max() * 0.015, 0.02)

for bar, rank in zip(
    bars,
    core["core_rank"]
):
    ax_rank.text(
        bar.get_width() + x_offset,
        bar.get_y() + bar.get_height() / 2,
        f"#{int(rank)}",
        va="center",
        ha="left",
        fontsize=8.5,
        color=TEXT_COLOR
    )

for tick_label in ax_rank.get_yticklabels():
    if tick_label.get_text() == TARGET_OG:
        tick_label.set_color(TARGET_COLOR)
        tick_label.set_fontweight("bold")

rank_columns = [
    "main_gain_common_rank",
    "main_shap_common_rank",
    "scenario_a_gain_common_rank",
    "scenario_a_shap_common_rank"
]

column_labels = [
    "Main\nGain",
    "Main\nSHAP",
    "Scenario A\nGain",
    "Scenario A\nSHAP"
]

rank_matrix = (
    core
    .set_index("orthogroup")
    .loc[core["orthogroup"], rank_columns]
    .astype(float)
)

heat_values = np.log2(
    rank_matrix.values + 1.0
)

vmin = heat_values.min()
vmax = heat_values.max()

image = ax_heat.imshow(
    heat_values,
    cmap="YlGnBu_r",
    aspect="auto",
    interpolation="nearest",
    norm=Normalize(vmin=vmin, vmax=vmax)
)

ax_heat.set_xticks(
    np.arange(len(column_labels))
)
ax_heat.set_xticklabels(
    column_labels,
    rotation=0
)

ax_heat.set_yticks(
    np.arange(n_candidates)
)
ax_heat.set_yticklabels(
    core["orthogroup"]
)

ax_heat.set_title(
    "B  Rank consistency across analyses",
    loc="left",
    fontweight="bold",
    color=TEXT_COLOR
)

for row_index in range(rank_matrix.shape[0]):
    for column_index in range(rank_matrix.shape[1]):
        rank_value = int(
            rank_matrix.iloc[
                row_index,
                column_index
            ]
        )

        background_value = heat_values[
            row_index,
            column_index
        ]

        midpoint = (vmin + vmax) / 2

        text_color = (
            "white"
            if background_value < midpoint
            else "#263238"
        )

        ax_heat.text(
            column_index,
            row_index,
            str(rank_value),
            ha="center",
            va="center",
            fontsize=9,
            color=text_color,
            fontweight=(
                "bold"
                if core.iloc[row_index]["orthogroup"]
                == TARGET_OG
                else "normal"
            )
        )

if TARGET_OG in core["orthogroup"].values:
    target_row = np.where(
        core["orthogroup"].values == TARGET_OG
    )[0][0]

    rectangle = plt.Rectangle(
        (-0.5, target_row - 0.5),
        width=len(rank_columns),
        height=1,
        fill=False,
        edgecolor=TARGET_COLOR,
        linewidth=2.0
    )

    ax_heat.add_patch(rectangle)

for tick_label in ax_heat.get_yticklabels():
    if tick_label.get_text() == TARGET_OG:
        tick_label.set_color(TARGET_COLOR)
        tick_label.set_fontweight("bold")

ax_heat.set_xticks(
    np.arange(-0.5, len(rank_columns), 1),
    minor=True
)

ax_heat.set_yticks(
    np.arange(-0.5, n_candidates, 1),
    minor=True
)

ax_heat.grid(
    which="minor",
    color="white",
    linestyle="-",
    linewidth=1.2
)

ax_heat.tick_params(
    which="minor",
    bottom=False,
    left=False
)

for spine in ax_heat.spines.values():
    spine.set_visible(False)

colorbar = fig.colorbar(
    ScalarMappable(
        norm=Normalize(vmin=vmin, vmax=vmax),
        cmap="YlGnBu_r"
    ),
    ax=ax_heat,
    fraction=0.045,
    pad=0.03
)

colorbar.set_label(
    r"$\log_2$(rank + 1); lower is better",
    rotation=90
)

fig.suptitle(
    "Robustly prioritized orthogroups associated with yeast-form capability",
    fontsize=14,
    fontweight="bold",
    color=TEXT_COLOR
)

fig.savefig(
    PNG_PATH,
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)

fig.savefig(
    PDF_PATH,
    bbox_inches="tight",
    facecolor="white"
)

fig.savefig(
    SVG_PATH,
    bbox_inches="tight",
    facecolor="white"
)

plt.show()
plt.close(fig)

print("Figures saved to:")
print(PNG_PATH)
print(PDF_PATH)
print(SVG_PATH)
