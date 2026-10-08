#!/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu


def plot_feature(df, feature, output_file):
    """
    Generate violin plot and perform Mann-Whitney U test
    """

    group_yeast = df.loc[df["trait"] == 1, feature].dropna().values
    group_fungi = df.loc[df["trait"] == 0, feature].dropna().values

    stat, p_value = mannwhitneyu(
        group_yeast,
        group_fungi,
        alternative="two-sided"
    )

    if p_value < 0.001:
        sig_label = "***"
    elif p_value < 0.01:
        sig_label = "**"
    elif p_value < 0.05:
        sig_label = "*"
    else:
        sig_label = "ns"

    colors = ["#ddb7aa", "#d0e1e3"]

    plt.figure(figsize=(6, 6))

    ax = sns.violinplot(
        x="group",
        y=feature,
        data=df,
        inner="box",
        palette=colors,
        cut=0
    )

    sns.stripplot(
        x="group",
        y=feature,
        data=df,
        color="k",
        size=4,
        jitter=True,
        alpha=0.6
    )

    ax.set_xlabel("")
    ax.set_ylabel(feature, fontsize=12)
    ax.tick_params(axis="both", labelsize=11)

    y_min = df[feature].min()
    y_max = df[feature].max()
    y_range = y_max - y_min

    if y_range == 0:
        y_range = 1

    ax.set_ylim(
        y_min - 0.05 * y_range,
        y_max + 0.25 * y_range
    )

    y_sig = y_max + 0.03 * y_range
    h = 0.02 * y_range
    x1, x2 = 0, 1

    ax.plot(
        [x1, x1, x2, x2],
        [y_sig - h, y_sig, y_sig, y_sig - h],
        lw=1.2,
        c="k"
    )

    ax.text(
        (x1 + x2) / 2,
        y_sig + h * 0.2,
        f"{sig_label}\n(p = {p_value:.3g})",
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold"
    )

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(
        f"{feature}: p = {p_value:.3g} ({sig_label}) -> {output_file}"
    )


def main():

    genome_stat_file = "./genome_stat_updated_with_gf_counts.tsv"
    trait_file = Path(__file__).resolve().parent / "related_files" / "76_fungi_trait.txt"

    print("Loading datasets...")

    genome_df = pd.read_csv(
        genome_stat_file,
        sep="\t"
    )

    trait_df = pd.read_csv(
        trait_file,
        sep="\t"
    )

    print("Merging datasets...")

    merged_df = pd.merge(
        genome_df,
        trait_df[["names", "trait"]],
        left_on="Genome",
        right_on="names",
        how="inner"
    )

    merged_df.drop(
        columns=["names"],
        inplace=True
    )

    merged_df["group"] = merged_df["trait"].map(
        {
            1: "General Yeast",
            0: "Multicellular Fungi"
        }
    )

    sns.set_theme(style="whitegrid")

    features = [
        "Size(Mb)",
        "GC(%)",
        "TE content(%)",
        "Gene(K)",
        "Gene family number"
    ]

    for feature in features:

        safe_name = (
            feature.replace("%", "pct")
                   .replace("(", "_")
                   .replace(")", "")
                   .replace(" ", "_")
        )

        output_file = f"{safe_name}_violin.png"

        plot_feature(
            merged_df,
            feature,
            output_file
        )

    print("All plots finished.")


if __name__ == "__main__":
    main()
