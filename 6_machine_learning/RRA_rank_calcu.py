#!/usr/bin/env python3

from pathlib import Path
from math import comb

import numpy as np
import pandas as pd

MAIN_GAIN_PATH = Path(
    "results/ml_output_main_xgb/"
    "main_binary_xgbrf_gain_importance.tsv"
)

MAIN_SHAP_PATH = Path(
    "results/ml_output_main_xgb/"
    "main_binary_xgbrf_shap_importance.tsv"
)

SENS_A_GAIN_PATH = Path(
    "results/ml_output_sens_xgb/sens_A/"
    "xgbrf_gain_importance.tsv"
)

SENS_A_SHAP_PATH = Path(
    "results/ml_output_sens_xgb/sens_A/"
    "xgbrf_shap_importance.tsv"
)

OUTPUT_DIR = Path(
    "results/ml_output_rra"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TOP_K = 30
MIN_NONZERO_FRACTION = 0.90
REQUIRE_POSITIVE_MEDIAN = True

REQUIRE_TOP_K_IN_ALL_LISTS = True

def read_importance_table(
    file_path,
    prefix,
    score_column,
    median_column
):
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    df = pd.read_csv(file_path, sep="\t")

    required_columns = {
        "orthogroup",
        score_column,
        median_column,
        "nonzero_fold_fraction"
    }

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(
            f"{file_path.name} is missing columns: {sorted(missing_columns)}"
        )

    if df["orthogroup"].duplicated().any():
        duplicated = df.loc[
            df["orthogroup"].duplicated(),
            "orthogroup"
        ].tolist()

        raise ValueError(
            f"{file_path.name} contains duplicate OGs, including: "
            f"{duplicated[:10]}"
        )

    numeric_columns = [
        score_column,
        median_column,
        "nonzero_fold_fraction"
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    if df[numeric_columns].isna().any().any():
        bad_columns = df[numeric_columns].columns[
            df[numeric_columns].isna().any()
        ].tolist()

        raise ValueError(
            f"{file_path.name} has missing or nonnumeric values in: "
            f"{bad_columns}"
        )

    df = df.sort_values(
        score_column,
        ascending=False,
        kind="mergesort"
    ).reset_index(drop=True)

    df[f"{prefix}_rank"] = (
        df[score_column]
        .rank(method="min", ascending=False)
        .astype(int)
    )

    rename_dict = {
        score_column: f"{prefix}_score",
        median_column: f"{prefix}_median",
        "nonzero_fold_fraction": f"{prefix}_nonzero_fraction"
    }

    std_candidates = [
        c for c in df.columns
        if c.startswith("std_")
    ]

    if std_candidates:
        rename_dict[std_candidates[0]] = f"{prefix}_std"

    df = df.rename(columns=rename_dict)
    df = df.set_index("orthogroup")

    keep_columns = [
        f"{prefix}_score",
        f"{prefix}_median",
        f"{prefix}_nonzero_fraction",
        f"{prefix}_rank"
    ]

    if f"{prefix}_std" in df.columns:
        keep_columns.insert(1, f"{prefix}_std")

    return df[keep_columns]

main_gain = read_importance_table(
    MAIN_GAIN_PATH,
    prefix="main_gain",
    score_column="mean_gain",
    median_column="median_gain"
)

main_shap = read_importance_table(
    MAIN_SHAP_PATH,
    prefix="main_shap",
    score_column="mean_abs_shap",
    median_column="median_abs_shap"
)

scenario_a_gain = read_importance_table(
    SENS_A_GAIN_PATH,
    prefix="scenario_a_gain",
    score_column="mean_gain",
    median_column="median_gain"
)

scenario_a_shap = read_importance_table(
    SENS_A_SHAP_PATH,
    prefix="scenario_a_shap",
    score_column="mean_abs_shap",
    median_column="median_abs_shap"
)

all_tables = {
    "main_gain": main_gain,
    "main_shap": main_shap,
    "scenario_a_gain": scenario_a_gain,
    "scenario_a_shap": scenario_a_shap
}

common_orthogroups = sorted(
    set.intersection(
        *(set(df.index) for df in all_tables.values())
    )
)

if not common_orthogroups:
    raise RuntimeError("No OGs are shared by the four importance tables.")

print("OG counts by table:")
for name, df in all_tables.items():
    print(f"  {name}: {len(df):,}")

print(f"\nOGs shared by all four tables: {len(common_orthogroups):,}")

combined = pd.concat(
    [
        main_gain.loc[common_orthogroups],
        main_shap.loc[common_orthogroups],
        scenario_a_gain.loc[common_orthogroups],
        scenario_a_shap.loc[common_orthogroups]
    ],
    axis=1
)

score_columns = {
    "main_gain": "main_gain_score",
    "main_shap": "main_shap_score",
    "scenario_a_gain": "scenario_a_gain_score",
    "scenario_a_shap": "scenario_a_shap_score"
}

common_rank_columns = []

for prefix, score_column in score_columns.items():
    rank_column = f"{prefix}_common_rank"
    percentile_column = f"{prefix}_percentile"

    combined[rank_column] = (
        combined[score_column]
        .rank(method="min", ascending=False)
        .astype(int)
    )

    combined[percentile_column] = (
        combined[rank_column] / len(combined)
    )

    common_rank_columns.append(rank_column)

def beta_order_statistic_cdf(p, k, m):
    """Evaluate Beta(k, m-k+1) using the equivalent binomial tail."""
    p = float(np.clip(p, 0.0, 1.0))

    return sum(
        comb(m, j) *
        (p ** j) *
        ((1.0 - p) ** (m - j))
        for j in range(k, m + 1)
    )

def conservative_rra_pvalue(normalized_ranks):
    """Apply Bonferroni correction to the minimum order-statistic probability."""
    ranks = np.sort(
        np.asarray(normalized_ranks, dtype=float)
    )

    m = len(ranks)

    beta_probabilities = [
        beta_order_statistic_cdf(
            p=ranks[k - 1],
            k=k,
            m=m
        )
        for k in range(1, m + 1)
    ]

    rho = min(beta_probabilities)

    return min(1.0, m * rho)

percentile_columns = [
    "main_gain_percentile",
    "main_shap_percentile",
    "scenario_a_gain_percentile",
    "scenario_a_shap_percentile"
]

combined["rra_p_conservative"] = combined[
    percentile_columns
].apply(
    lambda row: conservative_rra_pvalue(row.values),
    axis=1
)

def benjamini_hochberg(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)

    order = np.argsort(pvalues)
    sorted_p = pvalues[order]

    adjusted = sorted_p * n / np.arange(1, n + 1)

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.clip(adjusted, 0.0, 1.0)

    result = np.empty(n, dtype=float)
    result[order] = adjusted

    return result

combined["rra_q_bh"] = benjamini_hochberg(
    combined["rra_p_conservative"].values
)

combined["minus_log10_rra_q"] = -np.log10(
    np.clip(combined["rra_q_bh"], 1e-300, 1.0)
)

nonzero_columns = [
    "main_gain_nonzero_fraction",
    "main_shap_nonzero_fraction",
    "scenario_a_gain_nonzero_fraction",
    "scenario_a_shap_nonzero_fraction"
]

median_columns = [
    "main_gain_median",
    "main_shap_median",
    "scenario_a_gain_median",
    "scenario_a_shap_median"
]

combined["minimum_nonzero_fraction"] = combined[
    nonzero_columns
].min(axis=1)

combined["maximum_rank"] = combined[
    common_rank_columns
].max(axis=1)

combined["geometric_mean_rank"] = np.exp(
    np.log(
        combined[common_rank_columns].astype(float)
    ).mean(axis=1)
)

combined["passes_nonzero_threshold"] = (
    combined["minimum_nonzero_fraction"]
    >= MIN_NONZERO_FRACTION
)

combined["passes_positive_median"] = (
    combined[median_columns] > 0
).all(axis=1)

combined["passes_top_k"] = (
    combined[common_rank_columns] <= TOP_K
).all(axis=1)

selection_mask = combined[
    "passes_nonzero_threshold"
].copy()

if REQUIRE_POSITIVE_MEDIAN:
    selection_mask &= combined[
        "passes_positive_median"
    ]

if REQUIRE_TOP_K_IN_ALL_LISTS:
    selection_mask &= combined[
        "passes_top_k"
    ]

core_candidates = (
    combined.loc[selection_mask]
    .sort_values(
        [
            "rra_p_conservative",
            "maximum_rank",
            "minimum_nonzero_fraction",
            "geometric_mean_rank"
        ],
        ascending=[True, True, False, True]
    )
    .copy()
)

core_candidates.insert(
    0,
    "core_rank",
    np.arange(1, len(core_candidates) + 1)
)

all_output_path = (
    OUTPUT_DIR /
    "all_common_orthogroups_rra.tsv"
)

core_output_path = (
    OUTPUT_DIR /
    "core_candidates_rra.tsv"
)

all_results = (
    combined
    .sort_values(
        [
            "rra_p_conservative",
            "maximum_rank",
            "minimum_nonzero_fraction"
        ],
        ascending=[True, True, False]
    )
)

all_results.to_csv(
    all_output_path,
    sep="\t",
    index=True,
    index_label="orthogroup"
)

core_candidates.to_csv(
    core_output_path,
    sep="\t",
    index=True,
    index_label="orthogroup"
)

print("\nSelection settings:")
print(f"  TOP_K = {TOP_K}")
print(
    f"  MIN_NONZERO_FRACTION = "
    f"{MIN_NONZERO_FRACTION}"
)
print(
    f"  REQUIRE_POSITIVE_MEDIAN = "
    f"{REQUIRE_POSITIVE_MEDIAN}"
)
print(
    f"  REQUIRE_TOP_K_IN_ALL_LISTS = "
    f"{REQUIRE_TOP_K_IN_ALL_LISTS}"
)

print(
    f"\nCore candidates: "
    f"{len(core_candidates)}"
)

display_columns = [
    "core_rank",
    "main_gain_common_rank",
    "main_shap_common_rank",
    "scenario_a_gain_common_rank",
    "scenario_a_shap_common_rank",
    "minimum_nonzero_fraction",
    "maximum_rank",
    "geometric_mean_rank",
    "rra_p_conservative",
    "rra_q_bh"
]

if len(core_candidates) > 0:
    print("\nCore candidate ranking:")
    print(
        core_candidates[
            display_columns
        ].to_string()
    )
else:
    print(
        "\nNo OGs meet all selection criteria."
        ""
        ""
    )

print("\nOutput files:")
print(all_output_path)
print(core_output_path)
