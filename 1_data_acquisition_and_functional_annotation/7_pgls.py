#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compare genome features between groups using PGLS and plot distributions."""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from scipy.linalg import solve_triangular
from scipy.optimize import minimize_scalar
from scipy.stats import mannwhitneyu, t as student_t


BASE_DIR = Path(__file__).resolve().parent
GENOME_FILE = BASE_DIR / "related_files" / "genome_stat_updated_with_gf_counts.tsv"
TRAIT_FILE = BASE_DIR / "related_files" / "76_fungi_trait.txt"
TREE_FILE = BASE_DIR / "related_files" / "fungi_tree_named_rooted.nwk"
OUTPUT_DIR = BASE_DIR / "related_files"

FEATURES = [
    ("Size(Mb)", "Genome size (Mb)"),
    ("GC(%)", "GC content (%)"),
    ("TE content(%)", "TE content (%)"),
    ("Gene(K)", "Gene count (thousands)"),
    ("Gene family number", "Gene-family count"),
]

GROUP_LABELS = {
    1: "General Yeast",
    0: "Multicellular Fungi",
}
GROUP_ORDER = ["General Yeast", "Multicellular Fungi"]
PALETTE = {
    "General Yeast": "#ddb7aa",
    "Multicellular Fungi": "#d0e1e3",
}


def require_columns(df: pd.DataFrame, required: list[str], source: Path) -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {source}: {missing}")


def require_unique(df: pd.DataFrame, column: str, source: Path) -> None:
    duplicated = sorted(
        df.loc[df[column].duplicated(keep=False), column].astype(str).unique()
    )
    if duplicated:
        raise ValueError(
            f"Duplicate taxon names in {source}, column {column!r}: {duplicated}"
        )


def normalize_taxon_name(value: object) -> str:
    """Normalize whitespace and quotes in taxon names."""

    name = str(value).strip().strip("'\"")
    name = re.sub(r"\s+", "_", name)
    return name


def load_and_merge_data(
    genome_file: Path,
    trait_file: Path,
) -> pd.DataFrame:
    genome_df = pd.read_csv(genome_file, sep="\t")
    trait_df = pd.read_csv(trait_file, sep="\t")

    feature_names = [name for name, _ in FEATURES]
    require_columns(genome_df, ["Genome", *feature_names], genome_file)
    require_columns(trait_df, ["names", "trait"], trait_file)
    genome_df = genome_df.copy()
    trait_df = trait_df.copy()
    genome_df["Genome"] = genome_df["Genome"].map(normalize_taxon_name)
    trait_df["names"] = trait_df["names"].map(normalize_taxon_name)
    require_unique(genome_df, "Genome", genome_file)
    require_unique(trait_df, "names", trait_file)

    for feature in feature_names:
        genome_df[feature] = pd.to_numeric(genome_df[feature], errors="raise")
    trait_df["trait"] = pd.to_numeric(trait_df["trait"], errors="raise")

    if trait_df["trait"].isna().any():
        raise ValueError("The trait column contains missing values")
    observed_traits = set(trait_df["trait"].astype(int))
    if observed_traits != {0, 1}:
        raise ValueError(
            "The trait column must contain both 0 and 1 only; "
            f"observed values: {sorted(observed_traits)}"
        )
    if not np.all(trait_df["trait"] == trait_df["trait"].astype(int)):
        raise ValueError("The trait column contains non-integer values")
    trait_df["trait"] = trait_df["trait"].astype(int)

    genome_names = set(genome_df["Genome"])
    trait_names = set(trait_df["names"])
    if genome_names != trait_names:
        raise ValueError(
            "Genome/trait taxon mismatch. "
            f"Genome only: {sorted(genome_names - trait_names)}; "
            f"trait only: {sorted(trait_names - genome_names)}"
        )

    merged = genome_df.merge(
        trait_df[["names", "trait"]],
        left_on="Genome",
        right_on="names",
        how="inner",
        validate="one_to_one",
    ).drop(columns="names")
    merged["group"] = merged["trait"].map(GROUP_LABELS)
    return merged


def load_tree_and_match_taxa(
    tree_file: Path,
    data_names: list[str],
    output_dir: Path,
):
    tree = Phylo.read(str(tree_file), "newick")
    terminals = tree.get_terminals()
    tip_by_name = {}

    for tip in terminals:
        if tip.name is None:
            raise ValueError("The tree contains an unnamed terminal tip")
        normalized = normalize_taxon_name(tip.name)
        if normalized in tip_by_name:
            raise ValueError(
                "Tree-tip names collide after normalization: "
                f"{normalized!r}"
            )
        tip_by_name[normalized] = tip

    data_set = set(data_names)
    tree_set = set(tip_by_name)
    match_rows = []
    for name in sorted(data_set | tree_set):
        if name in data_set and name in tree_set:
            status = "matched"
        elif name in data_set:
            status = "data_only"
        else:
            status = "tree_only"
        match_rows.append({"taxon": name, "status": status})

    match_report = output_dir / "taxon_matching_report.tsv"
    pd.DataFrame(match_rows).to_csv(match_report, sep="\t", index=False)

    missing_from_tree = sorted(data_set - tree_set)
    if missing_from_tree:
        raise ValueError(
            "Some data taxa are absent from the tree. No taxa were silently "
            f"dropped. See {match_report}. Missing: {missing_from_tree}"
        )

    extra_tree_tips = sorted(tree_set - data_set)
    if extra_tree_tips:
        print(
            f"NOTE: {len(extra_tree_tips)} extra tree tips are not used; "
            f"see {match_report}"
        )

    bad_branches = []
    for clade in tree.find_clades(order="level"):
        if clade is tree.root:
            continue
        length = clade.branch_length
        if length is None or not math.isfinite(float(length)) or float(length) < 0:
            bad_branches.append((clade.name, length))
    if bad_branches:
        preview = bad_branches[:10]
        raise ValueError(
            "PGLS requires finite, non-negative branch lengths on every "
            f"non-root branch. Examples of invalid branches: {preview}"
        )

    tips_in_data_order = [tip_by_name[name] for name in data_names]
    return tree, tips_in_data_order, extra_tree_tips


def brownian_covariance(tree, ordered_tips: list) -> np.ndarray:
    """Return C_ij = root-to-MRCA path length in data-row order."""

    n_taxa = len(ordered_tips)
    covariance = np.zeros((n_taxa, n_taxa), dtype=float)

    for i, tip_i in enumerate(ordered_tips):
        depth = float(tree.distance(tree.root, tip_i))
        if not math.isfinite(depth) or depth <= 0:
            raise ValueError(
                f"Tip {tip_i.name!r} has invalid root-to-tip distance: {depth}"
            )
        covariance[i, i] = depth
        for j in range(i):
            mrca = tree.common_ancestor(tip_i, ordered_tips[j])
            shared_depth = float(tree.distance(tree.root, mrca))
            covariance[i, j] = shared_depth
            covariance[j, i] = shared_depth

    if not np.all(np.isfinite(covariance)):
        raise ValueError("The phylogenetic covariance matrix is not finite")

    # Rescale covariance to improve numerical conditioning.
    mean_tip_variance = float(np.mean(np.diag(covariance)))
    if mean_tip_variance <= 0:
        raise ValueError("Mean phylogenetic tip variance must be positive")
    return covariance / mean_tip_variance


def pagel_lambda_covariance(
    brownian_cov: np.ndarray,
    lambda_value: float,
    nugget: float = 1e-10,
) -> np.ndarray:
    diagonal = np.diag(np.diag(brownian_cov))
    transformed = diagonal + lambda_value * (brownian_cov - diagonal)
    return transformed + np.eye(transformed.shape[0]) * nugget


def whitened_fit(
    y: np.ndarray,
    design: np.ndarray,
    brownian_cov: np.ndarray,
    lambda_value: float,
):
    covariance = pagel_lambda_covariance(brownian_cov, lambda_value)
    chol = np.linalg.cholesky(covariance)
    y_white = solve_triangular(chol, y, lower=True, check_finite=False)
    x_white = solve_triangular(chol, design, lower=True, check_finite=False)
    information = x_white.T @ x_white
    beta = np.linalg.solve(information, x_white.T @ y_white)
    residual_white = y_white - x_white @ beta
    rss = float(residual_white @ residual_white)
    logdet_cov = 2.0 * float(np.log(np.diag(chol)).sum())
    return beta, rss, information, logdet_cov


def profile_log_likelihood(
    lambda_value: float,
    y: np.ndarray,
    design: np.ndarray,
    brownian_cov: np.ndarray,
) -> float:
    n = len(y)
    try:
        _, rss, _, logdet_cov = whitened_fit(
            y, design, brownian_cov, lambda_value
        )
    except np.linalg.LinAlgError:
        return -np.inf
    if not math.isfinite(rss) or rss <= 0:
        return -np.inf
    sigma2_ml = rss / n
    return -0.5 * (
        n * (math.log(2.0 * math.pi) + 1.0 + math.log(sigma2_ml))
        + logdet_cov
    )


def fit_pgls(
    y: np.ndarray,
    trait: np.ndarray,
    brownian_cov: np.ndarray,
) -> dict[str, float]:
    y = np.asarray(y, dtype=float)
    trait = np.asarray(trait, dtype=float)
    design = np.column_stack([np.ones(len(y)), trait])
    n, p = design.shape

    if n <= p:
        raise ValueError("Too few observations for PGLS")
    if set(np.unique(trait)) != {0.0, 1.0}:
        raise ValueError("Each feature must retain both trait groups")

    objective = lambda lam: -profile_log_likelihood(
        float(lam), y, design, brownian_cov
    )
    optimized = minimize_scalar(
        objective,
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1e-7, "maxiter": 500},
    )

    candidates = [
        (0.0, profile_log_likelihood(0.0, y, design, brownian_cov)),
        (1.0, profile_log_likelihood(1.0, y, design, brownian_cov)),
    ]
    if optimized.success and math.isfinite(float(optimized.fun)):
        optimized_lambda = float(np.clip(optimized.x, 0.0, 1.0))
        candidates.append(
            (
                optimized_lambda,
                profile_log_likelihood(
                    optimized_lambda, y, design, brownian_cov
                ),
            )
        )

    lambda_hat, log_likelihood = max(candidates, key=lambda item: item[1])
    if not math.isfinite(log_likelihood):
        raise RuntimeError("Pagel-lambda optimization failed for this feature")

    beta, rss, information, _ = whitened_fit(
        y, design, brownian_cov, lambda_hat
    )
    degrees_freedom = n - p
    sigma2_unbiased = rss / degrees_freedom
    beta_covariance = sigma2_unbiased * np.linalg.inv(information)
    standard_error = float(np.sqrt(beta_covariance[1, 1]))
    t_statistic = float(beta[1] / standard_error)
    p_value = float(
        2.0 * student_t.sf(abs(t_statistic), df=degrees_freedom)
    )
    critical_t = float(student_t.ppf(0.975, df=degrees_freedom))
    ci_low = float(beta[1] - critical_t * standard_error)
    ci_high = float(beta[1] + critical_t * standard_error)

    # Parameters counted for AIC: two betas, residual variance, and lambda.
    aic = float(2 * (p + 2) - 2 * log_likelihood)

    return {
        "beta_yeast_minus_multicellular": float(beta[1]),
        "standard_error": standard_error,
        "CI95_low": ci_low,
        "CI95_high": ci_high,
        "t": t_statistic,
        "df": int(degrees_freedom),
        "PGLS_p": p_value,
        "Pagel_lambda_ML": float(lambda_hat),
        "log_likelihood": float(log_likelihood),
        "AIC": aic,
    }


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    if np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("All P values must be finite and between 0 and 1")
    n = len(values)
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = ranked * n / np.arange(1, n + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0.0, 1.0)
    adjusted = np.empty(n, dtype=float)
    adjusted[order] = adjusted_ranked
    return adjusted


def significance_label(q_value: float) -> str:
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return "ns"


def format_probability(value: float) -> str:
    return f"{value:.3g}"


def safe_file_stem(feature: str) -> str:
    stem = feature.lower().replace("%", "pct")
    stem = re.sub(r"[^a-z0-9]+", "_", stem).strip("_")
    return stem


def plot_feature(
    data: pd.DataFrame,
    feature: str,
    y_label: str,
    p_value: float,
    q_value: float,
    output_dir: Path,
) -> None:
    plot_data = data.dropna(subset=[feature]).copy()
    plot_data["group"] = pd.Categorical(
        plot_data["group"], categories=GROUP_ORDER, ordered=True
    )

    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    sns.violinplot(
        data=plot_data,
        x="group",
        y=feature,
        hue="group",
        order=GROUP_ORDER,
        hue_order=GROUP_ORDER,
        palette=PALETTE,
        inner="box",
        cut=0,
        linewidth=1.0,
        legend=False,
        ax=ax,
    )
    sns.stripplot(
        data=plot_data,
        x="group",
        y=feature,
        order=GROUP_ORDER,
        color="black",
        size=4,
        jitter=0.18,
        alpha=0.6,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel(y_label, fontsize=12)
    ax.tick_params(axis="both", labelsize=11)

    y_min = float(plot_data[feature].min())
    y_max = float(plot_data[feature].max())
    y_range = y_max - y_min
    if y_range == 0:
        y_range = max(abs(y_max), 1.0)

    y_sig = y_max + 0.08 * y_range
    bracket_height = 0.025 * y_range
    ax.plot(
        [0, 0, 1, 1],
        [
            y_sig - bracket_height,
            y_sig,
            y_sig,
            y_sig - bracket_height,
        ],
        linewidth=1.2,
        color="black",
    )
    ax.text(
        0.5,
        y_sig + 0.01 * y_range,
        (
            f"PGLS p = {format_probability(p_value)}; "
            f"BH q = {format_probability(q_value)} "
            f"({significance_label(q_value)})"
        ),
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.set_ylim(y_min - 0.05 * y_range, y_max + 0.25 * y_range)
    sns.despine(ax=ax)
    fig.tight_layout()

    stem = safe_file_stem(feature)
    fig.savefig(output_dir / f"{stem}_violin_pgls.png", dpi=300)
    fig.savefig(output_dir / f"{stem}_violin_pgls.pdf")
    fig.savefig(output_dir / f"{stem}_violin_pgls.svg")
    plt.close(fig)


def main() -> None:
    genome_file = GENOME_FILE
    trait_file = TRAIT_FILE
    tree_file = TREE_FILE
    output_dir = OUTPUT_DIR

    for required_file in [genome_file, trait_file, tree_file]:
        if not required_file.is_file():
            raise FileNotFoundError(f"Required file not found: {required_file}")
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Genome metrics: {genome_file}")
    print(f"Traits: {trait_file}")
    print(f"Tree: {tree_file}")
    print(f"Output: {output_dir}")

    merged = load_and_merge_data(genome_file, trait_file)
    data_names = merged["Genome"].tolist()
    tree, ordered_tips, _ = load_tree_and_match_taxa(
        tree_file, data_names, output_dir
    )
    full_brownian_cov = brownian_covariance(tree, ordered_tips)

    results = []
    for feature, y_label in FEATURES:
        valid_mask = merged[feature].notna().to_numpy()
        subset = merged.loc[valid_mask].copy()
        covariance = full_brownian_cov[np.ix_(valid_mask, valid_mask)]

        yeast_values = subset.loc[subset["trait"] == 1, feature].to_numpy()
        multicellular_values = subset.loc[
            subset["trait"] == 0, feature
        ].to_numpy()
        mw_result = mannwhitneyu(
            yeast_values,
            multicellular_values,
            alternative="two-sided",
        )

        pgls_result = fit_pgls(
            subset[feature].to_numpy(dtype=float),
            subset["trait"].to_numpy(dtype=float),
            covariance,
        )
        results.append(
            {
                "metric": feature,
                "n_total": len(subset),
                "n_yeast": len(yeast_values),
                "n_multicellular": len(multicellular_values),
                **pgls_result,
                "Mann_Whitney_U_uncorrected": float(mw_result.statistic),
                "Mann_Whitney_p_uncorrected": float(mw_result.pvalue),
            }
        )

    results_df = pd.DataFrame(results)
    results_df["BH_FDR"] = benjamini_hochberg(results_df["PGLS_p"])
    results_df["FDR_significance"] = results_df["BH_FDR"].map(
        significance_label
    )

    results_file = output_dir / "genome_architecture_PGLS_results.tsv"
    results_df.to_csv(results_file, sep="\t", index=False, float_format="%.10g")
    merged.to_csv(
        output_dir / "genome_architecture_merged_data.tsv",
        sep="\t",
        index=False,
    )

    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })
    for feature, y_label in FEATURES:
        row = results_df.loc[results_df["metric"] == feature].iloc[0]
        plot_feature(
            merged,
            feature,
            y_label,
            float(row["PGLS_p"]),
            float(row["BH_FDR"]),
            output_dir,
        )

    display_columns = [
        "metric",
        "n_total",
        "beta_yeast_minus_multicellular",
        "CI95_low",
        "CI95_high",
        "Pagel_lambda_ML",
        "PGLS_p",
        "BH_FDR",
        "FDR_significance",
    ]
    print("\nPhylogenetically corrected results:")
    print(results_df[display_columns].to_string(index=False))
    print(f"\nSaved results and plots to: {output_dir}")


if __name__ == "__main__":
    main()
