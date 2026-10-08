#!/bin/env python3
# -*- coding: utf-8 -*-

"""Compare orthogroup counts and identify families absent from all yeasts."""

import sys
import os
from pathlib import Path

import pandas as pd
import numpy as np

ORTHO_FILE = "./Orthogroups.GeneCount_renamed.tsv"
MAPPING_FILE = Path(__file__).resolve().parents[1] / "related_files" / "76_fungi_trait.txt"
OUT_COUNTS = "./Orthogroups.presence_counts.tsv"
OUT_MISSING = "./Orthogroups.missing_in_unicellular.tsv"
OUT_SUMMARY = "./Orthogroups.group_summary.txt"

P_VALUE_THRESHOLD = 0.05
# For identifying "missing" OGs: fraction of multicellular species that must have the OG
MULTIC_FRAC_THRESHOLD = 0.5

use_scipy = True
try:
    from scipy.stats import mannwhitneyu
except Exception:
    use_scipy = False

def main():
    for p in (ORTHO_FILE, MAPPING_FILE):
        if not os.path.exists(p):
            print(f"ERROR: required file not found: {p}", file=sys.stderr)
            sys.exit(1)

    try:
        mapping_df = pd.read_csv(MAPPING_FILE, sep="\t", header=0, dtype=str)
    except Exception as e:
        print(f"ERROR: failed to read mapping file: {e}", file=sys.stderr)
        sys.exit(1)

    if 'names' not in mapping_df.columns or 'trait' not in mapping_df.columns:
        print("ERROR: mapping file must contain columns named 'names' and 'trait' (trait: 1=unicellular,0=multicellular).", file=sys.stderr)
        sys.exit(1)

    mapping_df['trait'] = mapping_df['trait'].astype(str).str.strip()
    mapping_df = mapping_df[mapping_df['trait'].isin(['0','1'])]
    mapping_df['trait'] = mapping_df['trait'].astype(int)

    try:
        ortho_df = pd.read_csv(ORTHO_FILE, sep="\t", header=0, dtype=str)
    except Exception as e:
        print(f"ERROR: failed to read orthogroup file: {e}", file=sys.stderr)
        sys.exit(1)

    ortho_id_col = ortho_df.columns[0]
    species_cols = [c for c in ortho_df.columns if c != ortho_id_col]

    species_in_mapping = [s for s in species_cols if s in set(mapping_df['names'])]
    if len(species_in_mapping) == 0:
        print("ERROR: no species columns in orthogroup file match mapping 'names'.", file=sys.stderr)
        sys.exit(1)

    missing_species = set(mapping_df['names']) - set(species_in_mapping)
    if len(missing_species) > 0:
        print(f"WARNING: {len(missing_species)} species from mapping are not present in orthogroup table columns. They will be ignored.", file=sys.stderr)

    ortho_sub = ortho_df[[ortho_id_col] + species_in_mapping].copy()
    # Convert species columns to numeric counts (non-numeric -> 0)
    for c in species_in_mapping:
        ortho_sub[c] = pd.to_numeric(ortho_sub[c], errors='coerce').fillna(0).astype(float)

    presence = (ortho_sub[species_in_mapping] > 0).astype(int)
    presence.insert(0, ortho_id_col, ortho_sub[ortho_id_col])

    per_species_counts = presence[species_in_mapping].sum(axis=0).reset_index()
    per_species_counts.columns = ['species', 'present_OG_count']
    trait_map = mapping_df.set_index('names')['trait'].to_dict()
    per_species_counts['trait'] = per_species_counts['species'].map(trait_map).astype(int)

    uni_counts = per_species_counts.loc[per_species_counts['trait'] == 1, 'present_OG_count'].values
    multi_counts = per_species_counts.loc[per_species_counts['trait'] == 0, 'present_OG_count'].values

    def describe(arr):
        if len(arr) == 0:
            return {'n':0, 'mean':np.nan, 'median':np.nan, 'std':np.nan}
        return {'n': int(len(arr)), 'mean': float(np.mean(arr)), 'median': float(np.median(arr)), 'std': float(np.std(arr, ddof=1)) if len(arr)>1 else 0.0}

    uni_desc = describe(uni_counts)
    multi_desc = describe(multi_counts)

    # Statistical test: Mann-Whitney U (one-sided: alternative='less' tests uni < multi)
    p_value = None
    test_used = None
    test_note = ""
    if len(uni_counts) >= 1 and len(multi_counts) >= 1:
        if use_scipy:
            try:
                stat, p_value = mannwhitneyu(uni_counts, multi_counts, alternative='less')
                test_used = "Mann-Whitney U (one-sided, alternative='less')"
            except TypeError:
                # older scipy may not support 'alternative'; fall back to two-sided and halve p-value if appropriate
                stat, p_two = mannwhitneyu(uni_counts, multi_counts)
                # compute one-sided p-value approximation
                p_value = p_two / 2.0
                test_used = "Mann-Whitney U (two-sided result converted to one-sided)"
                test_note = "Note: SciPy version did not support 'alternative' parameter; two-sided p-value halved for one-sided inference."
        else:
            print("ERROR: scipy not available; please install scipy to perform the statistical test (pip install scipy).", file=sys.stderr)
            sys.exit(1)
    else:
        print("ERROR: insufficient species in one of the groups to run the statistical test.", file=sys.stderr)
        sys.exit(1)

    significant = (p_value is not None) and (p_value < P_VALUE_THRESHOLD)

    missing_ogs_df = None
    if significant:
        uni_species = per_species_counts.loc[per_species_counts['trait']==1, 'species'].tolist()
        multi_species = per_species_counts.loc[per_species_counts['trait']==0, 'species'].tolist()
        n_multi = len(multi_species)
        if n_multi == 0:
            print("WARNING: no multicellular species found; cannot identify missing OGs.", file=sys.stderr)
        else:
            pres_only = presence.copy()
            pres_only.set_index(ortho_id_col, inplace=True)
            multi_presence_sum = pres_only[multi_species].sum(axis=1)
            uni_presence_sum = pres_only[uni_species].sum(axis=1) if len(uni_species)>0 else pd.Series(0, index=pres_only.index)

            # Condition: uni_presence_sum == 0 AND multi_presence_sum >= threshold_count
            threshold_count = max(1, int(np.ceil(MULTIC_FRAC_THRESHOLD * n_multi)))
            cond = (uni_presence_sum == 0) & (multi_presence_sum >= threshold_count)
            missing_series = cond[cond].index.tolist()

            if len(missing_series) > 0:
                rows = []
                for og in missing_series:
                    multic_count = int(multi_presence_sum.loc[og])
                    multic_frac = float(multic_count) / float(n_multi)
                    present_multic_species = pres_only.loc[og, multi_species][pres_only.loc[og, multi_species] == 1].index.tolist()
                    rows.append({
                        ortho_id_col: og,
                        'multic_present_count': multic_count,
                        'multic_present_fraction': multic_frac,
                        'multic_present_species': ";".join(present_multic_species)
                    })
                missing_ogs_df = pd.DataFrame(rows)
                missing_ogs_df.to_csv(OUT_MISSING, sep="\t", index=False)
            else:
                missing_ogs_df = pd.DataFrame(columns=[ortho_id_col, 'multic_present_count', 'multic_present_fraction', 'multic_present_species'])
                missing_ogs_df.to_csv(OUT_MISSING, sep="\t", index=False)

    per_species_counts.to_csv(OUT_COUNTS, sep="\t", index=False)

    with open(OUT_SUMMARY, 'w') as fh:
        fh.write("Gene family presence comparison summary\n")
        fh.write("======================================\n\n")
        fh.write(f"Input orthogroup file: {ORTHO_FILE}\n")
        fh.write(f"Mapping file: {MAPPING_FILE}\n\n")
        fh.write("Per-group descriptive statistics (present OG counts per species):\n")
        fh.write(f"Unicellular (trait=1): n={uni_desc['n']}, mean={uni_desc['mean']:.2f}, median={uni_desc['median']:.1f}, std={uni_desc['std']:.2f}\n")
        fh.write(f"Multicellular (trait=0): n={multi_desc['n']}, mean={multi_desc['mean']:.2f}, median={multi_desc['median']:.1f}, std={multi_desc['std']:.2f}\n\n")
        fh.write(f"Statistical test: {test_used}\n")
        if test_note:
            fh.write(test_note + "\n")
        fh.write(f"Test statistic: {stat}\n")
        fh.write(f"One-sided p-value (uni < multi): {p_value:.6g}\n")
        fh.write(f"Significance threshold: {P_VALUE_THRESHOLD}\n")
        fh.write(f"Result: {'SIGNIFICANT (unicellular < multicellular)' if significant else 'NOT SIGNIFICANT'}\n\n")

        if significant:
            fh.write(f"Missing OG identification criteria: present in >= {MULTIC_FRAC_THRESHOLD*100:.0f}% of multicellular species (threshold count = {threshold_count}), and absent in all unicellular species.\n")
            if missing_ogs_df is not None and len(missing_ogs_df)>0:
                fh.write(f"Number of OGs meeting criteria: {len(missing_ogs_df)}\n")
                fh.write(f"Saved detailed list to: {OUT_MISSING}\n")
            else:
                fh.write("No OGs met the criteria.\n")
        else:
            fh.write("Since the group comparison was not significant, missing OG identification was not performed.\n")

    print("Analysis completed.")
    print(f"Per-species present-OG counts saved to: {OUT_COUNTS}")
    print(f"Summary saved to: {OUT_SUMMARY}")
    if significant:
        print(f"Statistical test indicates unicellular counts are significantly smaller (p = {p_value:.6g}).")
        print(f"Missing OGs (present in multicellular but absent in unicellular) saved to: {OUT_MISSING}")
    else:
        print(f"No significant difference detected (p = {p_value:.6g}); missing-OG detection skipped.")

if __name__ == "__main__":
    main()
