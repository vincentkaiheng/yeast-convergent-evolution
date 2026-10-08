#!/usr/bin/env python3
"""Build an OG presence/absence matrix and exact-k intersections across six yeast lineages."""

import os
import glob
import pandas as pd
from typing import Dict, Set

# Example: ../../3_orthogroup_based_analysis_of_gene_family_evolution_across_independent_yeast_lineages/PAPS/run_PAPS(modified_from_Jordi_Paps)/Output/
input_dir = "/home/yankh/76_converge/paps/PAPS/Output"

expected_prefixes = [
    "Pucciniomycotina_yeast-atleast3_non_Pucciniomycotina_yeast-absent",
    "Ustilaginomycotina_yeast-atleast4_non_Ustilaginomycotina_yeast-absent",
    "Agaricomycotina_yeast-atleast3_non_Agaricomycotina_yeast-absent",
    "Taphrinomycotina_yeast-atleast4_non_Taphrinomycotina_yeast-absent",
    "Saccharomycotina_yeast-atleast6_non_Saccharomycotina_yeast-absent",
    "Pezizomycotina_yeast-atleast2_non_Pezizomycotina_yeast-absent"
]
expected_classes = [prefix.split("_yeast-atleast", 1)[0] for prefix in expected_prefixes]

file_glob_pattern = os.path.join(input_dir, "*MCL_genes_IDs.out")

output_dir = "./OGs_upset_outputs"
os.makedirs(output_dir, exist_ok=True)

og_offset = 1

upset_table_path = os.path.join(output_dir, "upset_input_OG.tsv")
intersection_prefix = os.path.join(output_dir, "intersect_exact_k_")

def extract_class_from_basename(basename: str) -> str:
    """Return the lineage matching a PAPS query filename."""

    for prefix in expected_prefixes:
        if basename.startswith(prefix):
            return prefix.split("_yeast-atleast", 1)[0]
    return None

def parse_hg_from_line(line: str) -> str:
    """Read the HG identifier from the first tab-separated field."""

    if not line:
        return None
    parts = line.strip().split("\t")
    if not parts:
        return None
    hg = parts[0].strip()
    if not hg:
        return None
    return hg

def hg_to_og(hg: str, offset: int = 1) -> str:
    """Map HG_n to OG(n - 1) for the preserved OrthoFinder row order."""

    if not hg.startswith("HG_"):
        raise ValueError(f"Unexpected HG format: '{hg}'")
    try:
        hg_num = int(hg.split("_", 1)[1])
    except Exception as e:
        raise ValueError(f"Cannot parse integer from HG '{hg}': {e}")
    og_num = hg_num - offset
    if og_num < 0:
        raise ValueError(f"Computed OG number negative for HG '{hg}' with offset {offset}")
    return f"OG{og_num:07d}"

def main():
    files = glob.glob(file_glob_pattern)
    if not files:
        raise FileNotFoundError(f"No files found with pattern: {file_glob_pattern}")

    class_hg_sets: Dict[str, Set[str]] = dict()

    for f in files:
        basename = os.path.basename(f)
        cls = extract_class_from_basename(basename)
        if not cls:
            continue

        hg_set = set()
        with open(f, "r") as fh:
            for line in fh:
                if not line.strip():
                    continue
                if line.startswith("#"):
                    continue
                hg = parse_hg_from_line(line)
                if hg:
                    hg_set.add(hg)

        class_hg_sets[cls] = hg_set

    missing = [c for c in expected_classes if c not in class_hg_sets]
    if missing:
        print(f"Warning: the following expected classes were not found and will be included as empty: {missing}")
        for c in missing:
            class_hg_sets[c] = set()

    def hg_key(hg_str: str) -> int:
        try:
            return int(hg_str.split("_", 1)[1])
        except:
            return float("inf")

    all_hgs = sorted({hg for s in class_hg_sets.values() for hg in s}, key=hg_key)

    df = pd.DataFrame(index=all_hgs)
    for cls in expected_classes:
        presence = [1 if hg in class_hg_sets.get(cls, set()) else 0 for hg in all_hgs]
        df[cls] = presence

    df.insert(0, "HG", df.index)

    og_list = []
    for hg in df["HG"].tolist():
        try:
            og = hg_to_og(hg, offset=og_offset)
        except Exception as e:
            raise RuntimeError(f"Error converting HG '{hg}' to OG: {e}")
        og_list.append(og)
    df.insert(0, "OG", og_list)

    df.to_csv(upset_table_path, sep="\t", index=False)
    print(f"Wrote UpSet input table (OG + presence/absence) to: {upset_table_path}")

    presence_cols = expected_classes

    df["n_sets"] = df[presence_cols].sum(axis=1).astype(int)

    for k in range(2, 7):
        df_k = df[df["n_sets"] == k].copy()
        out_k = f"{intersection_prefix}{k}.tsv"

        cols_to_write = ["OG", "HG"] + presence_cols
        df_k.to_csv(out_k, sep="\t", index=False, columns=cols_to_write)
        print(f"Wrote {len(df_k)} entries in exactly {k} sets to: {out_k}")

    print("Processing complete.")

if __name__ == "__main__":
    main()
