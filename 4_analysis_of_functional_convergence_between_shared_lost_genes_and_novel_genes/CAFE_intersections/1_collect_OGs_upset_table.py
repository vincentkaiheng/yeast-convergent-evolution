#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Build six-lineage CAFE intersections.")
parser.add_argument("--change", choices=["expanded", "contracted"], default="expanded")
change = parser.parse_args().change

node_mapping = {
    "node37": "Pucciniomycotina_yeast",
    "node45_asco": "Taphrinomycotina_yeast",
    "node45_basi": "Ustilaginomycotina_yeast",
    "node52": "Agaricomycotina_yeast",
    "node53": "Saccharomycotina_yeast",
    "node82": "Pezizomycotina_yeast"
}

output_dir = f"./OGs_{change}_upset_outputs"
os.makedirs(output_dir, exist_ok=True)

upset_table_path = os.path.join(output_dir, "og_upset_input.tsv")
intersection_prefix = os.path.join(output_dir, "og_intersect_exact_k_")

def main():
    class_og_sets = {}
    all_og_pool = set()

    print("Reading input files...")
    for node, class_name in node_mapping.items():
        filename = f"{node}.{change}"

        if not os.path.exists(filename):
            raise FileNotFoundError(filename)

        with open(filename, 'r') as f:
            ogs = {line.strip() for line in f if line.strip() and line.strip() != "FamilyID"}
            invalid = [og for og in ogs if not re.fullmatch(r"OG\d{7}", og)]
            if invalid:
                raise ValueError(f"Invalid OG identifiers in {filename}: {invalid[:5]}")
            class_og_sets[class_name] = ogs
            all_og_pool.update(ogs)
            print(f"  - {class_name}: Found {len(ogs)} OGs")

    sorted_ogs = sorted(list(all_og_pool))

    df = pd.DataFrame(index=sorted_ogs)

    ordered_classes = list(node_mapping.values())
    for cls in ordered_classes:
        df[cls] = [1 if og in class_og_sets[cls] else 0 for og in sorted_ogs]

    df.index.name = "OG"
    df.to_csv(upset_table_path, sep="\t")
    print(f"\nSaved UpSet input table to: {upset_table_path}")

    df["n_sets"] = df.sum(axis=1)

    for k in range(2, 7):
        df_k = df[df["n_sets"] == k].copy()
        if not df_k.empty:
            out_k = f"{intersection_prefix}{k}.tsv"

            df_k.drop(columns=["n_sets"]).to_csv(out_k, sep="\t")
            print(f"Wrote {len(df_k)} OGs found in exactly {k} sets to: {out_k}")
        else:
            print(f"No OGs found in exactly {k} sets.")

    print("\nProcessing complete.")

if __name__ == "__main__":
    main()
