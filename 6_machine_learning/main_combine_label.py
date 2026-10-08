#!/usr/bin/env python3

import os
import sys
import pandas as pd

orthogroups_path = '../3_orthogroup_based_analysis_of_gene_family_evolution_across_independent_yeast_lineages/orthofinder_result/Orthogroups.GeneCount.tsv'
ml_label_path   = './related_files/ml_55_matrix_main_label.txt'
output_path    = './related_files/ml_55_matrix_main_labeled.tsv'

EXPECTED_SPECIES_COUNT = 55

def main():
    ml = pd.read_csv(ml_label_path, sep="\t", header=0)

    required_columns = {"species_name", "main_label"}
    missing_columns = required_columns - set(ml.columns)
    if missing_columns:
        raise ValueError(
            f"Label file is missing columns: {', '.join(sorted(missing_columns))}"
        )

    ml = ml[["species_name", "main_label"]].copy()
    ml["species_name"] = ml["species_name"].astype(str).str.strip()
    ml["main_label"] = pd.to_numeric(ml["main_label"], errors="raise").astype(int)

    if len(ml) != EXPECTED_SPECIES_COUNT:
        raise ValueError(
            f"Main labels must contain {EXPECTED_SPECIES_COUNT} species; "
            f"found {len(ml)}."
        )

    if ml["species_name"].duplicated().any():
        duplicates = ml.loc[
            ml["species_name"].duplicated(), "species_name"
        ].tolist()
        raise ValueError(f"Duplicate species in label file: {duplicates}")

    og = pd.read_csv(orthogroups_path, sep="\t", header=0)

    if "Orthogroup" not in og.columns:
        raise ValueError("Orthogroup column is missing from the gene-count matrix.")

    og_t = og.set_index("Orthogroup").T.reset_index()
    og_t = og_t.rename(columns={"index": "species_full"})
    og_t["species_full"] = og_t["species_full"].astype(str).str.strip()

    og_t = og_t[
        og_t["species_full"].str.lower() != "total"
    ].reset_index(drop=True)

    feature_columns = [
        column for column in og_t.columns if column != "species_full"
    ]

    matched_rows = []
    unmatched = []
    ambiguous = []

    species_full_lower = og_t["species_full"].str.lower()

    for _, row in ml.iterrows():
        species_name = row["species_name"]
        species_name_lower = species_name.lower()

        mask = species_full_lower.str.contains(
            species_name_lower,
            case=False,
            regex=False,
            na=False,
        )

        match_count = mask.sum()

        if match_count == 1:
            counts = og_t.loc[mask, feature_columns].iloc[0].to_dict()

            matched_row = {
                "species_name": species_name,
                "main_label": row["main_label"],
            }
            matched_row.update(counts)
            matched_rows.append(matched_row)

        elif match_count == 0:
            unmatched.append(species_name)

        else:
            candidates = og_t.loc[mask, "species_full"].tolist()
            ambiguous.append((species_name, candidates))

    if unmatched or ambiguous:
        if unmatched:
            print("Unmatched species:")
            for species_name in unmatched:
                print(f"  - {species_name}")

        if ambiguous:
            print("\nSpecies matching multiple matrix columns:")
            for species_name, candidates in ambiguous:
                print(f"  - {species_name}: {candidates}")

        sys.exit("Species matching failed; no output was written.")

    if len(matched_rows) != EXPECTED_SPECIES_COUNT:
        sys.exit(
            f"Matched {len(matched_rows)} species; expected "
            f"{EXPECTED_SPECIES_COUNT}. No output was written."
        )

    result = pd.DataFrame(matched_rows)
    output_columns = ["species_name", "main_label"] + feature_columns
    result = result[output_columns]

    result.to_csv(output_path, sep="\t", index=False)

    print(f"Selected {len(result)} main-analysis species.")
    print(f"Labeled matrix saved to: {os.path.abspath(output_path)}")

if __name__ == "__main__":
    main()
