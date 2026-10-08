#!/usr/bin/env python3
"""Collect unique GO terms across each OG and lineage."""

import os
import glob
import sys
from collections import OrderedDict

# Example: ../../3_orthogroup_based_analysis_of_gene_family_evolution_across_independent_yeast_lineages/PAPS/run_PAPS(modified_from_Jordi_Paps)/Output/
INPUT_DIR = "/home/yankh/76_converge/paps/PAPS/Output"
ANNO_DIR = "/home/yankh/76_converge/76_faa/all_OGs_anno"
OUTPUT_DIR = "./GOs"

EXPECTED_PREFIXES = [
    "Pucciniomycotina_yeast-atleast3_non_Pucciniomycotina_yeast-absent",
    "Ustilaginomycotina_yeast-atleast4_non_Ustilaginomycotina_yeast-absent",
    "Agaricomycotina_yeast-atleast3_non_Agaricomycotina_yeast-absent",
    "Taphrinomycotina_yeast-atleast4_non_Taphrinomycotina_yeast-absent",
    "Saccharomycotina_yeast-atleast6_non_Saccharomycotina_yeast-absent",
    "Pezizomycotina_yeast-atleast2_non_Pezizomycotina_yeast-absent"
]

os.makedirs(OUTPUT_DIR, exist_ok=True)

pattern = os.path.join(INPUT_DIR, "*MCL_genes_IDs.out")
all_files = glob.glob(pattern)

files_to_process = []
file_to_taxa = {}

for p in all_files:
    basename = os.path.basename(p)
    for prefix in EXPECTED_PREFIXES:
        if basename.startswith(prefix):
            files_to_process.append(p)
            file_to_taxa[p] = prefix.split("_yeast-atleast", 1)[0]
            break

if not files_to_process:
    print("No matching input files found. Please check INPUT_DIR and file prefixes.", file=sys.stderr)
    sys.exit(1)

print(f"Found {len(files_to_process)} matching input files. Starting processing...")

for filepath in files_to_process:
    basename = os.path.basename(filepath)
    taxonomy = file_to_taxa[filepath]

    print(f"Processing file: {basename} -> Taxonomy: {taxonomy}")

    out_table_path = os.path.join(OUTPUT_DIR, f"{taxonomy}_OG_HG_GOs.tsv")
    out_summary_path = os.path.join(OUTPUT_DIR, f"{taxonomy}_GO_summary.txt")

    all_go_set = set()

    with open(filepath, "r") as fin, open(out_table_path, "w") as fout:
        fout.write("OG\tHG\tGOs\n")
        for lineno, line in enumerate(fin, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) == 0:
                continue
            hg = cols[0].strip()
            if not hg:
                continue

            if not hg.startswith("HG_"):
                print(f"Warning: Non-standard HG identifier found in {basename} at line {lineno}: '{hg}'. Skipping.", file=sys.stderr)
                continue
            try:
                hg_num = int(hg.split("_")[1])
            except Exception as e:
                print(f"Failed to parse HG number ({hg}): {e}", file=sys.stderr)
                continue

            og_num = hg_num - 1
            og_id = f"OG{og_num:07d}"

            anno_path = os.path.join(ANNO_DIR, f"{og_id}.emapper.annotations")
            og_go_set = OrderedDict()

            if os.path.isfile(anno_path):
                try:
                    with open(anno_path, "r") as af:
                        for a_line in af:
                            a_line = a_line.rstrip("\n")
                            if not a_line:
                                continue
                            if a_line.startswith("#"):
                                continue
                            a_cols = a_line.split("\t")
                            if len(a_cols) >= 10:
                                gos_field = a_cols[9].strip()
                                if gos_field and gos_field != "-":
                                    for go in gos_field.split(','):
                                        go = go.strip()
                                        if go:
                                            og_go_set[go] = None
                except Exception as e:
                    print(f"Failed to read annotation file: {anno_path}, Error: {e}", file=sys.stderr)
            else:
                print(f"Warning: Annotation file not found: {anno_path}", file=sys.stderr)

            og_go_list = list(og_go_set.keys())
            go_str = ",".join(og_go_list) if og_go_list else ""

            fout.write(f"{og_id}\t{hg}\t{go_str}\n")

            all_go_set.update(og_go_list)

    with open(out_summary_path, "w") as sf:
        for go in sorted(all_go_set):
            sf.write(go + "\n")

    print(f"Outputs successfully saved to: {out_table_path} and {out_summary_path}")

print("Processing complete.")
