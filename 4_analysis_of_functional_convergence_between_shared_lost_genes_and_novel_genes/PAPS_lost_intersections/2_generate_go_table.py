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

ALLOWED_TAXA = {
    "Pucciniomycotina",
    "Ustilaginomycotina",
    "Agaricomycotina",
    "Taphrinomycotina",
    "Saccharomycotina",
    "Pezizomycotina"
}

os.makedirs(OUTPUT_DIR, exist_ok=True)

pattern = os.path.join(INPUT_DIR, "*MCL_genes_IDs.out")
files = [p for p in glob.glob(pattern) if "yeast-absent_sister" in os.path.basename(p)]

if not files:
    print("No matching files found. Check INPUT_DIR and filename patterns.", file=sys.stderr)
    sys.exit(1)

print(f"Found {len(files)} matching input files.")

for filepath in files:
    basename = os.path.basename(filepath)
    parts = basename.split("_")
    taxonomy = None

    if "yeast-absent_si" in parts:
        idx = parts.index("yeast-absent_si")
        if idx >= 1:
            taxonomy = parts[idx - 1]

    if taxonomy is None:
        taxonomy = parts[0]

    print(f"Processing file: {basename} -> Lineage: {taxonomy}")

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
                print(f"Warning: in {basename} at line {lineno} invalid HG identifier: '{hg}'; skipping.", file=sys.stderr)
                continue
            try:
                hg_num = int(hg.split("_")[1])
            except Exception as e:
                print(f"Cannot parse HG identifier ({hg}): {e}", file=sys.stderr)
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
                    print(f"Cannot read annotation file: {anno_path}, error: {e}", file=sys.stderr)
            else:
                print(f"Warning: annotation file not found: {anno_path}", file=sys.stderr)

            og_go_list = list(og_go_set.keys())
            if og_go_list:
                go_str = ",".join(og_go_list)
            else:
                go_str = ""

            fout.write(f"{og_id}\t{hg}\t{go_str}\n")

            all_go_set.update(og_go_list)

    with open(out_summary_path, "w") as sf:
        for go in sorted(all_go_set):
            sf.write(go + "\n")

    print(f"Wrote: {out_table_path} and {out_summary_path}")

print("Processing complete.")
