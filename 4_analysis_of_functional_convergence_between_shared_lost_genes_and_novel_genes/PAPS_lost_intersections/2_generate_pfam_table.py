#!/usr/bin/env python3
"""Collect unique PFAM annotations across each OG and lineage."""

import os
import glob
import sys
from collections import OrderedDict

# Example: ../../3_orthogroup_based_analysis_of_gene_family_evolution_across_independent_yeast_lineages/PAPS/run_PAPS(modified_from_Jordi_Paps)/Output/
INPUT_DIR  = "/home/yankh/76_converge/paps/PAPS/Output"
ANNO_DIR   = "/home/yankh/76_converge/76_faa/all_OGs_anno"
OUTPUT_DIR = "./PFAMs"

os.makedirs(OUTPUT_DIR, exist_ok=True)

pattern = os.path.join(INPUT_DIR, "*MCL_genes_IDs.out")
files = [p for p in glob.glob(pattern) if "yeast-absent_si" in os.path.basename(p)]

if not files:
    print("No matching input files found. Check INPUT_DIR and filename pattern.", file=sys.stderr)
    sys.exit(1)

print(f"Found {len(files)} matching input files. Starting processing...")

for filepath in files:
    basename = os.path.basename(filepath)
    parts = basename.split("_")

    taxon = None
    if "yeast-absent_si" in parts:
        idx = parts.index("yeast-absent_si")
        if idx >= 1:
            taxon = parts[idx - 1]

    if taxon is None or not taxon:
        taxon = parts[0]

    print(f"Processing: {basename} -> taxon: {taxon}")

    table_path   = os.path.join(OUTPUT_DIR, f"{taxon}_OG_HG_PFAMs.tsv")
    summary_path = os.path.join(OUTPUT_DIR, f"{taxon}_PFAM_summary.txt")

    all_pfams = set()

    with open(filepath, "r") as fin, open(table_path, "w") as fout:
        fout.write("OG\tHG\tPFAMs\n")

        for lineno, line in enumerate(fin, start=1):
            line = line.rstrip("\n")
            if not line.strip():
                continue

            cols = line.split("\t")
            if not cols:
                continue

            hg = cols[0].strip()
            if not hg.startswith("HG_"):
                print(f"Warning: line {lineno} in {basename} has invalid HG format: '{hg}' -> skipped", file=sys.stderr)
                continue

            try:
                hg_num = int(hg.split("_")[1])
            except ValueError:
                print(f"Error: cannot parse HG number from '{hg}' (line {lineno}, {basename})", file=sys.stderr)
                continue

            og_num = hg_num - 1
            og_id  = f"OG{og_num:07d}"
            anno_file = os.path.join(ANNO_DIR, f"{og_id}.emapper.annotations")

            pfam_dict = OrderedDict()

            if os.path.isfile(anno_file):
                try:
                    with open(anno_file, "r") as af:
                        for a_line in af:
                            a_line = a_line.rstrip("\n")
                            if not a_line or a_line.startswith("#"):
                                continue
                            a_cols = a_line.split("\t")
                            if len(a_cols) < 21:
                                continue
                            pfam_field = a_cols[20].strip()
                            if pfam_field and pfam_field != "-" and pfam_field.lower() != "n/a":
                                for pfam in pfam_field.split(","):
                                    pfam = pfam.strip()
                                    if pfam:
                                        pfam_dict[pfam] = None
                except Exception as e:
                    print(f"Error reading annotation file {anno_file}: {e}", file=sys.stderr)
            else:
                print(f"Warning: annotation file not found: {anno_file}", file=sys.stderr)

            pfam_list = list(pfam_dict.keys())
            pfam_str  = ",".join(pfam_list) if pfam_list else ""
            fout.write(f"{og_id}\t{hg}\t{pfam_str}\n")

            all_pfams.update(pfam_list)

    with open(summary_path, "w") as sf:
        for pfam in sorted(all_pfams):
            sf.write(pfam + "\n")

    print(f"Finished: {table_path}\n          {summary_path}")

print("All files processed successfully!")
