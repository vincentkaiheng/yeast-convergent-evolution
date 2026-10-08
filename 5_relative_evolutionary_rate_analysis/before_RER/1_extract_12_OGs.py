#!/usr/bin/env python3

from __future__ import annotations
import os
import sys
import csv
import glob
import shutil

# Example: ../../3_orthogroup_based_analysis_of_gene_family_evolution_across_independent_yeast_lineages/orthofinder_result/Orthogroups.GeneCount.tsv
MATRIX_PATH = "/home/yankh/76_converge/76_faa/OrthoFinder/Results_Jun15_1/Orthogroups/Orthogroups.GeneCount.tsv"
SEQDIR = "/home/yankh/76_converge/76_faa/OrthoFinder/Results_Jun15_1/Orthogroup_Sequences"
OUTDIR = "./OGs_13filtered"
THRESHOLD = 12  # Retain OGs present in at least 13 species.

def count_species_presence(counts: list[str]) -> int:
    present = 0
    for c in counts:
        if c is None:
            continue
        s = str(c).strip()
        if s == "" or s == "0":
            continue

        try:
            val = int(float(s))
        except Exception:
            val = 1
        if val > 0:
            present += 1
    return present

def find_fasta_for_og(seqdir: str, og_id: str) -> str | None:
    candidates = [
        os.path.join(seqdir, og_id + ".fa"),
        os.path.join(seqdir, og_id + ".fasta"),
        os.path.join(seqdir, og_id + ".fa.gz"),
        os.path.join(seqdir, og_id + ".fasta.gz"),
        os.path.join(seqdir, og_id + "*.fa*"),
    ]
    for pat in candidates:
        matches = glob.glob(pat)
        if matches:
            return os.path.abspath(matches[0])
    return None

def validate_paths():
    if not os.path.isfile(MATRIX_PATH):
        print(f"[ERROR] Gene-count matrix not found: {MATRIX_PATH}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(SEQDIR):
        print(f"[ERROR] Sequence directory not found: {SEQDIR}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR, exist_ok=True)

def main():
    validate_paths()
    selected = []
    missing_fasta = []
    copied = []

    with open(MATRIX_PATH, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        try:
            header = next(reader)
        except StopIteration:
            print("[ERROR] Gene-count matrix is empty.", file=sys.stderr)
            sys.exit(1)

        species_columns = [i for i, name in enumerate(header) if i > 0 and name.strip().lower() != "total"]

        for row in reader:
            if not row:
                continue
            og = row[0].strip()
            if og == "":
                continue
            counts = [row[i] for i in species_columns]
            present_count = count_species_presence(counts)
            if present_count > THRESHOLD:
                selected.append((og, present_count))

    selected.sort(key=lambda x: x[1], reverse=True)
    print(f"Selected {len(selected)} OGs present in more than {THRESHOLD} species.")

    for og, cnt in selected:
        fasta = find_fasta_for_og(SEQDIR, og)
        if fasta is None:
            missing_fasta.append((og, cnt))
            print(f"[WARN] FASTA not found: {og} ({cnt} species)")
            continue
        dest = os.path.join(OUTDIR, os.path.basename(fasta))
        try:
            shutil.copy2(fasta, dest)
            copied.append((og, fasta, dest))

        except Exception as e:
            print(f"[ERROR] Copy failed: {fasta} -> {dest} : {e}", file=sys.stderr)

    print("\n=== Summary ===")
    print(f"Selected OGs: {len(selected)}")
    print(f"Files copied: {len(copied)}")
    if missing_fasta:
        print(f"Missing OG FASTA files: {len(missing_fasta)}. Up to 20 examples:")
        for og, cnt in missing_fasta[:20]:
            print(f"  {og}\t{cnt}")
        print("Check that FASTA filenames match the OG IDs.")
    else:
        print("FASTA files were found for all selected OGs.")

if __name__ == "__main__":
    main()
