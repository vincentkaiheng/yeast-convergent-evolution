#!/usr/bin/env python3

from typing import List, Tuple, Set
import os
import sys
import csv
import re

TABLE_PATH = "./related_files/cor_neg.txt"
ORTHO_DIR = "../before_RER/shared_orthologs"
OG_DIR = "../before_RER/OGs_13filtered"
OUT_DIR = "./neg2anno"

def read_genes_from_table(path: str) -> List[str]:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input table not found: {path}")
    genes = []
    with open(path, "r", encoding="utf-8") as fh:
        fh.seek(0)
        try:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames and any(fn.strip() == "gene" for fn in reader.fieldnames):
                for row in reader:
                    g = row.get("gene")
                    if g and g.strip():
                        genes.append(g.strip())
                return genes
        except Exception:
            pass

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0].strip() not in {"gene", "Rho", ""}:
                genes.append(parts[0].strip())
    return genes

def fasta_iter(handle):
    header = None
    seq_lines = []
    for raw in handle:
        line = raw.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_lines)
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        yield header, "".join(seq_lines)

def extract_ids_from_ortho(ortho_path: str) -> Set[str]:
    ids = set()
    with open(ortho_path, "r", encoding="utf-8") as fh:
        for header, _seq in fasta_iter(fh):
            if not header:
                continue
            seq_id = header.split()[0]
            ids.add(seq_id)
    return ids

def find_og_path_from_gene(gene: str, og_dir: str) -> str:
    core = re.sub(r"^aligned_", "", gene)
    core = re.sub(r"_pruned$", "", core)
    candidates = [os.path.join(og_dir, core + ext) for ext in (".fa", ".fasta")]
    for p in candidates:
        if os.path.isfile(p):
            return p

    m = re.search(r"(OG\d+)", core)
    if m:
        ogcore = m.group(1)
        candidates = [os.path.join(og_dir, ogcore + ext) for ext in (".fa", ".fasta")]
        for p in candidates:
            if os.path.isfile(p):
                return p
    return ""

def write_fasta(out_path: str, records: List[Tuple[str,str]]):
    with open(out_path, "w", encoding="utf-8") as outfh:
        for header, seq in records:
            outfh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                outfh.write(seq[i:i+80] + "\n")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    try:
        genes = read_genes_from_table(TABLE_PATH)
    except Exception as e:
        print(f"[ERROR] Cannot read table {TABLE_PATH}: {e}", file=sys.stderr)
        sys.exit(1)
    if not genes:
        print(f"[ERROR] No gene IDs found in {TABLE_PATH}.", file=sys.stderr)
        sys.exit(1)

    total = len(genes)
    count_written = 0
    missing_ortho = []
    missing_og = []
    no_matches = []

    for gene in genes:
        print(f"[INFO] Processing: {gene}")

        ortho_path = ""
        for ext in (".fa", ".fasta"):
            cand = os.path.join(ORTHO_DIR, gene + ext)
            if os.path.isfile(cand):
                ortho_path = cand
                break
        if not ortho_path:
            missing_ortho.append(gene)
            print(f"  [WARN] Ortholog alignment not found: {os.path.join(ORTHO_DIR, gene + '.fa')} (checked .fa and .fasta)")
            continue

        try:
            ids = extract_ids_from_ortho(ortho_path)
        except Exception as e:
            print(f"  [ERROR] Cannot parse ortholog alignment ({ortho_path}): {e}")
            continue
        if not ids:
            print(f"  [WARN] No sequence IDs found in: {ortho_path}")

        og_path = find_og_path_from_gene(gene, OG_DIR)
        if not og_path:
            missing_og.append(gene)
            print(f"  [WARN] Original OG file not found in {OG_DIR}")
            continue

        matched = []
        try:
            with open(og_path, "r", encoding="utf-8") as ogfh:
                for header, seq in fasta_iter(ogfh):
                    if not header:
                        continue
                    seq_id = header.split()[0]
                    if seq_id in ids:
                        matched.append((header, seq))
        except Exception as e:
            print(f"  [ERROR] Cannot parse OG file ({og_path}): {e}")
            continue

        if matched:
            out_path = os.path.join(OUT_DIR, gene + ".fa")
            write_fasta(out_path, matched)
            count_written += 1
            print(f"  [OK] Wrote {len(matched)} sequences -> {out_path}")
        else:
            no_matches.append(gene)
            print(f"  [WARN] No matching sequences in OG file: {og_path} (requested IDs: {len(ids)})")

    print("\n=== Summary ===")
    print(f"Ortholog sets: {total}")
    print(f"Files written: {count_written} (output directory: {os.path.abspath(OUT_DIR)})")
    if missing_ortho:
        print(f"Missing ortholog alignments: {len(missing_ortho)} (examples: {missing_ortho[:10]})")
    if missing_og:
        print(f"Missing original OG files: {len(missing_og)} (examples: {missing_og[:10]})")
    if no_matches:
        print(f"OG files without matching sequences: {len(no_matches)} (examples: {no_matches[:10]})")
    print("Done.")

if __name__ == "__main__":
    main()
