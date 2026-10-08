#!/bin/env python3
# -*- coding: utf-8 -*-
"""Split the gene-family table using the species in each subtree."""

import re
from pathlib import Path
import pandas as pd
import sys

INPUT_TABLE = Path('./cafe_related_files/gene_families_filtered.txt')
TREE_ASCO = Path('./asco/tree_asco.txt')
TREE_BASI = Path('./basi/tree_basi.txt')
OUT_ASCO = Path('./asco/gene_family_asco.txt')
OUT_BASI = Path('./basi/gene_family_basi.txt')

def parse_newick_species(newick_text):
    """Extract labels followed by branch lengths from Newick text."""
    labels = re.findall(r'([A-Za-z0-9_.\-]+):', newick_text)
    return set(labels)

def read_tree_file(path):
    """Read species labels from a tree file."""
    if not path.exists():
        raise FileNotFoundError(f"Tree file not found: {path}")
    text = path.read_text(encoding='utf-8').strip()
    if not text:
        return set()
    return parse_newick_species(text)

def select_and_write(df, species_set, out_path):
    """Write species columns for a subtree and report unmatched names."""
    df = df.copy()
    df.columns = df.columns.str.strip()
    required_cols = ['desc', 'Orthogroup']
    for rc in required_cols:
        if rc not in df.columns:
            raise KeyError(f"Required column '{rc}' not found in input table.")

    present_species = [c for c in df.columns if c in species_set]
    missing_species = sorted(species_set - set(df.columns))

    out_cols = required_cols + present_species
    df_out = df[out_cols]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_path, sep='\t', index=False, encoding='utf-8')
    return present_species, missing_species, out_path

def main():
    try:
        if not INPUT_TABLE.exists():
            print(f"ERROR: Input table not found: {INPUT_TABLE}", file=sys.stderr)
            sys.exit(1)
        df = pd.read_csv(INPUT_TABLE, sep='\t', header=0, dtype=str, keep_default_na=False)
        asco_species = read_tree_file(TREE_ASCO)
        basi_species = read_tree_file(TREE_BASI)

        asco_present, asco_missing, asco_out = select_and_write(df, asco_species, OUT_ASCO)
        print(f"Wrote ASCO output to: {asco_out}")
        print(f"ASCO: {len(asco_present)} species columns found and written.")
        if asco_missing:
            print(f"ASCO: {len(asco_missing)} species from tree not found in table. Example missing names: {asco_missing[:10]}")
        else:
            print("ASCO: all species from tree were found in the table.")

        basi_present, basi_missing, basi_out = select_and_write(df, basi_species, OUT_BASI)
        print(f"Wrote BASI output to: {basi_out}")
        print(f"BASI: {len(basi_present)} species columns found and written.")
        if basi_missing:
            print(f"BASI: {len(basi_missing)} species from tree not found in table. Example missing names: {basi_missing[:10]}")
        else:
            print("BASI: all species from tree were found in the table.")

    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        sys.exit(2)

if __name__ == '__main__':
    main()
