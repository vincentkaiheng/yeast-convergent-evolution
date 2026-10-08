#!/bin/env python3
# -*- coding: utf-8 -*-
"""Convert an OrthoFinder table to one line of PAPS sequence IDs per orthogroup."""

import os
from pathlib import Path
import csv
import re
import sys

INPUT_OG_FILE = str(Path(__file__).resolve().parents[1] / "orthofinder_result" / "Orthogroups.tsv")
MAPPING_FILE = (
    Path(__file__).resolve().parents[2] / "1_data_acquisition_and_functional_annotation"
    / "related_files" / "76_fungi_trait.txt"
)
OUTPUT_DIR = str(Path(__file__).resolve().parent)
OUTPUT_OG_FILE = os.path.join(OUTPUT_DIR, os.path.basename(INPUT_OG_FILE))


def make_label_from_species(species_name):
    """Build a label from the first two genus and first three species letters."""
    s = species_name.strip()
    if not s:
        return "Unknown"
    parts = re.split(r'[\s\-_]+', s)
    genus = parts[0] if len(parts) >= 1 else ''
    epithet = parts[1] if len(parts) >= 2 else ''
    genus_part = (genus[:2].capitalize()) if len(genus) >= 2 else genus.capitalize()
    epithet_part = (epithet[:3].capitalize()) if len(epithet) >= 3 else epithet.capitalize()
    label = (genus_part + epithet_part)
    label = re.sub(r'[^A-Za-z0-9]', '', label)
    if not label:
        return "Unknown"
    return label


def read_mapping_file(mapping_path):
    """Read filename and numeric-species mappings."""
    mapping_by_filename = {}
    mapping_by_number = {}
    if not os.path.exists(mapping_path):
        raise FileNotFoundError(f"Mapping file not found: {mapping_path}")

    with open(mapping_path, 'r', encoding='utf-8') as fh:
        reader = csv.reader(fh, delimiter='\t')
        rows = [r for r in reader if any(cell.strip() for cell in r)]

    if not rows:
        raise ValueError("Mapping file is empty or contains only blank lines: {}".format(mapping_path))

    header = rows[0]
    header_lower = [c.strip().lower() for c in header]
    has_header = ('names' in header_lower and 'numbering' in header_lower and 'file_names' in header_lower)

    if has_header:
        idx_names = header_lower.index('names')
        idx_numbering = header_lower.index('numbering')
        idx_file_names = header_lower.index('file_names')
        start = 1
    else:
        idx_names, idx_numbering, idx_file_names = 0, 1, 3
        start = 0

    for row in rows[start:]:
        if len(row) <= max(idx_names, idx_numbering, idx_file_names):
            continue
        species_name = row[idx_names].strip()
        numbering = str(row[idx_numbering]).strip()
        file_names_field = row[idx_file_names].strip() if len(row) > idx_file_names else ''
        if species_name and numbering:
            mapping_by_number[numbering] = species_name
        if file_names_field:
            mapping_by_filename[file_names_field] = (species_name, numbering)

    return mapping_by_filename, mapping_by_number


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def build_label_map(mapping_by_number):
    """Assign unique species labels, appending numeric IDs for duplicates."""
    label_by_number = {}
    used_labels = {}
    for numbering, species in mapping_by_number.items():
        base_label = make_label_from_species(species)
        if base_label in used_labels:
            unique_label = f"{base_label}{numbering}"
        else:
            unique_label = base_label
        used_labels[unique_label] = True
        label_by_number[str(numbering)] = unique_label
    return label_by_number


def process_token(token, label_by_number):
    """Replace a numeric species prefix with its PAPS label."""
    token = token.strip()
    if not token:
        return ''
    m = re.match(r'^(\d+)\|(.*)$', token)
    if m:
        num = m.group(1)
        rest = m.group(2).strip()
        label = label_by_number.get(str(num))
        if not label:
            label = f"Sp{num}"
        return f"{label}_{rest}"
    else:
        if '|' in token:
            return token.replace('|', '_', 1)
        return token


def process_cell(cell, label_by_number):
    """Convert comma- or semicolon-separated sequence identifiers."""
    if cell is None:
        return []
    cell_str = str(cell).strip()
    if cell_str == '' or cell_str == '-' or cell_str.lower() == 'na':
        return []
    tokens = re.split(r'\s*,\s*|\s*;\s*', cell_str)
    out_tokens = []
    for tok in tokens:
        if not tok:
            continue
        processed = process_token(tok, label_by_number)
        if processed:
            out_tokens.append(processed)
    return out_tokens


def process_orthogroups(input_path, output_path, mapping_path):
    ensure_dir(os.path.dirname(output_path))

    mapping_by_filename, mapping_by_number = read_mapping_file(mapping_path)
    label_by_number = build_label_map(mapping_by_number)

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input orthogroups file not found: {input_path}")

    with open(input_path, 'r', encoding='utf-8') as infh:
        reader = csv.reader(infh, delimiter='\t')
        rows = list(reader)

    if not rows or len(rows) <= 1:
        raise ValueError("Input Orthogroups.tsv appears empty or contains only a header: {}".format(input_path))

    data_rows = rows[1:]

    processed_lines = []
    for row in data_rows:
        if not row:
            processed_lines.append('')
            continue
        remaining = row[1:]
        all_tokens = []
        for cell in remaining:
            tokens = process_cell(cell, label_by_number)
            if tokens:
                all_tokens.extend(tokens)
        merged = ' '.join(all_tokens)
        processed_lines.append(merged)

    with open(output_path, 'w', encoding='utf-8') as outfh:
        for line in processed_lines:
            outfh.write(line + "\n")

    print(f"Processed orthogroups written to: {output_path}")


if __name__ == '__main__':
    try:
        process_orthogroups(INPUT_OG_FILE, OUTPUT_OG_FILE, MAPPING_FILE)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    sys.exit(0)
