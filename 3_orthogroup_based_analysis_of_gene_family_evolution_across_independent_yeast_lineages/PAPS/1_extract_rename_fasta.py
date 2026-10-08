#!/bin/env python3
# -*- coding: UTF-8 -*-

"""Rename protein files and replace numeric species prefixes with PAPS labels."""

import os
from pathlib import Path
import csv
import re
import shutil

# Example: ../../1_data_acquisition_and_functional_annotation/76_faa/
INPUT_DIR = "/home/yankh/76_converge/76_faa"
MAPPING_FILE = (
    Path(__file__).resolve().parents[2] / "1_data_acquisition_and_functional_annotation"
    / "related_files" / "76_fungi_trait.txt"
)
OUTPUT_DIR = str(Path(__file__).resolve().parent / "faa")
LABEL_TXT = str(Path(__file__).resolve().parent / "label.txt")

def safe_filename_from_species(species_name):
    """Build a FASTA filename from the species name."""
    name = species_name.strip()
    name = re.sub(r'\s+', '_', name)
    # Keep these filename stems consistent with order.txt.
    name = re.sub(r'[^A-Za-z0-9_]', '', name)
    return name + ".faa"

def make_label_from_species(species_name):
    """Build a label from the first two genus and first three species letters."""
    s = species_name.strip()
    parts = re.split(r'[\s\-_]+', s)
    if len(parts) >= 2:
        genus = parts[0]
        epithet = parts[1]
    elif len(parts) == 1:
        genus = parts[0]
        epithet = ''
    else:
        genus = ''
        epithet = ''
    genus_part = (genus[:2].capitalize()) if len(genus) >= 2 else genus.capitalize()
    epithet_part = (epithet[:3].capitalize()) if len(epithet) >= 3 else epithet.capitalize()
    label = (genus_part + epithet_part)
    label = re.sub(r'[^A-Za-z0-9]', '', label)
    if label == "":
        label = "Unknown"
    return label

def read_mapping_file(mapping_path):
    """Read filename and numeric-species mappings."""
    mapping_by_filename = {}
    mapping_by_number = {}
    with open(mapping_path, 'r', newline='', encoding='utf-8') as fh:
        sample = fh.read(4096)
        fh.seek(0)
        reader = csv.reader(fh, delimiter='\t')
        rows = list(reader)
        if not rows:
            raise ValueError("Mapping file is empty: {}".format(mapping_path))
        header = rows[0]
        has_header = False
        header_lower = [c.lower() for c in header]
        if 'file_names' in header_lower or 'names' in header_lower and 'numbering' in header_lower:
            has_header = True

        if has_header:
            idx_names = header_lower.index('names') if 'names' in header_lower else None
            idx_numbering = header_lower.index('numbering') if 'numbering' in header_lower else None
            idx_file_names = header_lower.index('file_names') if 'file_names' in header_lower else None
            start_row = 1
        else:
            idx_names, idx_numbering, idx_file_names = 0, 1, 3
            start_row = 0

        for row in rows[start_row:]:
            if len(row) <= max(idx_names or 0, idx_numbering or 0, idx_file_names or 0):
                continue
            try:
                species_name = row[idx_names].strip()
                numbering = str(row[idx_numbering]).strip()
                file_name_field = row[idx_file_names].strip()
            except Exception:
                continue
            mapping_by_filename[file_name_field] = (species_name, numbering)
            mapping_by_number[numbering] = species_name
    return mapping_by_filename, mapping_by_number

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def process_faa_files():
    ensure_dir(OUTPUT_DIR)
    ensure_dir(os.path.dirname(LABEL_TXT))
    print("Reading mapping file: {}".format(MAPPING_FILE))
    mapping_by_filename, mapping_by_number = read_mapping_file(MAPPING_FILE)

    label_by_number = {}
    used_labels = {}
    for number, species in mapping_by_number.items():
        base_label = make_label_from_species(species)
        if base_label in used_labels:
            unique_label = "{}{}".format(base_label, number)
        else:
            unique_label = base_label
        used_labels[unique_label] = True
        label_by_number[str(number)] = unique_label

    label_by_species = {}
    for num, species in mapping_by_number.items():
        label_by_species[species] = label_by_number[str(num)]

    processed_label_records = []

    input_files = sorted(os.listdir(INPUT_DIR))
    for fname in input_files:
        if not fname.endswith('.faa'):
            continue
        in_path = os.path.join(INPUT_DIR, fname)
        base = fname[:-4]

        mapped = mapping_by_filename.get(base)
        if mapped:
            species_from_filename, numbering_for_file = mapped
        else:
            species_from_filename, numbering_for_file = None, None
            print("WARNING: no mapping found for filename base '{}', file skipped: {}".format(base, fname))
            continue

        output_filename_base = safe_filename_from_species(species_from_filename)
        output_path = os.path.join(OUTPUT_DIR, output_filename_base)

        if os.path.exists(output_path):
            output_filename_base_noext = output_filename_base[:-4]
            output_filename_base = "{}_{}.faa".format(output_filename_base_noext, numbering_for_file)
            output_path = os.path.join(OUTPUT_DIR, output_filename_base)

        with open(in_path, 'r', encoding='utf-8') as infh, open(output_path, 'w', encoding='utf-8') as outfh:
            for line in infh:
                if not line:
                    continue
                if line.startswith('>'):
                    m = re.match(r'^>(\d+)\|(.*)', line.rstrip('\n'))
                    if m:
                        number = m.group(1)
                        rest = m.group(2)
                        label = label_by_number.get(str(number))
                        if not label:
                            if numbering_for_file is not None and str(numbering_for_file) in label_by_number:
                                label = label_by_number[str(numbering_for_file)]
                            else:
                                label = make_label_from_species(species_from_filename)
                        new_header = ">{}_{}".format(label, rest)
                        outfh.write(new_header + "\n")
                    else:
                        new_line = line.replace('|', '_', 1) if '|' in line else line
                        outfh.write(new_line)
                else:
                    outfh.write(line)
        species_label = label_by_species.get(species_from_filename, make_label_from_species(species_from_filename))
        processed_label_records.append((output_filename_base, species_label))
        print("Processed: {} -> {}, label={}".format(fname, output_filename_base, species_label))

    with open(LABEL_TXT, 'w', encoding='utf-8') as labfh:
        for outfname, lbl in processed_label_records:
            labfh.write("{}\t{}\n".format(os.path.splitext(outfname)[0], lbl))
    print("Label file written to: {}".format(LABEL_TXT))
    print("Processing complete. Output directory: {}".format(OUTPUT_DIR))

if __name__ == "__main__":
    process_faa_files()
