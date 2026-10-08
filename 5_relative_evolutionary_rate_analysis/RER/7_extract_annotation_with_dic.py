#!/usr/bin/env python3

import os
import glob
import io
import pandas as pd
from Bio import SeqIO

def load_annotation_tables(anno_dir):
    anno_mapping = {}
    pattern = os.path.join(anno_dir, '*.emapper.annotations')
    for filepath in glob.glob(pattern):
        filtered_lines = []
        with open(filepath, 'r', encoding='utf-8') as fh:
            for raw_line in fh:
                if raw_line.startswith('##'):
                    continue
                filtered_lines.append(raw_line)
        if not filtered_lines:
            continue

        content = io.StringIO(''.join(filtered_lines))
        try:
            df = pd.read_csv(content, sep='\t', header=0, dtype=str)
        except Exception:
            continue

        if df.shape[1] < 21:
            continue

        for _, row in df.iterrows():
            key = row.iloc[0]
            if pd.isna(key):
                continue

            row_list = [x if pd.notna(x) else '' for x in row.tolist()]
            anno_mapping[str(key)] = row_list

    return anno_mapping

def annotate_fasta_folder(fasta_dir, anno_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    anno_mapping = load_annotation_tables(anno_dir)

    columns = [
        "#query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
        "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
        "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
        "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
    ]

    for fasta_path in glob.glob(os.path.join(fasta_dir, '*.fa')):
        basename = os.path.basename(fasta_path)
        name, _ = os.path.splitext(basename)
        out_file = os.path.join(output_dir, f"{name}.emapper.annotations")

        rows = []
        for record in SeqIO.parse(fasta_path, 'fasta'):
            header = record.id
            if header in anno_mapping:
                row = anno_mapping[header]
                if len(row) < len(columns):
                    row = row + [''] * (len(columns) - len(row))
                elif len(row) > len(columns):
                    row = row[:len(columns)]

                safe_row = [str(x) if x is not None else '' for x in row]
                rows.append(safe_row)
            else:
                rows.append([header] + [''] * (len(columns) - 1))

        with open(out_file, 'w', encoding='utf-8') as w:
            w.write("\t".join(columns) + "\n")
            for row in rows:
                cleaned = [str(x).replace('\t', ' ').replace('\n', ' ') if x is not None else '' for x in row]
                w.write("\t".join(cleaned) + "\n")

if __name__ == '__main__':
    fasta_directory = './pos2anno'
    # Example: ../../1_data_acquisition_and_functional_annotation/76_anno/
    annotation_dir = '/home/yankh/76_converge/76_anno'
    output_directory = './pos_anno'

    annotate_fasta_folder(fasta_directory, annotation_dir, output_directory)
