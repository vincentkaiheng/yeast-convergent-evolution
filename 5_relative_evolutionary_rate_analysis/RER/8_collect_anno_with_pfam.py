#!/usr/bin/env python3

import os
import pandas as pd
from collections import Counter

def read_sequence_list(file_path):
    df = pd.read_csv(file_path, sep='\t', header=0, dtype=str)
    sequences = df['gene'].tolist()
    return sequences

def get_name_description_pfam_from_annotation(sequence, annotation_dir):
    annotation_file = os.path.join(annotation_dir, f"{sequence}.emapper.annotations")

    if not os.path.exists(annotation_file):
        return None, None, None

    # Use the most frequent informative value; ties follow file order.
    names = []
    descriptions = []
    pfams = []

    with open(annotation_file, 'r', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue

            columns = line.rstrip('\n').split('\t')

            if len(columns) > 7 and columns[7].strip().lower() not in {'', '-', 'na', 'nan', 'none'}:
                descriptions.append(columns[7].strip())

            if len(columns) > 8 and columns[8].strip().lower() not in {'', '-', 'na', 'nan', 'none'}:
                names.append(columns[8].strip())

            if len(columns) > 20 and columns[20].strip().lower() not in {'', '-', 'na', 'nan', 'none'}:
                pfams.append(columns[20].strip())

    most_common_name = None
    most_common_description = None
    most_common_pfam = None

    if names:
        most_common_name = Counter(names).most_common(1)[0][0]

    if descriptions:
        most_common_description = Counter(descriptions).most_common(1)[0][0]

    if pfams:
        most_common_pfam = Counter(pfams).most_common(1)[0][0]

    return most_common_name, most_common_description, most_common_pfam

def generate_gene_description_name_pfam_file(sequence_list, annotation_dir, output_file):
    records = []

    for seq in sequence_list:
        name, description, pfam = get_name_description_pfam_from_annotation(
            seq, annotation_dir
        )

        if description is None:
            description = '-'

        if name is None:
            name = '-'

        if pfam is None:
            pfam = '-'

        records.append((seq, description, name, pfam))

    result_df = pd.DataFrame(
        records,
        columns=['Gene', 'Description', 'Name', 'PFAM']
    )

    result_df.to_csv(
        output_file,
        sep='\t',
        index=False
    )

def main():
    sequence_list_file = 'related_files/cor_pos.txt'
    annotation_dir = 'pos_anno'
    output_file = 'related_files/pos_gene.txt'

    sequences = read_sequence_list(sequence_list_file)

    generate_gene_description_name_pfam_file(
        sequences,
        annotation_dir,
        output_file
    )

if __name__ == "__main__":
    main()
