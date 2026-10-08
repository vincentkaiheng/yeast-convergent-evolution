#!/usr/bin/env python3

import os
import pandas as pd
from collections import Counter

INPUT_FILE = '../OGs_upset_outputs/upset_input_OG.tsv'
ANNOTATION_DIR = '/home/yankh/76_converge/76_faa/all_OGs_anno'
OUTPUT_FILE = './OGs_detail.tsv'

def get_annotation_info(og_id, annotation_dir):
    """Return modal nonempty names, descriptions and PFAM strings."""

    annotation_file = os.path.join(
        annotation_dir,
        f"{og_id}.emapper.annotations"
    )

    if not os.path.exists(annotation_file):
        return 'No_annotation', 'No_description', 'No_PFAM'

    names = []
    descriptions = []
    pfams = []

    try:
        with open(annotation_file, 'r', encoding='utf-8') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue

                columns = line.rstrip('\n').split('\t')

                if len(columns) > 7:
                    value = columns[7].strip()
                    if value and value != '-':
                        descriptions.append(value)

                if len(columns) > 8:
                    value = columns[8].strip()
                    if value and value != '-':
                        names.append(value)

                if len(columns) > 20:
                    value = columns[20].strip()
                    if value and value != '-':
                        pfams.append(value)

    except Exception as e:
        print(
            f"Warning: Could not read annotation file "
            f"{annotation_file}. Error: {e}"
        )
        return 'Read_Error', 'Read_Error', 'Read_Error'

    most_common_name = (
        Counter(names).most_common(1)[0][0]
        if names else 'No_annotation'
    )

    most_common_description = (
        Counter(descriptions).most_common(1)[0][0]
        if descriptions else 'No_description'
    )

    most_common_pfam = (
        Counter(pfams).most_common(1)[0][0]
        if pfams else 'No_PFAM'
    )

    return most_common_name, most_common_description, most_common_pfam

def process_upset_input(input_file, annotation_dir, output_file):
    """Annotate OGs shared by at least two lineages."""

    print(f"Reading input file: {input_file}")

    try:
        df = pd.read_csv(
            input_file,
            sep='\t',
            dtype=str
        )

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file}")
        return

    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    series_cols = df.columns[2:].tolist()

    result_data = []

    for _, row in df.iterrows():
        og_id = row['OG']

        present_series = []

        for col in series_cols:
            if str(row[col]) == '1':
                present_series.append(col)

        count = len(present_series)
        combo = '&'.join(present_series)

        if count <= 1:
            continue

        name, description, pfam = get_annotation_info(
            og_id,
            annotation_dir
        )

        result_data.append({
            'OG': og_id,
            'count': count,
            'combo': combo,
            'name': name,
            'description': description,
            'PFAMs': pfam
        })

    if not result_data:
        print(
            "No OG was found with a count > 1. "
            "Output file will not be created."
        )
        return

    result_df = pd.DataFrame(
        result_data,
        columns=[
            'OG',
            'count',
            'combo',
            'name',
            'description',
            'PFAMs'
        ]
    )

    result_df = result_df.sort_values(
        by='count',
        ascending=False
    )

    print(f"Writing output file: {output_file}")

    result_df.to_csv(
        output_file,
        sep='\t',
        index=False
    )

    print("Process complete.")

if __name__ == "__main__":
    process_upset_input(
        INPUT_FILE,
        ANNOTATION_DIR,
        OUTPUT_FILE
    )
