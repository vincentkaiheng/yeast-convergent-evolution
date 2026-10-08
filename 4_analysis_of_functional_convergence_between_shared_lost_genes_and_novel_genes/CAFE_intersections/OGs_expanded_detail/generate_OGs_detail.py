#!/usr/bin/env python3

import os
import pandas as pd
from collections import Counter

INPUT_FILE = '../OGs_expanded_upset_outputs/og_upset_input.tsv'
ANNOTATION_DIR = '/home/yankh/76_converge/76_faa/all_OGs_anno'
OUTPUT_FILE = './OG_detail.tsv'

def get_name_and_description_from_annotation(og_id, annotation_dir):
    """Return the modal name and description from columns 9 and 8."""

    annotation_file = os.path.join(annotation_dir, f"{og_id}.emapper.annotations")

    if not os.path.exists(annotation_file):
        return 'No_annotation', 'No_description'

    names = []
    descriptions = []

    try:
        with open(annotation_file, 'r', encoding='utf-8') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue

                columns = line.rstrip('\n').split('\t')

                if len(columns) > 7 and columns[7].strip():
                    descriptions.append(columns[7].strip())
                if len(columns) > 8 and columns[8].strip():
                    names.append(columns[8].strip())
    except Exception as e:
        print(f"Warning: Could not read annotation file {annotation_file}. Error: {e}")
        return 'Read_Error', 'Read_Error'

    most_common_name = Counter(names).most_common(1)[0][0] if names else 'No_annotation'
    most_common_description = Counter(descriptions).most_common(1)[0][0] if descriptions else 'No_description'

    return most_common_name, most_common_description

def process_upset_input(input_file, annotation_dir, output_file):
    """Annotate OGs shared by at least two lineages."""

    print(f"Reading input file: {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t', dtype=str)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file}")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    series_cols = df.columns[1:].tolist()

    result_data = []

    for index, row in df.iterrows():
        og_id = row['OG']
        present_series = []

        for col in series_cols:
            if str(row[col]) == '1':
                present_series.append(col)

        count = len(present_series)
        combo = '&'.join(present_series)

        if count <= 1:
            continue

        name, description = get_name_and_description_from_annotation(og_id, annotation_dir)

        result_data.append({
            'OG': og_id,
            'count': count,
            'combo': combo,
            'name': name,
            'description': description
        })

    if not result_data:
        print("No OG was found with a count > 1. Output file will be empty or not created.")
        return

    result_df = pd.DataFrame(result_data, columns=['OG', 'count', 'combo', 'name', 'description'])

    result_df = result_df.sort_values(by='count', ascending=False)

    print(f"Writing output file: {output_file}")
    result_df.to_csv(output_file, sep='\t', index=False)
    print("Process complete.")

if __name__ == "__main__":
    process_upset_input(INPUT_FILE, ANNOTATION_DIR, OUTPUT_FILE)
