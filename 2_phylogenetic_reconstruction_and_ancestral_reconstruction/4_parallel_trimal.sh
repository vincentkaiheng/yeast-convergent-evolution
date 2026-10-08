#!/bin/bash

input_folder="aligned_sequences"
output_folder="trimmed_sequences"

mkdir -p "$output_folder"

parallel_jobs=15

find "$input_folder" -type f -name '*.faa' | xargs -P $parallel_jobs -I {} bash -c '
    file="{}"
    base_name=$(basename "$file")

    output_file="'$output_folder'/$base_name"

    trimal -in "$file" -out "$output_file" -automated1

'

echo "trimAl processing complete."
