#!/bin/bash

input_directory="./filtered_single_copy_sequences"
output_directory="./aligned_sequences"

mkdir -p "$output_directory"

parallel_jobs=60

find "$input_directory" -type f -name "*.faa" | xargs -n 1 -P $parallel_jobs -I {} bash -c '
  input_file="{}"
  filename=$(basename "$input_file")

  output_file="'$output_directory'/aligned_$filename"

  mafft --localpair --maxiterate 1000 "$input_file" > "$output_file"
'

echo "MAFFT alignments complete."
