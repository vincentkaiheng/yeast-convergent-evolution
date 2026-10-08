#!/bin/bash

SEQUENCE_DIR="/home/yankh/76_converge/76_faa/OrthoFinder/Results_Jun15_1/Orthogroup_Sequences"

while read node; do
    col=$((node + 1))

    contracted_ids_file="node${node}.contracted"
    output_dir="node${node}_contracted_genes"

    mkdir -p "$output_dir"

    awk -v c="$col" 'NR > 1 && $c < 0 {print $1}' Gamma_change.tab > "$contracted_ids_file"

    echo "Node $node: $(wc -l < "$contracted_ids_file") contracted gene families"

    while read og; do
        fasta_file="${SEQUENCE_DIR}/${og}.fa"
        if [ -f "$fasta_file" ]; then
            cp "$fasta_file" "$output_dir/"
        else
            echo "Warning: $fasta_file not found for node $node"
        fi
    done < "$contracted_ids_file"

done < node.txt
