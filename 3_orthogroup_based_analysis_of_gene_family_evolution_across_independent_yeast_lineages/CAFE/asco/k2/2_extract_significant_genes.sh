#!/bin/bash

SEQUENCE_DIR="/home/yankh/76_converge/76_faa/OrthoFinder/Results_Jun15_1/Orthogroup_Sequences"

grep "y" Gamma_family_results.txt | cut -f1 > p0.05.significant

grep -f p0.05.significant Gamma_change.tab > Gamma_p0.05change.tab

while read node; do
    # Column 1 holds the family ID; node IDs follow the table column order.
    col=$((node + 1))

    contracted_ids_file="node${node}significant.contracted"
    output_dir="node${node}significant_contracted_genes"

    mkdir -p "$output_dir"

    # Negative changes are contractions; use positive changes for expansions.
    cut -f1,$col Gamma_p0.05change.tab | grep -P "\t-" | cut -f1 > "$contracted_ids_file"

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
