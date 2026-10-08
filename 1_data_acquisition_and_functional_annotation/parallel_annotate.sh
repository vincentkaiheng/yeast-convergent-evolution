#!/bin/bash

# Example: ./76_faa/
input_dir="/home/yankh/76_converge/76_faa"
# Example: ./76_anno/
output_dir="/home/yankh/76_converge/76_anno"
dmnd_db="/share/ME-58T/yankh/eggnog_database/fungi.dmnd"
data_dir="/share/ME-58T/yankh/eggnog_database"
cpu_per_job=15
max_jobs=5

mkdir -p "$output_dir"

run_emapper() {
    local fa_file=$1
    local output_prefix="${fa_file##*/}"
    output_prefix="${output_prefix%.faa}"

    echo "Running command for $fa_file..."

    emapper.py -m diamond \
        --itype proteins \
        -i "$fa_file" \
        -o "$output_prefix" \
        --output_dir "$output_dir" \
        --cpu "$cpu_per_job" \
        --dbmem \
        --dmnd_db "$dmnd_db" \
        --data_dir "$data_dir"
}

export -f run_emapper
export input_dir
export output_dir
export dmnd_db
export data_dir
export cpu_per_job


find "$input_dir" -name "*.faa" | parallel --jobs "$max_jobs" run_emapper
