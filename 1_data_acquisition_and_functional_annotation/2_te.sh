#!/bin/bash


set -euo pipefail

# Example: ./76_fna/
INPUT_DIR="/home/yankh/76_converge/genomic_data/76_fna"
# Example output: ./TE_analysis/Aureobasidium_melanogenum_P16/
OUTPUT_DIR="/home/yankh/76_converge/genomic_data/TE_analysis"

THREADS_RM=10
THREADS_MASKER=5
PARALLEL_JOBS=5

mkdir -p "${OUTPUT_DIR}"

run_te_pipeline() {

    fna="$1"

    species=$(basename "${fna}" .fna)

    workdir="${OUTPUT_DIR}/${species}"

    mkdir -p "${workdir}"

    cd "${workdir}"

    echo "[$(date)] Start ${species}"

    BuildDatabase \
        -name "${species}.DB" \
        "${fna}"

    RepeatModeler \
        -database "${species}.DB" \
        -threads ${THREADS_RM}

    RepeatMasker \
        -lib "${species}.DB-families.fa" \
        -xsmall \
        -gff \
        -s \
        -pa ${THREADS_MASKER} \
        "${fna}"

    echo "[$(date)] Finished ${species}"
}

export -f run_te_pipeline
export THREADS_RM
export THREADS_MASKER
export OUTPUT_DIR

find "${INPUT_DIR}" -maxdepth 1 -name "*.fna" | \
parallel -j ${PARALLEL_JOBS} run_te_pipeline {}

echo "All jobs completed."
