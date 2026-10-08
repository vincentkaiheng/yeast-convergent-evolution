#!/bin/bash
# Align proteins, trim alignments and estimate trees for paralog pruning.
set -o pipefail

INPUT_DIR="./OGs_13filtered"
ALIGNED_DIR="./aligned_OGs"
TRIMMED_DIR="./trimmed_OGs"
MAFFT_JOBS=30
TRIMAL_JOBS=15
FASTTREE_JOBS=10
LOG_DIR="./pipeline_logs"
MAFFT_OPTS="--localpair --maxiterate 1000"
TRIMAL_OPTS="-automated1"
FASTTREE_OPTS="-lg -gamma"

mkdir -p "$ALIGNED_DIR" "$TRIMMED_DIR" "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
MAIN_LOG="$LOG_DIR/pipeline_${TIMESTAMP}.log"
echo "Pipeline start: $(date -u +"%Y-%m-%d %H:%M:%S %Z")" | tee -a "$MAIN_LOG"

missing_tools=()
command -v mafft >/dev/null 2>&1 || missing_tools+=("mafft")
command -v trimal >/dev/null 2>&1 || missing_tools+=("trimal")
command -v FastTreeMP >/dev/null 2>&1 || missing_tools+=("FastTreeMP")

if [ ${#missing_tools[@]} -ne 0 ]; then
  echo "[ERROR] Missing tools: ${missing_tools[*]}" | tee -a "$MAIN_LOG"
  echo "Install these tools and add them to PATH." | tee -a "$MAIN_LOG"
  exit 1
fi

mapfile -t FA_FILES < <(find "$INPUT_DIR" -type f -name "*.fa" | sort)
NUM_INPUT=${#FA_FILES[@]}
echo "Detected ${NUM_INPUT} .fa files in ${INPUT_DIR}" | tee -a "$MAIN_LOG"
if [ "$NUM_INPUT" -eq 0 ]; then
  echo "[WARN] No input .fa files found." | tee -a "$MAIN_LOG"
  exit 0
fi

echo -e "\n=== Step 1: MAFFT align ===" | tee -a "$MAIN_LOG"
MAFFT_FAIL_LOG="$LOG_DIR/mafft_failures_${TIMESTAMP}.log"
: > "$MAFFT_FAIL_LOG"

find "$INPUT_DIR" -type f -name "*.fa" | sort | \
  xargs -n 1 -P "$MAFFT_JOBS" -I {} bash -c '
    infile="{}"
    filename=$(basename "$infile")
    outfile="'"$ALIGNED_DIR"'/aligned_${filename}"
    task_log="'"$LOG_DIR"'/mafft_${filename}.log"
    echo "[MAFFT] start ${filename} at $(date +%FT%T)" > "$task_log"
    if mafft '"$MAFFT_OPTS"' "$infile" > "$outfile" 2>>"$task_log"; then
      echo "[MAFFT] success ${filename} at $(date +%FT%T)" >> "$task_log"
    else
      echo "[MAFFT] FAIL ${filename} at $(date +%FT%T)" >> "$task_log"
      echo "${filename}" >> "'"$MAFFT_FAIL_LOG"'"
      [ -f "$outfile" ] && rm -f "$outfile"
    fi
'

echo "MAFFT finished. Failures (if any) listed in: $MAFFT_FAIL_LOG" | tee -a "$MAIN_LOG"
MAFFT_FAIL_COUNT=$(wc -l < "$MAFFT_FAIL_LOG" 2>/dev/null || echo 0)
echo "MAFFT failure count: ${MAFFT_FAIL_COUNT}" | tee -a "$MAIN_LOG"

echo -e "\n=== Step 2: trimal trim ===" | tee -a "$MAIN_LOG"
TRIMAL_FAIL_LOG="$LOG_DIR/trimal_failures_${TIMESTAMP}.log"
: > "$TRIMAL_FAIL_LOG"

find "$ALIGNED_DIR" -type f -name "aligned_*.fa" | sort | \
  xargs -n 1 -P "$TRIMAL_JOBS" -I {} bash -c '
    infile="{}"
    base=$(basename "$infile")
    outfile="'"$TRIMMED_DIR"'/${base}"
    task_log="'"$LOG_DIR"'/trimal_${base}.log"
    echo "[TRIMAL] start ${base} at $(date +%FT%T)" > "$task_log"
    if trimal -in "$infile" -out "$outfile" '"$TRIMAL_OPTS"' >> "$task_log" 2>&1; then
      echo "[TRIMAL] success ${base} at $(date +%FT%T)" >> "$task_log"
    else
      echo "[TRIMAL] FAIL ${base} at $(date +%FT%T)" >> "$task_log"
      echo "${base}" >> "'"$TRIMAL_FAIL_LOG"'"
      [ -f "$outfile" ] && rm -f "$outfile"
    fi
'

echo "trimal finished. Failures (if any) listed in: $TRIMAL_FAIL_LOG" | tee -a "$MAIN_LOG"
TRIMAL_FAIL_COUNT=$(wc -l < "$TRIMAL_FAIL_LOG" 2>/dev/null || echo 0)
echo "trimal failure count: ${TRIMAL_FAIL_COUNT}" | tee -a "$MAIN_LOG"

echo -e "\n=== Step 3: FastTreeMP build trees ===" | tee -a "$MAIN_LOG"
FASTTREE_FAIL_LOG="$LOG_DIR/fasttree_failures_${TIMESTAMP}.log"
: > "$FASTTREE_FAIL_LOG"

find "$TRIMMED_DIR" -type f -name "aligned_*.fa" | sort | \
  xargs -n 1 -P "$FASTTREE_JOBS" -I {} bash -c '
    infile="{}"
    base=$(basename "$infile" .fa)
    treefile="'"$TRIMMED_DIR"'/${base}.tre"
    task_log="'"$LOG_DIR"'/fasttree_${base}.log"
    echo "[FastTree] start ${base} at $(date +%FT%T)" > "$task_log"
    if FastTreeMP '"$FASTTREE_OPTS"' "$infile" > "$treefile" 2>>"$task_log"; then
      echo "[FastTree] success ${base} at $(date +%FT%T)" >> "$task_log"
    else
      echo "[FastTree] FAIL ${base} at $(date +%FT%T)" >> "$task_log"
      echo "${base}" >> "'"$FASTTREE_FAIL_LOG"'"
      [ -f "$treefile" ] && rm -f "$treefile"
    fi
'

echo "FastTreeMP finished. Failures (if any) listed in: $FASTTREE_FAIL_LOG" | tee -a "$MAIN_LOG"
FASTTREE_FAIL_COUNT=$(wc -l < "$FASTTREE_FAIL_LOG" 2>/dev/null || echo 0)
echo "FastTreeMP failure count: ${FASTTREE_FAIL_COUNT}" | tee -a "$MAIN_LOG"

echo -e "\n=== Summary ===" | tee -a "$MAIN_LOG"
echo "Input fasta count: ${NUM_INPUT}" | tee -a "$MAIN_LOG"
echo "MAFFT failures: ${MAFFT_FAIL_COUNT}" | tee -a "$MAIN_LOG"
echo "trimal failures: ${TRIMAL_FAIL_COUNT}" | tee -a "$MAIN_LOG"
echo "FastTree failures: ${FASTTREE_FAIL_COUNT}" | tee -a "$MAIN_LOG"
echo "Detailed logs are in: $LOG_DIR" | tee -a "$MAIN_LOG"
echo "Pipeline end: $(date -u +"%Y-%m-%d %H:%M:%S %Z")" | tee -a "$MAIN_LOG"

total_failures=$((MAFFT_FAIL_COUNT + TRIMAL_FAIL_COUNT + FASTTREE_FAIL_COUNT))
if [ "$total_failures" -gt 0 ]; then
  echo "[WARN] pipeline finished with ${total_failures} failed items." | tee -a "$MAIN_LOG"
  exit 2
else
  echo "[OK] pipeline finished without per-file failures." | tee -a "$MAIN_LOG"
  exit 0
fi
