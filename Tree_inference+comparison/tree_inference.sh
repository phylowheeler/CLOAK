#!/bin/bash

# ========= CONFIGURATION =========
ALIGN_DIR=
OUTPUT_DIR=
MODEL_FILE=
THREADS_PER_JOB=1             

mkdir -p "$OUTPUT_DIR"
module load anaconda
module load parallel
source ~/.bashrc && conda activate
conda activate iqtree3

export ALIGN_DIR OUTPUT_DIR THREADS_PER_JOB MODEL_FILE

run_iqtree() {
    ID="$1"
    ALIGN_FILE="${ALIGN_DIR}/${ID}.nex"
    OUT_PREFIX="${OUTPUT_DIR}/${ID}"

    # Check if alignment file exists
    if [[ ! -f "$ALIGN_FILE" ]]; then
        echo "[$(date)] WARNING: Alignment file not found for $ID" >> "$OUTPUT_DIR/missing_files.log"
        return
    fi

    # Skip if output already exists
    if [[ -f "${OUT_PREFIX}.treefile" ]]; then
        echo "[$(date)] Skipping $ID (already completed)"
        return
    fi

    # Build model 
    MODEL_emperical="${MODEL_FILE}"

    echo "[$(date)] Starting $ID with model ${MODEL_emperical}"

    iqtree3 -s "$ALIGN_FILE" \
            -m "$MODEL_emperical" \
            -mfreq F
            -mrate E,I,G,I+G,R
            -T "$THREADS_PER_JOB" \
            -seed 42 \
            -keep-ident \
            -redo \
            -quiet \
            -pre "$OUT_PREFIX" \
        || echo "[$(date)] ERROR: $ID failed" >> "$OUTPUT_DIR/failed_jobs.log"

    # Clean up intermediate files (only if successful)
    if [[ -f "${OUT_PREFIX}.treefile" ]]; then
        rm -f "${OUT_PREFIX}.best_model.nex" \
              "${OUT_PREFIX}.bionj" \
              "${OUT_PREFIX}.ckp.gz" \
              "${OUT_PREFIX}.mldist"
        echo "[$(date)] Finished $ID"
    fi
}

export -f run_iqtree

# Number of jobs to run in parallel = total cpus / threads per job
MAX_JOBS=$((90 / THREADS_PER_JOB))

echo "[$(date)] Starting batch with $MAX_JOBS parallel IQ-TREE runs..."
parallel --jobs "$MAX_JOBS" run_iqtree
echo "[$(date)] All jobs complete."
