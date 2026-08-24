#!/bin/bash

# ========= CONFIGURATION =========
ALIGN_DIR=
OUTPUT_DIR=
MODEL_FILE=
THREADS_PER_JOB=1             

mkdir -p "$OUTPUT_DIR"
module load anaconda
source ~/.bashrc && conda activate
conda activate iqtree2

export ALIGN_DIR OUTPUT_DIR THREADS_PER_JOB MODEL_FILE

run_iqtree() {
    local ALIGN_FILE=$1
    OUT_PREFIX="${OUTPUT_DIR}/$(basename "$ALIGN_FILE" ".${ALIGN_FILE##*.}")"

    # Check if alignment file exists
    if [[ ! -f "$ALIGN_FILE" ]]; then
        echo "[$(date)] WARNING: Alignment file not found for $ALIGN_FILE" >> "$OUTPUT_DIR/missing_files.log"
        return
    fi

    # Skip if output already exists
    if [[ -f "${OUT_PREFIX}.treefile" ]]; then
        echo "[$(date)] Skipping $ALIGN_FILE (already completed)"
        return
    fi

    # Build model 
    MODEL_emperical="${MODEL_FILE}"

    echo "[$(date)] Starting $ALIGN_FILE with model ${MODEL_emperical}"

    iqtree3 -s "$ALIGN_FILE" \
            -m "$MODEL_emperical" \
            -mfreq F \
            -mrate E,I,G,I+G,R \
            -T "$THREADS_PER_JOB" \
            -seed 42 \
            -keep-ident \
            -redo \
            -quiet \
            -pre "$OUT_PREFIX" \
        || echo "[$(date)] ERROR: $ALIGN_FILE failed" >> "$OUTPUT_DIR/failed_jobs.log"

    # Clean up intermediate files (only if successful)
    if [[ -f "${OUT_PREFIX}.treefile" ]]; then
        rm -f "${OUT_PREFIX}.best_model.nex" \
              "${OUT_PREFIX}.bionj" \
              "${OUT_PREFIX}.ckp.gz" \
              "${OUT_PREFIX}.mldist"
        echo "[$(date)] Finished $ALIGN_FILE"
    fi
}

export -f run_iqtree

for file in "$ALIGN_DIR"/*; 
do [ -f "$file" ] && 
run_iqtree "$file" ; done

echo "[$(date)] All jobs complete."
