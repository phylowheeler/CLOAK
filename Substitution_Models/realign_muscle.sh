#!/bin/bash


set -euo pipefail

###############################################################################
# One-node MUSCLE5 realignment driver using GNU Parallel
#
# Design:
#   - 47 concurrent jobs
#   - 2 CPU threads per MUSCLE job
#   - total CPU allocation = 94
#
# Behavior on failures:
#   - DOES NOT stop on first failed locus
#   - Continues through all files
#   - Writes failed job details and failed input list under logs/
#
# Usage:
#   sbatch realign_muscle_parallel.sh INPUT_DIR OUTPUT_DIR MUSCLE_EXE PY_SCRIPT
###############################################################################

if [[ $# -lt 4 ]]; then
    echo "Usage: sbatch $0 INPUT_DIR OUTPUT_DIR MUSCLE_EXE PY_SCRIPT" >&2
    exit 1
fi

INPUT_DIR="$(realpath "$1")"
OUTPUT_DIR="$(realpath -m "$2")"
MUSCLE_EXE="$(realpath "$3")"
PY_SCRIPT="$(realpath "$4")"

JOBS=47
THREADS_PER_JOB=2

# Keep nested libraries from oversubscribing cores.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

LOG_DIR="$OUTPUT_DIR/logs"
SUMMARY_DIR="$OUTPUT_DIR/summaries"
ALIGNED_DIR="$OUTPUT_DIR/aligned"
mkdir -p "$LOG_DIR" "$SUMMARY_DIR" "$ALIGNED_DIR"

module load anaconda
module load parallel
module load python/3.9/3.9.10  # (OK if redundant with conda)

# Initialize conda and activate your env (as you had it)
source "$(conda info --base)"/etc/profile.d/conda.sh
conda activate muscledivver

JOBLOG="$LOG_DIR/parallel_joblog.tsv"
FAILED_JOBLOG="$LOG_DIR/failed_jobs.tsv"
FAILED_INPUTS="$LOG_DIR/failed_inputs.txt"
PARALLEL_RESULTS_DIR="$LOG_DIR/parallel_results"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: INPUT_DIR does not exist: $INPUT_DIR" >&2
    exit 1
fi

if [[ ! -x "$MUSCLE_EXE" ]]; then
    echo "ERROR: MUSCLE executable is not executable: $MUSCLE_EXE" >&2
    exit 1
fi

if [[ ! -f "$PY_SCRIPT" ]]; then
    echo "ERROR: Python script not found: $PY_SCRIPT" >&2
    exit 1
fi

run_one() {
    local infile="$1"
    local base stem outfile summary

    base="$(basename "$infile")"
    stem="${base%.*}"
    outfile="$ALIGNED_DIR/${stem}_realigned.fasta"
    summary="$SUMMARY_DIR/${stem}.summary.tsv"

    srun --exclusive -N1 -n1 -c "$THREADS_PER_JOB" \
        python3 "$PY_SCRIPT" \
            --input "$infile" \
            --output "$outfile" \
            --muscle "$MUSCLE_EXE" \
            --threads "$THREADS_PER_JOB" \
            --summary "$summary"
}
export -f run_one
export PY_SCRIPT MUSCLE_EXE THREADS_PER_JOB ALIGNED_DIR SUMMARY_DIR

declare -a files=()
while IFS= read -r -d '' f; do
    files+=("$f")
done < <(
    find "$INPUT_DIR" -maxdepth 1 -type f \
        \( -iname '*.fa' -o -iname '*.faa' -o -iname '*.fasta' -o -iname '*.fas' -o -iname '*.afa' \) \
        -print0 | sort -z
)

if [[ ${#files[@]} -eq 0 ]]; then
    echo "ERROR: No FASTA-like files found in $INPUT_DIR" >&2
    exit 1
fi

echo "Input directory   : $INPUT_DIR"
echo "Output directory  : $OUTPUT_DIR"
echo "MUSCLE executable : $MUSCLE_EXE"
echo "Python script     : $PY_SCRIPT"
echo "Files found       : ${#files[@]}"
echo "Parallel jobs      = $JOBS"
echo "Threads per MUSCLE = $THREADS_PER_JOB"
echo

set +e
printf '%s\0' "${files[@]}" | \
    parallel -0 \
             --jobs "$JOBS" \
             --joblog "$JOBLOG" \
             --results "$PARALLEL_RESULTS_DIR" \
             run_one {}
parallel_exit=$?
set -e

if [[ -f "$JOBLOG" ]]; then
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($7 != 0 || $8 != 0){print}' "$JOBLOG" > "$FAILED_JOBLOG"
    awk 'BEGIN{FS=OFS="\t"} NR>1 && ($7 != 0 || $8 != 0){sub(/^run_one /, "", $9); print $9}' "$JOBLOG" > "$FAILED_INPUTS"
else
    : > "$FAILED_JOBLOG"
    : > "$FAILED_INPUTS"
fi

failed_count=0
if [[ -s "$FAILED_INPUTS" ]]; then
    failed_count=$(wc -l < "$FAILED_INPUTS")
fi

success_count=$(( ${#files[@]} - failed_count ))

echo
echo "Finished."
echo "Realigned FASTA files : $ALIGNED_DIR"
echo "Per-file summaries    : $SUMMARY_DIR"
echo "GNU Parallel logs     : $LOG_DIR"
echo "Successful loci       : $success_count"
echo "Failed loci           : $failed_count"

if [[ $failed_count -gt 0 ]]; then
    echo
    echo "Some loci failed, but the job continued through the full directory."
    echo "Failed input list : $FAILED_INPUTS"
    echo "Failed job log    : $FAILED_JOBLOG"
    echo "GNU Parallel exit : $parallel_exit"
fi

exit 0
