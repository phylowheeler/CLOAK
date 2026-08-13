#!/usr/bin/env bash

# Usage:
#   ./run_lrm_distance.sh GENES_DIR SPECIES_TREE OUTFILE [GLOB_PATTERN]
#
# Example:
#   ./run_lrm_distance.sh gene_trees species_tree.nwk lrm_out.txt "*.treefile"

GENES_DIR="${1:-}"
SPECIES_TREE="${2:-}"
OUTFILE="${3:-}"
GLOB_PATTERN="${4:-*.treefile}"

if [[ -z "$GENES_DIR" || -z "$SPECIES_TREE" || -z "$OUTFILE" ]]; then
  echo "Usage: $0 GENES_DIR SPECIES_TREE OUTFILE [GLOB_PATTERN]" >&2
  exit 1
fi

if [[ ! -d "$GENES_DIR" ]]; then
  echo "[error] GENES_DIR not a directory: $GENES_DIR" >&2
  exit 1
fi
if [[ ! -f "$SPECIES_TREE" ]]; then
  echo "[error] SPECIES_TREE not found: $SPECIES_TREE" >&2
  exit 1
fi

# How many cores to use in parallel (default to 4 if not under SLURM)
NPROCS="${SLURM_CPUS_PER_TASK:-4}"

module load anaconda
module load parallel          # make sure this module exists on HPC
source ~/.bashrc && conda activate lrm_distance

echo "[info] Using $NPROCS parallel jobs"
echo "[info] Writing output to: $OUTFILE"

# Truncate output file
: > "$OUTFILE"

shopt -s nullglob
files=( "$GENES_DIR"/$GLOB_PATTERN )
shopt -u nullglob

if [[ ${#files[@]} -eq 0 ]]; then
  echo "[warn] No files matched pattern '$GLOB_PATTERN' in '$GENES_DIR'." >&2
  exit 0
fi

# Export variables needed inside parallel
export SPECIES_TREE

# Run Python in parallel over the list of files.
# parallel preserves input order by default and we redirect all stdout to OUTFILE.
parallel -j "$NPROCS" '
  gene="{}"
  # compute distance; suppress stderr; if it fails, skip (or print ERROR)
  if dist=$(python3 calculate_lrm_distance.py --gene "$gene" --species "$SPECIES_TREE" 2>/dev/null); then
      printf "%s: %s\n" "$(basename "$gene")" "$dist"
  else
      # If you want to record failures, uncomment:
      # printf "%s: ERROR\n" "$(basename "$gene")"
      :
  fi
' ::: "${files[@]}" > "$OUTFILE"

echo "[done] Processed ${#files[@]} gene trees into: $OUTFILE"
