#!/usr/bin/env bash

#script for sampling sequence files from a directory
# Check for required arguments: input a source directory to sample from, an output directory, and the desired sample size
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <source_dir> <dest_dir> <sample_size>"
    exit 1
fi

SOURCE_DIR="$1"
DEST_DIR="$2"
SAMPLE_SIZE="$3"

# Validate source directory
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory '$SOURCE_DIR' does not exist."
    exit 1
fi

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Sample files safely without replacement
find "$SOURCE_DIR" -maxdepth 1 -type f -print0 \
    | shuf -z -n "$SAMPLE_SIZE" \
    | xargs -0 -I {} cp {} "$DEST_DIR/"

echo "Successfully sampled $SAMPLE_SIZE files into '$DEST_DIR'."