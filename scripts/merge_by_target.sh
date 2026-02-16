#!/bin/bash
set -euo pipefail

INPUT_DIR="/moosefs/guarracino/HPRCv2/aln-ba8599f934840b4b1a836e83fee938d9c0bc5e77"
OUTPUT_DIR="/scratch/hprcv2-25k"

mkdir -p "$OUTPUT_DIR"

# Get unique name2 values
find "$INPUT_DIR" -maxdepth 1 -name "*.paf.gz" -printf "%f\n" | \
    sed 's/.*-vs-//; s/\.p95\.Pinf\.aln\.paf\.gz$//' | \
    sort -u | \
while read -r name2; do
    output_file="${OUTPUT_DIR}/${name2}.merged.paf.gz"

    # Skip if already exists
    if [[ -f "$output_file" ]]; then
        echo "Skipping $name2 (already exists)"
        continue
    fi

    echo "Merging files for target: $name2"

    # Concatenate all gzipped files with this name2
    cat "$INPUT_DIR"/*-vs-"${name2}".p95.Pinf.aln.paf.gz > "$output_file"
done

echo "Done! Merged files in: $OUTPUT_DIR"
