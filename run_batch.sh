#!/usr/bin/env bash

# Simple script to run the report generation commands
# Usage: ./run_batch.sh <start_N> <type>

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <start_N> <type>"
    echo "Example: $0 61 standard"
    exit 1
fi

START_N=$1
TYPE=$2

# Array of configurations: "records error length"
configs=(
    "10000 0.01 100"
    "10000 0.01 1000" 
    "10000 0.01 10000"
    "10000 0.05 100"
    "10000 0.05 1000"
    "10000 0.05 10000"
    "10000 0.1 100"
    "10000 0.1 1000"
    "10000 0.1 10000"
    "10000 0.2 100"
    "10000 0.2 1000"
    "10000 0.2 10000"
)

N=$START_N

for config in "${configs[@]}"; do
    read -r records error length <<< "$config"
    
    # Format error for filename (remove decimal point)
    error_str=$(echo "$error" | sed 's/0\.//' | sed 's/\.//')
    
    output_file="reports/decode-runtime-${TYPE}-R${records}-L${length}-E${error_str}-N${N}.md"
    
    echo "Running: N=$N, R=$records, L=$length, E=$error -> $output_file"
    
    ./run_steps_report.sh -r "$records" -E "$error" -L "$length" -N "$N" -T "$TYPE" > "$output_file"
    
    echo "Completed: $output_file"
    
    ((N++))
done

echo "All runs completed!"
