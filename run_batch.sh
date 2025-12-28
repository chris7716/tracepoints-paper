#!/usr/bin/env bash

# Simple script to run the report generation commands
# Usage: ./run_batch.sh <start_N> <type> [complexity_metric] [threads]

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <start_N> <type> [complexity_metric] [threads]"
    echo ""
    echo "Parameters:"
    echo "  start_N              Starting N ID number"
    echo "  type                 Tracepoint type (standard, mixed, variable, fastga)"
    echo "  complexity_metric    Optional: complexity metric (default: diagonal-distance)"
    echo "  threads              Optional: number of threads (default: 8)"
    echo ""
    echo "Note: Max complexity values are defined in the configs array within the script"
    echo ""
    echo "Examples:"
    echo "  $0 61 standard"
    echo "  $0 61 standard edit-distance"
    echo "  $0 61 standard diagonal-distance 16"
    exit 1
fi

START_N=$1
TYPE=$2
COMPLEXITY_METRIC=${3:-diagonal-distance}  # Default to diagonal-distance
THREADS=${4:-8}                            # Default to 8

# Validate inputs
if ! [[ "$START_N" =~ ^[0-9]+$ ]]; then
    echo "Error: start_N must be a positive integer, got: '$START_N'" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: threads must be a positive integer, got: '$THREADS'" >&2
    exit 1
fi

echo "Batch run configuration:"
echo "  Type: $TYPE"
echo "  Complexity Metric: $COMPLEXITY_METRIC"
echo "  Threads: $THREADS"
echo "  Starting N: $START_N"
echo ""

# Array of configurations: "records error length max_complexity"
configs=( 
    "10000 0.1 10000 10"
    "10000 0.1 10000 25"
    "10000 0.1 10000 32"
    "10000 0.1 10000 50"
    "10000 0.1 10000 100"
    "10000 0.1 10000 200"
)

N=$START_N

for config in "${configs[@]}"; do
    read -r records error length max_complexity <<< "$config"
    
    # Format error for filename (remove decimal point)
    error_str=$(echo "$error" | sed 's/0\.//' | sed 's/\.//')
    
    # Create more descriptive filename that includes complexity settings
    output_file="reports/decode-runtime-${TYPE}-${COMPLEXITY_METRIC}-mc${max_complexity}-R${records}-L${length}-E${error_str}-N${N}.md"
    
    echo "Running: N=$N, R=$records, L=$length, E=$error, MC=$max_complexity (${COMPLEXITY_METRIC}, t=${THREADS}) -> $output_file"
    
    # Run with all parameters
    if ! ./run_steps_report.sh -r "$records" -E "$error" -L "$length" -N "$N" -T "$TYPE" -C "$COMPLEXITY_METRIC" -M "$max_complexity" -t "$THREADS" > "$output_file"; then
        echo "Error: Failed to generate report for N=$N" >&2
        echo "Continuing with next configuration..."
        continue
    fi
    
    echo "Completed: $output_file"
    
    ((N++))
done

echo ""
echo "Batch run summary:"
echo "  Total configurations: ${#configs[@]}"
echo "  Type: $TYPE"
echo "  Complexity Metric: $COMPLEXITY_METRIC"
echo "  Threads: $THREADS"
echo "  N range: $START_N - $((N-1))"
echo "All runs completed!"