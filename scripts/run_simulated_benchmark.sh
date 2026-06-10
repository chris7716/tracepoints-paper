#!/bin/bash
set -euo pipefail

# =============================================================================
# Simulated Data Benchmarking Pipeline
# =============================================================================
# This script runs the complete simulated data generation and benchmarking
# pipeline from code.md. It generates sequences, aligns them, and benchmarks
# tracepoint encoding/decoding across multiple parameter combinations.
#
# Usage:
#   ./run_simulated_benchmark.sh [--generate] [--benchmark] [--plot]
#
# Options:
#   --generate   Run sequence generation and alignment (Step 1-2)
#   --benchmark  Run the full benchmark pipeline (Step 3)
#   --plot       Generate plots from results (Step 4)
#   (no args)    Run all steps
#
# Before running, configure the paths in the "CONFIGURATION" section below.
# =============================================================================

# =============================================================================
# CONFIGURATION - Adjust these paths to your environment
# =============================================================================

# Base directories
dir_base="${DIR_BASE:-/moosefs/guarracino/tracepoints}"
scratch_dir="${SCRATCH_DIR:-/scratch}"

# Tool binaries - adjust paths as needed
generate_dataset="${GENERATE_DATASET:-/moosefs/guarracino/git/WFA2-lib/bin/generate_dataset}"
align_benchmark="${ALIGN_BENCHMARK:-/moosefs/guarracino/git/WFA2-lib/bin/align_benchmark}"
out_to_paf="${OUT_TO_PAF:-$dir_base/scripts/out-to-paf.py}"
pafcheck="${PAFCHECK:-/moosefs/guarracino/git/pafcheck/target/release/pafcheck}"
cigzip="${CIGZIP:-$scratch_dir/cigzip/target/release/cigzip}"

# Benchmark parameters
NUM_RECORDS="${NUM_RECORDS:-10000}"
THREADS="${THREADS:-$(nproc)}"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

check_tool() {
    local tool="$1"
    local path="$2"
    # For Python scripts, just check if file exists and is readable
    if [[ "$path" == *.py ]]; then
        if [[ ! -r "$path" ]]; then
            echo "ERROR: $tool not found at: $path"
            echo "Please install $tool or update the path in this script."
            exit 1
        fi
    else
        if [[ ! -x "$path" ]]; then
            echo "ERROR: $tool not found or not executable at: $path"
            echo "Please install $tool or update the path in this script."
            exit 1
        fi
    fi
}

parse_time_log() {
    local log="$1"
    runtime_raw=$(grep "Elapsed (wall clock)" "$log" | sed 's/.*: //')
    echo "$runtime_raw" | awk -F: '{
        if (NF == 3) print $1 * 3600 + $2 * 60 + $3
        else if (NF == 2) print $1 * 60 + $2
        else print $1
    }'
}

parse_memory_log() {
    grep "Maximum resident set size" "$1" | sed 's/.*: //'
}

# =============================================================================
# STEP 1 & 2: GENERATE SEQUENCES AND ALIGNMENTS
# =============================================================================

generate_data() {
    log "Starting data generation..."

    # Check required tools
    check_tool "generate_dataset" "$generate_dataset"
    check_tool "align_benchmark" "$align_benchmark"
    check_tool "out-to-paf.py" "$out_to_paf"
    check_tool "pafcheck" "$pafcheck"

    # Create directories
    mkdir -p "$dir_base/simulated-data/alns"
    mkdir -p "$dir_base/simulated-data/seqs"

    for l in 100 1000 10000 100000; do
        for e in 0.001; do
            log "Generating: length=$l, error=$e"

            l_nodot=$(echo $l | sed 's/\.//g')
            e_nodot=$(echo $e | sed 's/\.//g')
            prefix=set_${l_nodot}_${e_nodot}

            # Generate sequences
            "$generate_dataset" -n "$NUM_RECORDS" -l "$l" -e "$e" \
                -o "$dir_base/simulated-data/seqs/$prefix.seq"

            # Align (single thread to preserve order)
            "$align_benchmark" \
                -i "$dir_base/simulated-data/seqs/$prefix.seq" \
                -a edit-wfa --wfa-memory high \
                --output "$dir_base/simulated-data/alns/$prefix.out" \
                -t 1

            # Convert SEQ to FASTA
            # WFA2 generate_dataset writes > = seqshort (PATTERN/query) and < = seqlong (TEXT/target)
            awk '
              BEGIN { counter=1 }
              /^>/ {
                  print ">query_" counter
                  print substr($0, 2)
              }
              /^</ {
                  print ">target_" counter
                  print substr($0, 2)
                  counter++
              }' "$dir_base/simulated-data/seqs/$prefix.seq" \
                > "$dir_base/simulated-data/seqs/$prefix.fa"

            # Convert OUT to PAF
            python3 "$out_to_paf" \
                "$dir_base/simulated-data/alns/$prefix.out" \
                "$dir_base/simulated-data/alns/$prefix.paf"

            # Verify
            "$pafcheck" \
                --paf "$dir_base/simulated-data/alns/$prefix.paf" \
                --sequence-files "$dir_base/simulated-data/seqs/$prefix.fa" \
                --threads "$THREADS"

            log "Completed: $prefix"
        done
    done

    log "Data generation complete!"
}

# =============================================================================
# STEP 3: RUN BENCHMARK PIPELINE
# =============================================================================

run_benchmark() {
    log "Starting benchmark pipeline..."

    check_tool "cigzip" "$cigzip"

    # Create output directories
    mkdir -p "$dir_base/simulated-data"/{encode,compress,decompress,decode}
    mkdir -p "$scratch_dir/benchmark_tmp"

    REPORT="$dir_base/simulated-data/benchmark.results.fastga-no-diff.tsv"
    TMP_DIR="$scratch_dir/benchmark_tmp"

    # -------------------------------------------------------------------------
    # Precompute CIGAR sizes and baseline alignment times
    # -------------------------------------------------------------------------
    log "Precomputing CIGAR sizes and baseline times..."

    > "$TMP_DIR/cigar_sizes.txt"

    for l in 100 1000 10000 100000; do
        for e in 0.001; do
            l_nodot=$(echo $l | sed 's/\.//g')
            e_nodot=$(echo $e | sed 's/\.//g')
            prefix=set_${l_nodot}_${e_nodot}
            input_paf="$dir_base/simulated-data/alns/$prefix.paf"

            log "Processing baseline: $prefix"

            # Compute CIGAR size (12 cols + cg:Z: only)
            size_cigar=$(awk -F'\t' '{
                out = $1
                for (i = 2; i <= 12; i++) out = out "\t" $i
                for (i = 13; i <= NF; i++) {
                    if ($i ~ /^cg:Z:/) { out = out "\t" $i; break }
                }
                print out
            }' "$input_paf" | wc -c)

            if [ "$size_cigar" -le 1 ]; then
                echo "ERROR: cg:Z: tag not found in $input_paf" >&2
                exit 1
            fi

            # Compress CIGAR-only PAF with bgzip
            bgzip_file="$dir_base/simulated-data/compress/$prefix.cigar.paf.gz"
            awk -F'\t' '{
                out = $1
                for (i = 2; i <= 12; i++) out = out "\t" $i
                for (i = 13; i <= NF; i++) {
                    if ($i ~ /^cg:Z:/) { out = out "\t" $i; break }
                }
                print out
            }' "$input_paf" | bgzip -l9 -c > "$bgzip_file"
            size_cigar_bgzip=$(wc -c < "$bgzip_file")

            # Benchmark bgzip decompression
            bgzip_decomp_log="$dir_base/simulated-data/decompress/$prefix.bgzip_decompress.log"
            \time -v bgzip -d -@ "$THREADS" -c "$bgzip_file" > "$scratch_dir/${prefix}.cigar.paf" 2> "$bgzip_decomp_log"
            rm -f "$scratch_dir/${prefix}.cigar.paf"

            bgzip_decomp_runtime=$(parse_time_log "$bgzip_decomp_log")
            bgzip_decomp_memory=$(parse_memory_log "$bgzip_decomp_log")

            # Compute baseline alignment time using cigzip with single tracepoint per alignment
            baseline_tp_paf="$TMP_DIR/$prefix.baseline.tp.paf"
            baseline_decode_paf="$TMP_DIR/$prefix.baseline.decoded.paf"
            baseline_log="$TMP_DIR/$prefix.baseline.log"

            # Encode with very high max-complexity = 1 tracepoint per alignment
            "$cigzip" encode \
                --paf "$input_paf" \
                --type standard \
                --complexity-metric edit-distance \
                --max-complexity 9999999 \
                > "$baseline_tp_paf"

            # Decode and measure
            \time -v "$cigzip" decode \
                --paf "$baseline_tp_paf" \
                --sequence-files "$dir_base/simulated-data/seqs/$prefix.fa" \
                --type standard \
                --complexity-metric edit-distance \
                --max-complexity 9999999 \
                --distance edit \
                --memory-mode high \
                -t "$THREADS" \
                > "$baseline_decode_paf" 2> "$baseline_log"

            align_runtime=$(parse_time_log "$baseline_log")
            align_memory=$(parse_memory_log "$baseline_log")

            mv "$baseline_log" "$dir_base/simulated-data/alns/$prefix.baseline.log"
            rm -f "$baseline_tp_paf" "$baseline_decode_paf"

            echo "${l}_${e} $size_cigar $size_cigar_bgzip $bgzip_decomp_runtime $bgzip_decomp_memory $align_runtime $align_memory" >> "$TMP_DIR/cigar_sizes.txt"
        done
    done

    # -------------------------------------------------------------------------
    # Run benchmarks for all parameter combinations
    # -------------------------------------------------------------------------
    log "Running benchmark combinations..."

    # Generate jobs
    > "$TMP_DIR/jobs.txt"
    for l in 100 1000 10000 100000; do
        for e in 0.001; do
            for tp_type in fastga-no-diff; do
                if [ "$tp_type" = "fastga" ] || [ "$tp_type" = "fastga-no-diff" ]; then
                    cm_list="none"
                    mc_list="100"
                else
                    cm_list="edit-distance diagonal-distance"
                    mc_list="32 64 128 256 512 1024"
                fi

                for cm in $cm_list; do
                    for mc in $mc_list; do
                        for memory_mode in high; do
                            echo "$l $e $tp_type $cm $mc $memory_mode" >> "$TMP_DIR/jobs.txt"
                        done
                    done
                done
            done
        done
    done

    mkdir -p "$scratch_dir/jobs"

    # Process each job
    while read -r l e tp_type cm mc memory_mode; do
        l_nodot=$(echo $l | sed 's/\.//g')
        e_nodot=$(echo $e | sed 's/\.//g')
        prefix=set_${l_nodot}_${e_nodot}
        input_paf="$dir_base/simulated-data/alns/$prefix.paf"
        seq_file="$dir_base/simulated-data/seqs/$prefix.fa"

        # Look up precomputed values
        size_cigar=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f2)
        size_cigar_bgzip=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f3)
        bgzip_decomp_runtime=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f4)
        bgzip_decomp_memory=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f5)
        align_runtime=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f6)
        align_memory=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f7)

        prefix2=$tp_type.$cm.$mc.$memory_mode
        full_prefix="$prefix.$prefix2"

        log "Benchmarking: $full_prefix"

        [ "$tp_type" = "fastga" ] || [ "$tp_type" = "fastga-no-diff" ] && cmd_args="" || cmd_args="--complexity-metric $cm"
        if [ "$tp_type" = "fastga" ]; then strategy_args="--strategy rice;rice"
        elif [ "$tp_type" = "fastga-no-diff" ]; then strategy_args="--strategy rice"
        else strategy_args=""
        fi

        # Determine number of repeats
        if [ "$mc" -le 64 ] || [ "$tp_type" = "fastga" ] || [ "$tp_type" = "fastga-no-diff" ]; then
            num_repeats=10
        else
            num_repeats=1
        fi

        job_scratch="$scratch_dir/jobs/$full_prefix"
        mkdir -p "$job_scratch"

        # File paths
        tp_paf_full="$job_scratch/$full_prefix.tp.full.paf"
        tp_paf="$job_scratch/$full_prefix.tp.paf"
        tpa_file="$job_scratch/$full_prefix.tpa"
        decomp_paf="$job_scratch/$full_prefix.decomp.tp.paf"
        decode_paf="$job_scratch/$full_prefix.paf"
        decode_heur_paf="$job_scratch/$full_prefix.heuristic.paf"
        encode_log="$job_scratch/$full_prefix.encode.log"
        compress_log="$job_scratch/$full_prefix.compress.log"
        decompress_log="$job_scratch/$full_prefix.decompress.log"
        decode_log="$job_scratch/$full_prefix.decode.log"
        decode_heur_log="$job_scratch/$full_prefix.decode.heuristic.log"

        # AWK script for tp:Z: extraction
        extract_tp_cols='BEGIN { FS=OFS="\t" }
        {
            out = $1
            for (i = 2; i <= 12; i++) out = out OFS $i
            for (i = 13; i <= NF; i++) {
                if ($i ~ /^tp:Z:/) { out = out OFS $i; break }
            }
            print out
        }'

        # --- ENCODE ---
        encode_runtime_sum=0
        encode_memory_sum=0
        for _rep in $(seq 1 $num_repeats); do
            \time -v "$cigzip" encode \
                --paf "$input_paf" \
                --type $tp_type \
                --max-complexity $mc \
                $cmd_args \
                -t "$THREADS" \
                > "$tp_paf_full" 2> "$encode_log"
            _rep_runtime=$(parse_time_log "$encode_log")
            _rep_memory=$(parse_memory_log "$encode_log")
            encode_runtime_sum=$(echo "$encode_runtime_sum + $_rep_runtime" | bc -l)
            encode_memory_sum=$((encode_memory_sum + _rep_memory))
        done
        encode_runtime=$(echo "scale=6; $encode_runtime_sum / $num_repeats" | bc -l | sed 's/^\./0./')
        encode_memory=$((encode_memory_sum / num_repeats))

        # Subset to 12 cols + tp:Z:
        awk "$extract_tp_cols" "$tp_paf_full" > "$tp_paf"
        size_tp=$(wc -c < "$tp_paf")

        num_tps=$(awk -F'\t' '{
            for (i = 13; i <= NF; i++) {
                if ($i ~ /^tp:Z:/) {
                    val = substr($i, 6)
                    if (val != "") { n = split(val, arr, ";"); sum += n }
                    break
                }
            }
        } END { print sum + 0 }' "$tp_paf")

        # --- COMPRESS ---
        compress_runtime_sum=0
        compress_memory_sum=0
        for _rep in $(seq 1 $num_repeats); do
            \time -v "$cigzip" compress \
                --input "$tp_paf" \
                --output "$tpa_file" \
                --type $tp_type \
                --max-complexity $mc \
                $cmd_args \
                $strategy_args \
                --distance edit \
                -t "$THREADS" \
                2> "$compress_log"
            _rep_runtime=$(parse_time_log "$compress_log")
            _rep_memory=$(parse_memory_log "$compress_log")
            compress_runtime_sum=$(echo "$compress_runtime_sum + $_rep_runtime" | bc -l)
            compress_memory_sum=$((compress_memory_sum + _rep_memory))
        done
        compress_runtime=$(echo "scale=6; $compress_runtime_sum / $num_repeats" | bc -l | sed 's/^\./0./')
        compress_memory=$((compress_memory_sum / num_repeats))
        size_tpa=$(wc -c < "$tpa_file")

        # --- DECOMPRESS ---
        decompress_runtime_sum=0
        decompress_memory_sum=0
        for _rep in $(seq 1 $num_repeats); do
            \time -v "$cigzip" decompress \
                --input "$tpa_file" \
                --output "$decomp_paf" \
                2> "$decompress_log"
            _rep_runtime=$(parse_time_log "$decompress_log")
            _rep_memory=$(parse_memory_log "$decompress_log")
            decompress_runtime_sum=$(echo "$decompress_runtime_sum + $_rep_runtime" | bc -l)
            decompress_memory_sum=$((decompress_memory_sum + _rep_memory))
        done
        decompress_runtime=$(echo "scale=6; $decompress_runtime_sum / $num_repeats" | bc -l | sed 's/^\./0./')
        decompress_memory=$((decompress_memory_sum / num_repeats))

        # Correctness check
        diff_count=$(diff "$tp_paf" "$decomp_paf" | wc -l)
        [ "$diff_count" -eq 0 ] && decompress_correct="true" || decompress_correct="false"

        # --- DECODE (--no-banded) ---
        decode_runtime_sum=0
        decode_memory_sum=0
        for _rep in $(seq 1 $num_repeats); do
            \time -v "$cigzip" decode \
                --paf "$decomp_paf" \
                --sequence-files "$seq_file" \
                --type $tp_type \
                --max-complexity $mc \
                $cmd_args \
                --distance edit \
                --memory-mode $memory_mode \
                --no-banded \
                -t "$THREADS" \
                > "$decode_paf" 2> "$decode_log"
            _rep_runtime=$(parse_time_log "$decode_log")
            _rep_memory=$(parse_memory_log "$decode_log")
            decode_runtime_sum=$(echo "$decode_runtime_sum + $_rep_runtime" | bc -l)
            decode_memory_sum=$((decode_memory_sum + _rep_memory))
        done
        decode_runtime=$(echo "scale=6; $decode_runtime_sum / $num_repeats" | bc -l | sed 's/^\./0./')
        decode_memory=$((decode_memory_sum / num_repeats))

        # --- DECODE (banded, default) ---
        decode_heur_runtime_sum=0
        decode_heur_memory_sum=0
        for _rep in $(seq 1 $num_repeats); do
            \time -v "$cigzip" decode \
                --paf "$decomp_paf" \
                --sequence-files "$seq_file" \
                --type $tp_type \
                --max-complexity $mc \
                $cmd_args \
                --distance edit \
                --memory-mode $memory_mode \
                -t "$THREADS" \
                > "$decode_heur_paf" 2> "$decode_heur_log"
            _rep_runtime=$(parse_time_log "$decode_heur_log")
            _rep_memory=$(parse_memory_log "$decode_heur_log")
            decode_heur_runtime_sum=$(echo "$decode_heur_runtime_sum + $_rep_runtime" | bc -l)
            decode_heur_memory_sum=$((decode_heur_memory_sum + _rep_memory))
        done
        decode_heur_runtime=$(echo "scale=6; $decode_heur_runtime_sum / $num_repeats" | bc -l | sed 's/^\./0./')
        decode_heur_memory=$((decode_heur_memory_sum / num_repeats))

        # Write result
        echo -e "$l\t$e\t$tp_type\t$cm\t$mc\t$memory_mode\t$size_cigar\t$size_cigar_bgzip\t$bgzip_decomp_runtime\t$bgzip_decomp_memory\t$align_runtime\t$align_memory\t$size_tp\t$num_tps\t$encode_runtime\t$encode_memory\t$size_tpa\t$compress_runtime\t$compress_memory\t$decompress_runtime\t$decompress_memory\t$decompress_correct\t$decode_runtime\t$decode_memory\t$decode_heur_runtime\t$decode_heur_memory" > "$TMP_DIR/$full_prefix.result.tsv"

        # Move outputs to final locations
        [ -f "$tp_paf" ] && mv "$tp_paf" "$dir_base/simulated-data/encode/"
        [ -f "$encode_log" ] && mv "$encode_log" "$dir_base/simulated-data/encode/$full_prefix.log"
        [ -f "$tpa_file" ] && mv "$tpa_file" "$dir_base/simulated-data/compress/"
        [ -f "$compress_log" ] && mv "$compress_log" "$dir_base/simulated-data/compress/$full_prefix.log"
        [ -f "$decomp_paf" ] && mv "$decomp_paf" "$dir_base/simulated-data/decompress/$full_prefix.tp.paf"
        [ -f "$decompress_log" ] && mv "$decompress_log" "$dir_base/simulated-data/decompress/$full_prefix.log"
        [ -f "$decode_paf" ] && mv "$decode_paf" "$dir_base/simulated-data/decode/"
        [ -f "$decode_heur_paf" ] && mv "$decode_heur_paf" "$dir_base/simulated-data/decode/"
        [ -f "$decode_log" ] && mv "$decode_log" "$dir_base/simulated-data/decode/$full_prefix.log"
        [ -f "$decode_heur_log" ] && mv "$decode_heur_log" "$dir_base/simulated-data/decode/$full_prefix.heuristic.log"

        rm -f "$tp_paf_full"
        rmdir "$job_scratch" 2>/dev/null || true

    done < "$TMP_DIR/jobs.txt"

    # -------------------------------------------------------------------------
    # Merge results
    # -------------------------------------------------------------------------
    log "Merging results..."

    # Header
    echo -e "l\te\ttp_type\tcm\tmc\tmemory_mode\tsize_cigar_bytes\tsize_cigar_bgzip_bytes\tbgzip_decompress_runtime_sec\tbgzip_decompress_memory_kb\talign_runtime_sec\talign_memory_kb\tsize_tp_bytes\tnum_tracepoints\tencode_runtime_sec\tencode_memory_kb\tsize_tpa_bytes\tcompress_runtime_sec\tcompress_memory_kb\tdecompress_runtime_sec\tdecompress_memory_kb\tdecompress_correct\tdecode_runtime_sec\tdecode_memory_kb\tdecode_heuristic_runtime_sec\tdecode_heuristic_memory_kb" > "$REPORT"

    # Merge all results (sorted)
    cat "$TMP_DIR"/*.result.tsv | sort -t$'\t' -k1,1n -k2,2n -k3,3 -k4,4 -k5,5n -k6,6 >> "$REPORT"

    # Cleanup
    rm -rf "$TMP_DIR"
    rm -rf "$scratch_dir/jobs"

    log "Benchmark complete! Results in: $REPORT"
}

# =============================================================================
# STEP 4: GENERATE PLOTS
# =============================================================================

generate_plots() {
    log "Generating plots..."

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    # Main figures
    if [ -f "$SCRIPT_DIR/plotting/plot-main-figures.R" ]; then
        log "Running plot-main-figures.R..."
        Rscript "$SCRIPT_DIR/plotting/plot-main-figures.R"
    else
        log "WARNING: plot-main-figures.R not found"
    fi

    # Supplementary figures
    if [ -f "$SCRIPT_DIR/plotting/plot-supp-figures.R" ]; then
        log "Running plot-supp-figures.R..."
        Rscript "$SCRIPT_DIR/plotting/plot-supp-figures.R"
    else
        log "WARNING: plot-supp-figures.R not found"
    fi

    log "Plotting complete!"
}

# =============================================================================
# MAIN
# =============================================================================

main() {
    local run_generate=false
    local run_benchmark=false
    local run_plot=false

    # Parse arguments
    if [ $# -eq 0 ]; then
        # No arguments: run all steps
        run_generate=true
        run_benchmark=true
        run_plot=true
    else
        for arg in "$@"; do
            case $arg in
                --generate)
                    run_generate=true
                    ;;
                --benchmark)
                    run_benchmark=true
                    ;;
                --plot)
                    run_plot=true
                    ;;
                --help|-h)
                    echo "Usage: $0 [--generate] [--benchmark] [--plot]"
                    echo ""
                    echo "Options:"
                    echo "  --generate   Run sequence generation and alignment"
                    echo "  --benchmark  Run the full benchmark pipeline"
                    echo "  --plot       Generate plots from results"
                    echo "  (no args)    Run all steps"
                    echo ""
                    echo "Environment variables for configuration:"
                    echo "  DIR_BASE          Base output directory"
                    echo "  SCRATCH_DIR       Fast scratch storage"
                    echo "  GENERATE_DATASET  Path to WFA2 generate_dataset binary"
                    echo "  ALIGN_BENCHMARK   Path to WFA2 align_benchmark binary"
                    echo "  OUT_TO_PAF        Path to out-to-paf.py script"
                    echo "  PAFCHECK          Path to pafcheck binary"
                    echo "  CIGZIP            Path to cigzip binary"
                    echo "  NUM_RECORDS       Number of sequence pairs (default: 10000)"
                    echo "  THREADS           Number of threads (default: nproc)"
                    exit 0
                    ;;
                *)
                    echo "Unknown option: $arg"
                    echo "Use --help for usage information"
                    exit 1
                    ;;
            esac
        done
    fi

    log "Starting simulated data benchmarking pipeline"
    log "Configuration:"
    log "  DIR_BASE=$dir_base"
    log "  SCRATCH_DIR=$scratch_dir"
    log "  NUM_RECORDS=$NUM_RECORDS"
    log "  THREADS=$THREADS"

    if $run_generate; then
        generate_data
    fi

    if $run_benchmark; then
        run_benchmark
    fi

    if $run_plot; then
        generate_plots
    fi

    log "Pipeline complete!"
}

main "$@"
