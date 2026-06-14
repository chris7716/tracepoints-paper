# Tracepoints

Paths:

```shell
dir_base=/moosefs/guarracino/tracepoints
scratch_dir=/scratch

generate_dataset=/moosefs/guarracino/git/WFA2-lib/bin/generate_dataset
align_benchmark=/moosefs/guarracino/git/WFA2-lib/bin/align_benchmark
# To built these two on our HPC
# unset LIBRARY_PATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH;
# CC=/usr/bin/gcc BUILD_EXAMPLES=0
# LD_FLAGS='-Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 -Wl,-rpath,/lib/x86_64-linux-gnu'
# make

out_to_paf=$dir_base/scripts/out-to-paf.py
pafcheck=/moosefs/guarracino/git/pafcheck/target/release/pafcheck

cigzip=$scratch_dir/cigzip/target/release/cigzip
```

## Simulated data

Generate:

```shell
mkdir -p $dir_base/simulated-data/alns
mkdir -p $dir_base/simulated-data/seqs

num_records=10000

for l in 100 1000 10000 100000; do
  for e in 0.01 0.05 0.10 0.20; do
    echo $l $e
    l_nodot=$(echo $l | sed 's/\.//g')
    e_nodot=$(echo $e | sed 's/\.//g')
    prefix=set_${l_nodot}_${e_nodot}
    
    # Generate sequences
    $generate_dataset -n $num_records -l $l -e $e -o $dir_base/simulated-data/seqs/$prefix.seq

    # Generate ordered output with single thread (parallel mode doesn't preserve order)
    $align_benchmark -i $dir_base/simulated-data/seqs/$prefix.seq -a edit-wfa --wfa-memory high --output $dir_base/simulated-data/alns/$prefix.out -t 1

    # Convert SEQ to FASTA
    awk '
      BEGIN { counter=1 }
      /^>/ {
          print ">target_" counter
          print substr($0, 2)
      }
      /^</ {
          print ">query_" counter
          print substr($0, 2)
          counter++
      }' $dir_base/simulated-data/seqs/$prefix.seq > $dir_base/simulated-data/seqs/$prefix.fa

    # Convert OUT to PAF
    python3 $out_to_paf $dir_base/simulated-data/alns/$prefix.out $dir_base/simulated-data/alns/$prefix.paf
    $pafcheck --paf $dir_base/simulated-data/alns/$prefix.paf --sequence-files $dir_base/simulated-data/seqs/$prefix.fa --threads 48
  done
done
```

Simulation - PAF CIGAR -> PAF TRACEPOINTS -> TPA -> PAF TRACEPOINTS -> PAF CIGAR:

```shell
# Create all output directories
mkdir -p $dir_base/simulated-data/{encode,compress,decompress,decode}
mkdir -p $scratch_dir/benchmark_tmp

# Unified report (final location)
REPORT="$dir_base/simulated-data/benchmark.results.tsv"
TMP_DIR="$scratch_dir/benchmark_tmp"

# Precompute CIGAR sizes (fast, run sequentially)
for l in 100 1000 10000 100000; do
  for e in 0.01 0.05 0.10 0.20; do
    l_nodot=$(echo $l | sed 's/\.//g')
    e_nodot=$(echo $e | sed 's/\.//g')
    prefix=set_${l_nodot}_${e_nodot}
    input_paf="$dir_base/simulated-data/alns/$prefix.paf"
    
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

    # Compress the CIGAR-only PAF with bgzip -l9 and save to compress folder
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
    \time -v bgzip -d -@ $(nproc) -c "$bgzip_file" > $scratch_dir/${prefix}.cigar.paf 2> "$bgzip_decomp_log"
    rm -f $scratch_dir/${prefix}.cigar.paf

    bgzip_decomp_runtime_raw=$(grep "Elapsed (wall clock)" "$bgzip_decomp_log" | sed 's/.*: //')
    bgzip_decomp_runtime=$(echo "$bgzip_decomp_runtime_raw" | awk -F: '{
        if (NF == 3) print $1 * 3600 + $2 * 60 + $3
        else if (NF == 2) print $1 * 60 + $2
        else print $1
    }')
    bgzip_decomp_memory=$(grep "Maximum resident set size" "$bgzip_decomp_log" | sed 's/.*: //')

    # Compute baseline alignment time using cigzip with single tracepoint per alignment
    # This ensures fair comparison (same tool for baseline and tracepoint decode)
    baseline_tp_paf="$TMP_DIR/$prefix.baseline.tp.paf"
    baseline_decode_paf="$TMP_DIR/$prefix.baseline.decoded.paf"
    baseline_log="$TMP_DIR/$prefix.baseline.log"

    # Encode with very high max-complexity = 1 tracepoint per alignment
    $cigzip encode \
        --paf "$input_paf" \
        --type standard \
        --complexity-metric edit-distance \
        --max-complexity 9999999 \
        > "$baseline_tp_paf"

    # Decode (compute full optimal alignment) and measure runtime/memory
    # Repeat 10 times for l<=10000 to get stable averages; once otherwise
    if [ "$l" -le 10000 ]; then
        baseline_repeats=1
    else
        baseline_repeats=1
    fi

    align_runtime_sum=0
    align_memory_sum=0
    for _rep in $(seq 1 $baseline_repeats); do
        \time -v $cigzip decode \
            --paf "$baseline_tp_paf" \
            --sequence-files "$dir_base/simulated-data/seqs/$prefix.fa" \
            --type standard \
            --complexity-metric edit-distance \
            --max-complexity 9999999 \
            --distance edit \
            --memory-mode high \
            -t $(nproc) \
            > "$baseline_decode_paf" 2> "$baseline_log"

        _rep_runtime_raw=$(grep "Elapsed (wall clock)" "$baseline_log" | sed 's/.*: //')
        _rep_runtime=$(echo "$_rep_runtime_raw" | awk -F: '{
            if (NF == 3) print $1 * 3600 + $2 * 60 + $3
            else if (NF == 2) print $1 * 60 + $2
            else print $1
        }')
        _rep_memory=$(grep "Maximum resident set size" "$baseline_log" | sed 's/.*: //')
        align_runtime_sum=$(echo "$align_runtime_sum + $_rep_runtime" | bc -l)
        align_memory_sum=$((align_memory_sum + _rep_memory))
    done
    align_runtime=$(echo "scale=6; $align_runtime_sum / $baseline_repeats" | bc -l | sed 's/^\./0./')
    align_memory=$((align_memory_sum / baseline_repeats))

    # Move baseline log to permanent location, cleanup temporary files
    mv "$baseline_log" "$dir_base/simulated-data/alns/$prefix.baseline.log"
    rm -f "$baseline_tp_paf" "$baseline_decode_paf"

    echo "${l}_${e} $size_cigar $size_cigar_bgzip $bgzip_decomp_runtime $bgzip_decomp_memory $align_runtime $align_memory" >> "$TMP_DIR/cigar_sizes.txt"
  done
done

# AWK script to extract 12 mandatory columns + tp:Z: tag
EXTRACT_TP_COLS='BEGIN { FS=OFS="\t" }
{
    out = $1
    for (i = 2; i <= 12; i++) out = out OFS $i
    for (i = 13; i <= NF; i++) {
        if ($i ~ /^tp:Z:/) { out = out OFS $i; break }
    }
    print out
}'

# Export the benchmark function
run_benchmark() {
    local l="$1"
    local e="$2"
    local tp_type="$3"
    local cm="$4"
    local mc="$5"
    local memory_mode="$6"
    local dir_base="$7"
    local cigzip="$8"
    local size_cigar="$9"
    local size_cigar_bgzip="${10}"
    local bgzip_decomp_runtime="${11}"
    local bgzip_decomp_memory="${12}"
    local align_runtime="${13}"
    local align_memory="${14}"
    local TMP_DIR="${15}"
    
    l_nodot=$(echo $l | sed 's/\.//g')
    e_nodot=$(echo $e | sed 's/\.//g')
    prefix=set_${l_nodot}_${e_nodot}
    input_paf="$dir_base/simulated-data/alns/$prefix.paf"
    seq_file="$dir_base/simulated-data/seqs/$prefix.fa"
    
    prefix2=$tp_type.$cm.$mc.$memory_mode
    full_prefix="$prefix.$prefix2"
    [ "$tp_type" = "fastga" ] && cmd_args="" || cmd_args="--complexity-metric $cm"
    # Use no compression layer for all methods. The empty strategy (standard) defaults to automatic, which selects nocomp.
    # For fastga we name the rice strategy explicitly, so we must add -nocomp; otherwise an explicit strategy defaults to a zstd wrapper.
    [ "$tp_type" = "fastga" ] && strategy_args="--strategy rice-nocomp;rice-nocomp" || strategy_args=""

    # Repeat 10 times for mc<=64 or fastga to get stable averages; once otherwise
    if [ "$mc" -le 64 ] || [ "$tp_type" = "fastga" ]; then
        num_repeats=10
    else
        num_repeats=1
    fi

    # Use unique scratch directory per job
    job_scratch="$scratch_dir/jobs/$full_prefix"
    mkdir -p "$job_scratch"
    
    # File paths (all in scratch)
    tp_paf_full="$job_scratch/$full_prefix.tp.full.paf"      # Full output from encode
    tp_paf="$job_scratch/$full_prefix.tp.paf"                # Subsetted: 12 cols + tp:Z:
    tpa_file="$job_scratch/$full_prefix.tpa"
    decomp_paf="$job_scratch/$full_prefix.decomp.tp.paf"
    decode_paf="$job_scratch/$full_prefix.paf"
    decode_heur_paf="$job_scratch/$full_prefix.heuristic.paf"
    encode_log="$job_scratch/$full_prefix.encode.log"
    compress_log="$job_scratch/$full_prefix.compress.log"
    decompress_log="$job_scratch/$full_prefix.decompress.log"
    decode_log="$job_scratch/$full_prefix.decode.log"
    decode_heur_log="$job_scratch/$full_prefix.decode.heuristic.log"
    
    # AWK script to extract 12 mandatory columns + tp:Z: tag
    extract_tp_cols='BEGIN { FS=OFS="\t" }
    {
        out = $1
        for (i = 2; i <= 12; i++) out = out OFS $i
        for (i = 13; i <= NF; i++) {
            if ($i ~ /^tp:Z:/) { out = out OFS $i; break }
        }
        print out
    }'
    
    # Helper functions
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
    
    # --- ENCODE ---
    encode_runtime_sum=0
    encode_memory_sum=0
    for _rep in $(seq 1 $num_repeats); do
        \time -v $cigzip encode \
            --paf "$input_paf" \
            --type $tp_type \
            --max-complexity $mc \
            $cmd_args \
            -t $(nproc) \
            > "$tp_paf_full" 2> "$encode_log"
        _rep_runtime=$(parse_time_log "$encode_log")
        _rep_memory=$(parse_memory_log "$encode_log")
        encode_runtime_sum=$(echo "$encode_runtime_sum + $_rep_runtime" | bc -l)
        encode_memory_sum=$((encode_memory_sum + _rep_memory))
    done
    encode_runtime=$(echo "scale=6; $encode_runtime_sum / $num_repeats" | bc -l | sed 's/^\./0./')
    encode_memory=$((encode_memory_sum / num_repeats))

    # --- SUBSET TO 12 COLS + tp:Z: FOR FAIR COMPARISON ---
    awk "$extract_tp_cols" "$tp_paf_full" > "$tp_paf"
    
    # Measure subsetted TP size (this is what we compare against CIGAR size)
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
    
    # --- COMPRESS (using subsetted PAF) ---
    compress_runtime_sum=0
    compress_memory_sum=0
    for _rep in $(seq 1 $num_repeats); do
        \time -v $cigzip compress \
            --input "$tp_paf" \
            --output "$tpa_file" \
            --type $tp_type \
            --max-complexity $mc \
            $cmd_args \
            $strategy_args \
            --distance edit \
            -t $(nproc) \
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
        \time -v $cigzip decompress \
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

    # Correctness check (compare subsetted original vs decompressed)
    diff_count=$(diff "$tp_paf" "$decomp_paf" | wc -l)
    [ "$diff_count" -eq 0 ] && decompress_correct="true" || decompress_correct="false"
    
    # --- DECODE (--no-banded) ---
    decode_runtime_sum=0
    decode_memory_sum=0
    for _rep in $(seq 1 $num_repeats); do
        \time -v $cigzip decode \
            --paf "$decomp_paf" \
            --sequence-files "$seq_file" \
            --type $tp_type \
            --max-complexity $mc \
            $cmd_args \
            --distance edit \
            --memory-mode $memory_mode \
            --no-banded \
            -t $(nproc) \
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
        \time -v $cigzip decode \
            --paf "$decomp_paf" \
            --sequence-files "$seq_file" \
            --type $tp_type \
            --max-complexity $mc \
            $cmd_args \
            --distance edit \
            --memory-mode $memory_mode \
            -t $(nproc) \
            > "$decode_heur_paf" 2> "$decode_heur_log"
        _rep_runtime=$(parse_time_log "$decode_heur_log")
        _rep_memory=$(parse_memory_log "$decode_heur_log")
        decode_heur_runtime_sum=$(echo "$decode_heur_runtime_sum + $_rep_runtime" | bc -l)
        decode_heur_memory_sum=$((decode_heur_memory_sum + _rep_memory))
    done
    decode_heur_runtime=$(echo "scale=6; $decode_heur_runtime_sum / $num_repeats" | bc -l | sed 's/^\./0./')
    decode_heur_memory=$((decode_heur_memory_sum / num_repeats))
    
    # --- WRITE RESULT TO TEMP FILE ---
    echo -e "$l\t$e\t$tp_type\t$cm\t$mc\t$memory_mode\t$size_cigar\t$size_cigar_bgzip\t$bgzip_decomp_runtime\t$bgzip_decomp_memory\t$align_runtime\t$align_memory\t$size_tp\t$num_tps\t$encode_runtime\t$encode_memory\t$size_tpa\t$compress_runtime\t$compress_memory\t$decompress_runtime\t$decompress_memory\t$decompress_correct\t$decode_runtime\t$decode_memory\t$decode_heur_runtime\t$decode_heur_memory" > "$TMP_DIR/$full_prefix.result.tsv"
    
    # Remove full PAF (keep only subsetted version)
    rm -f "$tp_paf_full"
    
    echo "Completed: $full_prefix"
}

export -f run_benchmark

# Generate all job parameters
> "$TMP_DIR/jobs.txt"
for l in 100 1000 10000 100000; do
  for e in 0.01 0.05 0.10 0.20; do
    # Look up precomputed values
    size_cigar=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f2)
    size_cigar_bgzip=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f3)
    bgzip_decomp_runtime=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f4)
    bgzip_decomp_memory=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f5)
    align_runtime=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f6)
    align_memory=$(grep "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt" | cut -d' ' -f7)

    for tp_type in fastga standard; do
      if [ "$tp_type" = "fastga" ]; then
          cm_list="none"
          mc_list="100"
      else
          cm_list="edit-distance diagonal-distance"
          mc_list="32 64 128 256 512 1024"
      fi

      for cm in $cm_list; do
        for mc in $mc_list; do
          for memory_mode in high; do
            echo "$l $e $tp_type $cm $mc $memory_mode $dir_base $cigzip $size_cigar $size_cigar_bgzip $bgzip_decomp_runtime $bgzip_decomp_memory $align_runtime $align_memory $TMP_DIR" >> "$TMP_DIR/jobs.txt"
          done
        done
      done
    done
  done
done

# Create scratch jobs directory
mkdir -p $scratch_dir/jobs

# Run in parallel (adjust -j for number of concurrent jobs)
parallel -j 1 --colsep ' ' run_benchmark {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} < "$TMP_DIR/jobs.txt"

# --- MOVE ALL OUTPUT FILES FROM SCRATCH TO FINAL DESTINATIONS ---
echo "Moving output files to final destinations..."

for job_dir in $scratch_dir/jobs/*/; do
    [ -d "$job_dir" ] || continue
    
    full_prefix=$(basename "$job_dir")
    
    # Encode outputs (subsetted PAF)
    [ -f "$job_dir/$full_prefix.tp.paf" ] && mv "$job_dir/$full_prefix.tp.paf" "$dir_base/simulated-data/encode/"
    [ -f "$job_dir/$full_prefix.encode.log" ] && mv "$job_dir/$full_prefix.encode.log" "$dir_base/simulated-data/encode/$full_prefix.log"
    
    # Compress outputs
    [ -f "$job_dir/$full_prefix.tpa" ] && mv "$job_dir/$full_prefix.tpa" "$dir_base/simulated-data/compress/"
    [ -f "$job_dir/$full_prefix.compress.log" ] && mv "$job_dir/$full_prefix.compress.log" "$dir_base/simulated-data/compress/$full_prefix.log"
    
    # Decompress outputs
    [ -f "$job_dir/$full_prefix.decomp.tp.paf" ] && mv "$job_dir/$full_prefix.decomp.tp.paf" "$dir_base/simulated-data/decompress/$full_prefix.tp.paf"
    [ -f "$job_dir/$full_prefix.decompress.log" ] && mv "$job_dir/$full_prefix.decompress.log" "$dir_base/simulated-data/decompress/$full_prefix.log"
    
    # Decode outputs
    [ -f "$job_dir/$full_prefix.paf" ] && mv "$job_dir/$full_prefix.paf" "$dir_base/simulated-data/decode/"
    [ -f "$job_dir/$full_prefix.heuristic.paf" ] && mv "$job_dir/$full_prefix.heuristic.paf" "$dir_base/simulated-data/decode/"
    [ -f "$job_dir/$full_prefix.decode.log" ] && mv "$job_dir/$full_prefix.decode.log" "$dir_base/simulated-data/decode/$full_prefix.log"
    [ -f "$job_dir/$full_prefix.decode.heuristic.log" ] && mv "$job_dir/$full_prefix.decode.heuristic.log" "$dir_base/simulated-data/decode/$full_prefix.heuristic.log"
    
    # Remove job directory
    rmdir "$job_dir" 2>/dev/null
done

# Header
echo -e "l\te\ttp_type\tcm\tmc\tmemory_mode\tsize_cigar_bytes\tsize_cigar_bgzip_bytes\tbgzip_decompress_runtime_sec\tbgzip_decompress_memory_kb\talign_runtime_sec\talign_memory_kb\tsize_tp_bytes\tnum_tracepoints\tencode_runtime_sec\tencode_memory_kb\tsize_tpa_bytes\tcompress_runtime_sec\tcompress_memory_kb\tdecompress_runtime_sec\tdecompress_memory_kb\tdecompress_correct\tdecode_runtime_sec\tdecode_memory_kb\tdecode_heuristic_runtime_sec\tdecode_heuristic_memory_kb" > "$REPORT"

# Merge all results into final report (sorted for reproducibility)
cat "$TMP_DIR"/*.result.tsv | sort -t$'\t' -k1,1n -k2,2n -k3,3 -k4,4 -k5,5n -k6,6 >> "$REPORT"

# Cleanup scratch
rm -rf "$TMP_DIR"
rm -rf $scratch_dir/jobs

echo "Benchmark complete. Results in: $REPORT"
```

### Indel-skewed check

The simulator draws each edit as mismatch/insertion/deletion with equal probability, so simulated indels are balanced (I/D ~ 1.0). The `--indels NUM,LENGTH` option adds extra deletions into the text, which appear as insertions in the query->target CIGAR, giving an insertion-skewed dataset to test DB-TP's sensitivity to an indel bias. We run at 10 Kb so the tracepoint streams (not per-record metadata) dominate the TPA, making the storage ratios meaningful. `--indels 2200` gives I/D ~ 5 at 10 Kb (calibrate per length; the simulator is unseeded so exact values vary):

```shell
mkdir -p $scratch_dir/indeltest && cd $scratch_dir/indeltest

# Balanced (I/D ~ 1.0) vs insertion-skewed (I/D ~ 5)
$generate_dataset -n 1000 -l 10000 -e 0.20                 -o bal.seq
$generate_dataset -n 1000 -l 10000 -e 0.20 --indels 2200,1 -o skew.seq

for f in bal skew; do
  $align_benchmark -i $f.seq -a edit-wfa --wfa-memory high --output $f.out -t 1
  python3 $out_to_paf $f.out $f.paf
done

# Per-condition I/D, tracepoint counts, TPA size, and ratio (TPA/PAF), mc=32.
# Representative run (simulator unseeded; exact values vary):
#   bal  (I/D 1.00): DB-TP tps=3216   ratio=0.0125 | EB-TP tps=53678  ratio=0.0277
#   skew (I/D 4.92): DB-TP tps=65919  ratio=0.0291 | EB-TP tps=75132  ratio=0.0265
# => under skew DB-TP's ratio rises (0.012 -> 0.029) and overtakes EB-TP (~0.027), which is unchanged.
for f in bal skew; do
  paf=$(stat -c%s $f.paf)
  id=$(awk -F'\t' '{for(i=13;i<=NF;i++) if(substr($i,1,5)=="cg:Z:"){s=substr($i,6);num="";
        for(j=1;j<=length(s);j++){c=substr(s,j,1);if(c~/[0-9]/)num=num c;
          else{v=num+0;if(c=="I")I+=v;else if(c=="D")D+=v;num=""}}}} END{printf "%.2f",I/D}' $f.paf)
  echo "=== $f (I/D=$id, PAF=$paf B) ==="
  for m in diagonal-distance edit-distance; do
    $cigzip encode --paf $f.paf --type standard --complexity-metric $m --max-complexity 32 --distance edit -o $f.$m.tp.paf
    tps=$(grep -o 'tp:Z:[^[:space:]]*' $f.$m.tp.paf | awk -F';' '{s+=NF} END{print s}')
    $cigzip compress --input $f.$m.tp.paf --type standard --complexity-metric $m --max-complexity 32 --output $f.$m.tpa
    tpa=$(stat -c%s $f.$m.tpa)
    awk -v m=$m -v tp=$tps -v t=$tpa -v p=$paf 'BEGIN{printf "  %-18s tps=%d TPA=%d B ratio=%.4f\n",m,tp,t,t/p}'
  done
done
```

## Real data

### Download

Human pangenome:

```shell
mkdir $scratch_dir/hprcv2-25k
cd $scratch_dir/hprcv2-25k

# Alignments
aws s3 sync s3://garrisonlab/hprcv2/pafs/all-vs-1/ . --no-sign-request
```

T2T ape pangenome:

```shell
mkdir -p $scratch_dir/t2t-ape-pangenome
cd $scratch_dir/t2t-ape-pangenome

# Sequences
wget -c https://garrisonlab.s3.amazonaws.com/t2t-primates/primates16.20240512.fa.gz

# Alignments
aws s3 sync s3://garrisonlab/t2t-primates/wfmash-v0.13.0/alignments_fixedSiamang/ . --no-sign-request
```

### Statistics

Sizes:

```shell
cd $scratch_dir/hprcv2-25k
find . -maxdepth 1 -name "*.paf" -print0 | xargs -0 cat | wc -c | awk '{print $1/1024/1024/1024}' # 940.11 GBs
find . -maxdepth 1 -name "*.paf.gz" -print0 | xargs -0 cat | wc -c | awk '{print $1/1024/1024/1024}' # 296.507 GBs
# CIGAR strings only (cg:Z: tag content, excluding the "cg:Z:" prefix)
find . -maxdepth 1 -name "*.paf" -print0 | xargs -0 cat | awk -F'\t' '{ for (i=13; i<=NF; i++) if (substr($i,1,5)=="cg:Z:") total += length($i) - 5 } END { print total/1024/1024/1024 " GBs" }' # 873.496 GBs

cd $scratch_dir/t2t-ape-pangenome
# excluding HPRCy1-vs-*
find . -maxdepth 1 -name "*.paf" ! -name "HPRCy1-vs-*" -print0 | xargs -0 cat | wc -c | awk '{print $1/1024/1024/1024}' # 78.9188 GBs
find . -maxdepth 1 -name "*.paf.gz" ! -name "HPRCy1-vs-*" -print0 | xargs -0 cat | wc -c | awk '{print $1/1024/1024/1024}' # 21.243 GBs
# CIGAR strings only (cg:Z: tag content, excluding the "cg:Z:" prefix)
find . -maxdepth 1 -name "*.paf" ! -name "HPRCy1-vs-*" -print0 | xargs -0 cat | awk -F'\t' '{ for (i=13; i<=NF; i++) if (substr($i,1,5)=="cg:Z:") total += length($i) - 5 } END { print total/1024/1024/1024 " GBs" }' # 78.8078 GBs
# including HPRCy1-vs-*
find . -maxdepth 1 -name "*.paf" -print0 | xargs -0 cat | wc -c | awk '{print $1/1024/1024/1024}' # 545.929 GBs
find . -maxdepth 1 -name "*.paf.gz" -print0 | xargs -0 cat | wc -c | awk '{print $1/1024/1024/1024}' # 147.705 GBs
# CIGAR strings only (cg:Z: tag content, excluding the "cg:Z:" prefix)
find . -maxdepth 1 -name "*.paf" -print0 | xargs -0 cat | awk -F'\t' '{ for (i=13; i<=NF; i++) if (substr($i,1,5)=="cg:Z:") total += length($i) - 5 } END { print total/1024/1024/1024 " GBs" }' # 545.126 GBs
```

Alignment statistics:

```shell
python3 $dir_base/scripts/paf_alignment_stats.py -t 16 --min-identity 0.0 *.paf.gz > t2t-ape-pangenome.aln-stats.tsv
python3 $dir_base/scripts/paf_alignment_stats.py -t 48 --min-identity 0.0 *.paf > hprcv2-25k.aln-stats.tsv
```

TPA random-access index size vs tracepoint data (one varint byte-offset per alignment record; the .tpa.idx file is excluded from reported storage):

```shell
# Point tpa_dir at the directory holding the *.tpa and matching *.tpa.idx files
for tpa_dir in /scratch/hprcv2-25k-tpas /scratch/t2t-ape-pangenome-tpas; do
  echo "## $tpa_dir"
  tpa=$(find "$tpa_dir" -name "*.tpa"     -print0 | xargs -0 stat -c%s | awk '{s+=$1} END {print s}')
  idx=$(find "$tpa_dir" -name "*.tpa.idx" -print0 | xargs -0 stat -c%s | awk '{s+=$1} END {print s}')
  awk -v tpa="$tpa" -v idx="$idx" 'BEGIN {
    printf "tpa = %.3f GiB, idx = %.3f GiB, idx/tpa = %.2f%%\n", tpa/1024^3, idx/1024^3, idx/tpa*100
  }'
done
# Human  (EB-TP delta=128): tpa = 23.348 GiB, idx = 1.443 GiB, idx/tpa = 6.18%
# Primate (EB-TP delta=128): tpa = 0.700 GiB, idx = 0.002 GiB, idx/tpa = 0.29%
```

Score distributions (gap-affine2p scores from CIGARs):

```shell
# Uses huge max-complexity to get 1 tracepoint per alignment (fast, just computes sc:i: scores)

# HPRCv2
cd $scratch_dir/hprcv2-25k
find . -name "*.paf.gz" | \
  parallel -j $(nproc) \
    "$cigzip encode --paf {} \
      --type standard --complexity-metric edit-distance \
      --max-complexity 9999999 \
      --distance gap-affine2p --penalties 5,8,2,24,1 \
      -t 1 2>/dev/null \
    | awk -F'\t' '{for(i=13;i<=NF;i++) if(\$i~/^sc:i:/){print substr(\$i,6)+0;next}}'" \
  | sort -n | uniq -c \
  | awk 'BEGIN{OFS="\t"; print "score","frequency"}{print $2,$1}' \
  > $dir_base/hprcv2-25k/hprcv2-25k.input-score-frequency.tsv

# T2T ape pangenome
cd $scratch_dir/t2t-ape-pangenome
find . -name "*.paf.gz" | \
  parallel -j $(nproc) \
    "$cigzip encode --paf {} \
      --type standard --complexity-metric edit-distance \
      --max-complexity 9999999 \
      --distance gap-affine2p --penalties 5,8,2,24,1 \
      -t 1 2>/dev/null \
    | awk -F'\t' '{for(i=13;i<=NF;i++) if(\$i~/^sc:i:/){print substr(\$i,6)+0;next}}'" \
  | sort -n | uniq -c \
  | awk 'BEGIN{OFS="\t"; print "score","frequency"}{print $2,$1}' \
  > $dir_base/t2t-ape-pangenome/t2t-ape-pangenome.input-score-frequency.tsv

# T2T ape pangenome - score distribution by target genome
cd $scratch_dir/t2t-ape-pangenome
find . -name "*.paf.gz" | \
  parallel -j $(nproc) \
    "$cigzip encode --paf {} \
      --type standard --complexity-metric edit-distance \
      --max-complexity 9999999 \
      --distance gap-affine2p --penalties 5,8,2,24,1 \
      -t 1 2>/dev/null \
    | awk -F'\t' '{target=\$6; sub(/#.*/, \"\", target); for(i=13;i<=NF;i++) if(\$i~/^sc:i:/){print substr(\$i,6)+0 \"\\t\" target; next}}'" \
  | sort -t$'\t' -k2,2 -k1,1n | uniq -c \
  | awk 'BEGIN{OFS="\t"; print "score","frequency","target"}{print $2,$1,$3}' \
  > $dir_base/t2t-ape-pangenome/t2t-ape-pangenome.input-score-frequency.by-target.tsv
```

Tracepoint statistics:

```shell
cd .../encode

python3 $dir_base/scripts/paf_tracepoint_stats.py -t 48 *.diagonal-distance.32.tp.paf > t2t-ape-pangenome.diagonal-distance-32.tp-stats.tsv
python3 $dir_base/scripts/paf_tracepoint_stats.py -t 48 *.edit-distance.32.tp.paf > t2t-ape-pangenome.edit-distance-32.tp-stats.tsv

python3 $dir_base/scripts/paf_tracepoint_stats.py -t 48 *.diagonal-distance.32.tp.paf > hprcv2-25k.diagonal-distance-32.tp-stats.tsv
python3 $dir_base/scripts/paf_tracepoint_stats.py -t 48 *.edit-distance.32.tp.paf > hprcv2-25k.edit-distance-32.tp-stats.tsv
```

### Benchmarking

Real data - PAF CIGAR -> PAF TRACEPOINTS -> TPA -> PAF TRACEPOINTS -> PAF CIGAR:

```shell
# Key differences from simulated data:
# - Score-based correctness check instead of CIGAR equality
#   (reconstructed score must be >= original score)
# - File sizes measured using only 12 mandatory PAF columns + relevant tag
#   (cg:Z: for CIGAR, tp:Z: for tracepoints) - excludes so:i:, sc:i:, etc.

# Benchmark parameters
MC_LIST="16 32 64 128 256"      # EB-TP/DB-TP max-complexity values (16 = bootstrap warm-up, excluded from reporting)
FASTGA_MC_LIST="500 1000 2000"  # FL-TP trace spacings (fastga scaling)
TMP_DIR="$scratch_dir/real_benchmark_tmp"
mkdir -p $TMP_DIR
HEADER="dataset\tpaf_file\tcm\tmc\tdecode_mode\tsize_cigar_bytes\tsize_tp_bytes\tnum_tracepoints\tsize_tpa_bytes\tencode_runtime_sec\tencode_memory_kb\tcompress_runtime_sec\tcompress_memory_kb\tdecompress_runtime_sec\tdecompress_memory_kb\tdecode_runtime_sec\tdecode_memory_kb\tnum_alignments\tscore_identical\tscore_improved\tscore_degraded"

# Export benchmark function
run_benchmark_real() {
    local dataset="$1"
    local input_paf="$2"
    local seq_files="$3"
    local cm="$4"
    local dir_base="$5"
    local cigzip="$6"
    local TMP_DIR="$7"
    local MC="$8"

    paf_name=$(basename "$input_paf" .paf.gz)
    paf_name=$(basename "$paf_name" .paf)  # Handle both .paf.gz and .paf

    prefix="$paf_name.$cm.$MC"
    job_scratch="$scratch_dir/real_jobs/$prefix"
    mkdir -p "$job_scratch"

    # File paths
    tp_paf_full="$job_scratch/$prefix.tp.full.paf"
    tp_paf="$job_scratch/$prefix.tp.paf"
    tpa_file="$job_scratch/$prefix.tpa"
    decomp_paf="$job_scratch/$prefix.decomp.tp.paf"
    decode_paf="$job_scratch/$prefix.decoded.paf"
    degraded_paf="$job_scratch/$prefix.degraded.paf"

    encode_log="$job_scratch/encode.log"
    compress_log="$job_scratch/compress.log"
    decompress_log="$job_scratch/decompress.log"
    decode_log="$job_scratch/decode.log"

    # AWK scripts (must be defined inside function for export)
    extract_cg_cols='BEGIN { FS=OFS="\t" }
    {
        out = $1
        for (i = 2; i <= 12; i++) out = out OFS $i
        for (i = 13; i <= NF; i++) {
            if ($i ~ /^cg:Z:/) { out = out OFS $i; break }
        }
        print out
    }'

    extract_tp_cols='BEGIN { FS=OFS="\t" }
    {
        out = $1
        for (i = 2; i <= 12; i++) out = out OFS $i
        for (i = 13; i <= NF; i++) {
            if ($i ~ /^tp:Z:/) { out = out OFS $i; break }
        }
        print out
    }'

    # AWK to extract sc:i: score from a PAF line
    extract_score='{ for (i = 13; i <= NF; i++) if ($i ~ /^sc:i:/) { print substr($i, 6) + 0; next } print "" }'

    # Helper functions
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

    # Compute CIGAR size (12 cols + cg:Z: only - excludes other tags for fair comparison)
    size_cigar=$(awk "$extract_cg_cols" "$input_paf" | wc -c)

    # Skip empty PAF files
    if [[ "$size_cigar" -lt 10 ]]; then
        echo "Skipping empty PAF: $dataset/$paf_name ($cm, mc=$MC)"
        rm -rf "$job_scratch"
        return 0
    fi

    if [ "$cm" = "fastga" ]; then
        tp_type="fastga";  cm_args="";                         strategy_args="--strategy rice-nocomp;rice-nocomp"
    elif [ "$cm" = "edit-distance" ]; then
        tp_type="standard"; cm_args="--complexity-metric $cm"; strategy_args="--strategy rice-nocomp;2d-delta-nocomp"
    else  # diagonal-distance
        tp_type="standard"; cm_args="--complexity-metric $cm"; strategy_args="--strategy raw-nocomp;2d-delta-nocomp"
    fi
    dist_args="--distance gap-affine2p --penalties 5,8,2,24,1"

    # --- ENCODE ---
    \time -v $cigzip encode \
        --paf "$input_paf" \
        --type $tp_type \
        $cm_args \
        --max-complexity $MC \
        $dist_args \
        -t $(nproc) \
        > "$tp_paf_full" 2> "$encode_log"

    encode_runtime=$(parse_time_log "$encode_log")
    encode_memory=$(parse_memory_log "$encode_log")

    # Subset to 12 cols + tp:Z: (excludes other tags for fair size comparison)
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
    \time -v $cigzip compress \
        --input "$tp_paf" \
        --output "$tpa_file" \
        --type $tp_type \
        $cm_args \
        --max-complexity $MC \
        $strategy_args \
        $dist_args \
        -t $(nproc) \
        2> "$compress_log"

    compress_runtime=$(parse_time_log "$compress_log")
    compress_memory=$(parse_memory_log "$compress_log")
    size_tpa=$(wc -c < "$tpa_file")
    rm -f "$tp_paf"

    # --- DECOMPRESS ---
    \time -v $cigzip decompress \
        --input "$tpa_file" \
        --output "$decomp_paf" \
        2> "$decompress_log"

    decompress_runtime=$(parse_time_log "$decompress_log")
    decompress_memory=$(parse_memory_log "$decompress_log")
    rm -f "$tpa_file"

    # --- DECODE (banded, default) ---
    \time -v $cigzip decode \
        --paf "$decomp_paf" \
        --sequence-list $seq_files \
        --type $tp_type \
        $cm_args \
        --max-complexity $MC \
        $dist_args \
        --memory-mode high \
        -t $(nproc) \
        > "$decode_paf" 2> "$decode_log"

    decode_runtime=$(parse_time_log "$decode_log")
    decode_memory=$(parse_memory_log "$decode_log")
    rm -f "$decomp_paf"

    # Compare scores: original (from tp_paf_full) vs reconstructed (from decode_paf)
    # Sort both files by first 9 columns to ensure matching alignment order
    # (multi-threaded output may not preserve input order)
    # Lower score = better alignment (penalties are negative)
    # Also save original input alignments (with CIGAR) that lead to degraded scores
    sorted_encoded="$job_scratch/sorted_encoded.paf"
    sorted_decode="$job_scratch/sorted_decode.paf"

    # Sort encoded and decoded PAFs (input already pre-sorted by caller)
    sort -T /scratch --parallel=$(nproc) -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 -k7,7n -k8,8n -k9,9n "$tp_paf_full" > "$sorted_encoded" &
    sort -T /scratch --parallel=$(nproc) -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 -k7,7n -k8,8n -k9,9n "$decode_paf" > "$sorted_decode" &
    wait
    rm -f "$tp_paf_full" "$decode_paf"

    read num_aln score_identical score_improved score_degraded <<< $(
        paste "$input_paf" "$sorted_encoded" "$sorted_decode" | \
        awk -F'\t' -v degraded_file="$degraded_paf" '
        BEGIN { identical=0; improved=0; degraded=0; n=0 }
        {
            n++
            # Count fields in each file (input, encoded, decoded)
            # We pasted 3 files, need to find boundaries
            # Extract sc:i: score from encoded (second file)
            orig_score = ""
            recon_score = ""

            # Find where each file ends by counting fields
            # All 3 files have same alignment, so same first 12 mandatory cols
            # Scan for sc:i: tags - first occurrence is from encoded, second from decoded
            found_orig = 0
            for (i = 13; i <= NF; i++) {
                if ($i ~ /^sc:i:/) {
                    if (!found_orig) {
                        orig_score = substr($i, 6) + 0
                        found_orig = 1
                    } else {
                        recon_score = substr($i, 6) + 0
                        break
                    }
                }
            }

            if (orig_score == "" || recon_score == "") next
            if (recon_score == orig_score) identical++
            else if (recon_score > orig_score) improved++  # scores are negative, closer to 0 = better
            else {
                degraded++
                # Output: input_paf_line (12 cols + cg:Z:), orig_score, recon_score
                # Find cg:Z: tag and build output from first file
                out = $1
                for (i = 2; i <= 12; i++) out = out "\t" $i
                for (i = 13; i <= NF; i++) {
                    if ($i ~ /^cg:Z:/) { out = out "\t" $i; break }
                }
                print out "\t" orig_score "\t" recon_score >> degraded_file
            }
        }
        END { print n, identical, improved, degraded }'
    )
    rm -f "$sorted_encoded" "$sorted_decode"

    # Write result
    echo -e "$dataset\t$paf_name\t$cm\t$MC\theuristic\t$size_cigar\t$size_tp\t$num_tps\t$size_tpa\t$encode_runtime\t$encode_memory\t$compress_runtime\t$compress_memory\t$decompress_runtime\t$decompress_memory\t$decode_runtime\t$decode_memory\t$num_aln\t$score_identical\t$score_improved\t$score_degraded" > "$TMP_DIR/$prefix.result.tsv"

    # Determine output directory and move files immediately
    if [[ "$dataset" == "hprcv2-25k" ]]; then
        out_dir="$dir_base/hprcv2-25k"
    else
        out_dir="$dir_base/t2t-ape-pangenome"
    fi

    # Move encode outputs
    [ -f "$tp_paf" ] && bgzip -l 9 -@ 96 "$tp_paf" && mv "$tp_paf.gz" "$out_dir/encode/"
    [ -f "$encode_log" ] && mv "$encode_log" "$out_dir/encode/$prefix.log"

    # Move compress outputs
    [ -f "$tpa_file" ] && mv "$tpa_file" "$out_dir/compress/"
    [ -f "$compress_log" ] && mv "$compress_log" "$out_dir/compress/$prefix.log"

    # Move decompress outputs
    [ -f "$decomp_paf" ] && bgzip -l 9 -@ 96 "$decomp_paf" && mv "$decomp_paf.gz" "$out_dir/decompress/$prefix.tp.paf.gz"
    [ -f "$decompress_log" ] && mv "$decompress_log" "$out_dir/decompress/$prefix.log"

    # Move decode outputs
    [ -f "$decode_paf" ] && bgzip -l 9 -@ 96 "$decode_paf" && mv "$decode_paf.gz" "$out_dir/decode/$prefix.paf.gz"
    [ -f "$decode_log" ] && mv "$decode_log" "$out_dir/decode/$prefix.log"
    [ -f "$degraded_paf" ] && bgzip -l 9 -@ 96 "$degraded_paf" && mv "$degraded_paf.gz" "$out_dir/decode/$prefix.degraded.paf.gz"

    # Cleanup job scratch directory
    rm -rf "$job_scratch"

    echo "Completed: $dataset/$paf_name ($cm, mc=$MC)"
}

# Create jobs directory
mkdir -p $scratch_dir/real_jobs

# =============================================================================
# HPRCv2 (Human-vs-Human) Benchmark
# =============================================================================

HPRCV2_PAFS="$dir_base/hprcv2-25k-pafs"
HPRCV2_FASTA_DIR="$scratch_dir/hprcv2-25k-fasta"
mkdir -p $dir_base/hprcv2-25k/{encode,compress,decompress,decode}
HPRCV2_SEQS="$TMP_DIR/hprcv2_seqlist.txt"
ls "$HPRCV2_FASTA_DIR"/*.fa > "$HPRCV2_SEQS"
REPORT_HPRCV2="$dir_base/hprcv2-25k/benchmark.results.tsv"
echo -e "$HEADER" > "$REPORT_HPRCV2"
# Process HPRCv2: copy one PAF.gz at a time to /scratch, decompress, run all cm/mc combos, then clean up
echo "Running HPRCv2 benchmark..."
for paf in $HPRCV2_PAFS/*.paf.gz; do
    [ -f "$paf" ] || continue
    paf_basename=$(basename "$paf")
    scratch_paf="$scratch_dir/$paf_basename"
    input_paf="$scratch_dir/${paf_basename%.paf.gz}.paf"
    cp "$paf" "$scratch_paf"
    zcat "$scratch_paf" > "$input_paf"
    rm -f "$scratch_paf"
    # EB-TP runs the full sweep; DB-TP is capped at mc<=64 (larger b gives very large segments and costly WFA
    # reconstruction on real data). mc=16 is the bootstrap warm-up.
    for cm in edit-distance diagonal-distance; do
        for mc in $MC_LIST; do
            [ "$cm" = "diagonal-distance" ] && [ "$mc" -gt 64 ] && continue
            run_benchmark_real "hprcv2-25k" "$input_paf" "$HPRCV2_SEQS" "$cm" "$dir_base" "$cigzip" "$TMP_DIR" "$mc"
        done
    done
    # FL-TP (fastga): trace spacings $FASTGA_MC_LIST (fastga scaling, like EB/DB-TP). Reconstructs under the same
    # dual gap-affine model as EB/DB-TP (see run_benchmark_real), so it is a full Table 2 entry.
    #for fl_mc in $FASTGA_MC_LIST; do
    #    run_benchmark_real "hprcv2-25k" "$input_paf" "$HPRCV2_SEQS" "fastga" "$dir_base" "$cigzip" "$TMP_DIR" "$fl_mc"
    #done
    rm -f "$input_paf"
done
# Merge HPRCv2 results
cat "$TMP_DIR"/*.result.tsv 2>/dev/null | grep "^hprcv2-25k" | sort >> "$REPORT_HPRCV2"
rm -f "$TMP_DIR"/*.result.tsv
# =============================================================================

# =============================================================================
# Primates Benchmark
# =============================================================================

PRIMATES_PAFS="$dir_base/t2t-ape-pangenome-pafs"
PRIMATES_FASTA_DIR="$scratch_dir/t2t-ape-pangenome-fasta"
mkdir -p $dir_base/t2t-ape-pangenome/{encode,compress,decompress,decode}
PRIMATES_SEQS="$TMP_DIR/primates_seqlist.txt"
ls "$PRIMATES_FASTA_DIR"/*.fa > "$PRIMATES_SEQS"
REPORT_PRIMATES="$dir_base/t2t-ape-pangenome/benchmark.results.tsv"
echo -e "$HEADER" > "$REPORT_PRIMATES"
# Process Primates: copy one PAF.gz at a time to /scratch, decompress, run all cm/mc combos, then clean up
echo "Running Primates benchmark..."
for paf in $PRIMATES_PAFS/*.paf.gz; do
    [ -f "$paf" ] || continue
    paf_basename=$(basename "$paf")
    scratch_paf="$scratch_dir/$paf_basename"
    input_paf="$scratch_dir/${paf_basename%.paf.gz}.paf"
    cp "$paf" "$scratch_paf"
    zcat "$scratch_paf" > "$input_paf"
    rm -f "$scratch_paf"
    # EB-TP runs the full sweep; DB-TP is capped at mc<=64 (larger b gives very large segments and costly WFA
    # reconstruction on real data). mc=16 is the bootstrap warm-up.
    for cm in edit-distance diagonal-distance; do
        for mc in $MC_LIST; do
            [ "$cm" = "diagonal-distance" ] && [ "$mc" -gt 64 ] && continue
            run_benchmark_real "primates" "$input_paf" "$PRIMATES_SEQS" "$cm" "$dir_base" "$cigzip" "$TMP_DIR" "$mc"
        done
    done
    # FL-TP (fastga): trace spacings $FASTGA_MC_LIST (fastga scaling, like EB/DB-TP). Reconstructs under the same
    # dual gap-affine model as EB/DB-TP (see run_benchmark_real), so it is a full Table 2 entry.
    #for fl_mc in $FASTGA_MC_LIST; do
    #    run_benchmark_real "primates" "$input_paf" "$PRIMATES_SEQS" "fastga" "$dir_base" "$cigzip" "$TMP_DIR" "$fl_mc"
    #done
    rm -f "$input_paf"
done
# Merge Primates results
cat "$TMP_DIR"/*.result.tsv 2>/dev/null | grep "^primates" | sort >> "$REPORT_PRIMATES"
rm -f "$TMP_DIR"/*.result.tsv
# =============================================================================

# Cleanup
rm -rf "$TMP_DIR"
rm -rf $scratch_dir/real_jobs

echo "Real data benchmark complete."
echo "HPRCv2 results: $REPORT_HPRCV2"
echo "Primates results: $REPORT_PRIMATES"
```

Bgzip decompression benchmark:

```shell
for pangenome in hprcv2-25k t2t-ape-pangenome; do
    mkdir -p $dir_base/$pangenome/decompress
    log_dir=$dir_base/$pangenome/decompress

    cd $scratch_dir/$pangenome
    (echo -e "file\tsize_gz_bytes\tdecompress_runtime_sec\tdecompress_memory_kb"; for gz in *.paf.gz; do
        name=$(basename "$gz" .paf.gz)
        log="$log_dir/${name}.bgzip_decompress.log"
        size_gz=$(wc -c < "$gz")

        \time -v bgzip -d -@ $(nproc) -c "$gz" > "/tmp/${name}.paf" 2> "$log"
        rm -f "/tmp/${name}.paf"

        runtime_raw=$(grep "Elapsed (wall clock)" "$log" | sed 's/.*: //')
        runtime=$(echo "$runtime_raw" | awk -F: '{
            if (NF == 3) print $1 * 3600 + $2 * 60 + $3
            else if (NF == 2) print $1 * 60 + $2
            else print $1
        }')
        memory=$(grep "Maximum resident set size" "$log" | sed 's/.*: //')

        echo -e "${name}\t${size_gz}\t${runtime}\t${memory}"
    done) > $pangenome.bgzip-decompression.benchmark.tsv
done
```

## Plotting

All figures are saved directly to `paper/figures/`.

Main figures (compression ratios + decoding costs, from simulated data):

```shell
# Input:  data/simulated-data/benchmark.results.tsv
# Output: compression_ratios.png/.pdf, decoding_cost.png/.pdf
# Deps:   tidyverse, ggplot2, patchwork, cowplot, grid, gridExtra
Rscript scripts/plotting/plot-main-figures.R
```

Supplementary figures S1-S3 (delta effect, tracepoint counts, score distributions):

```shell
# Input:  data/simulated-data/benchmark.results.tsv
#         data/real-data/hprcv2-25k.input-score-frequency.tsv
#         data/real-data/t2t-ape-pangenome.input-score-frequency.by-target.tsv
# Output: figS1_delta_effect.png/.pdf, figS2_tracepoint_counts.png/.pdf,
#         figS3_score_distributions.png/.pdf
# Deps:   tidyverse, ggplot2, patchwork, cowplot, ggh4x
Rscript scripts/plotting/plot-supp-figures.R
```

Supplementary figure S4 (primate pangenome performance across 16 targets):

```shell
# Input:  data/real-data/t2t-ape-pangenome.benchmark.results.tsv
# Output: figS4_primate_performance.png/.pdf
# Deps:   tidyverse, ggplot2, patchwork
Rscript scripts/plotting/plot-supp-figure-primates.R
```

CIGAR dot-plot with tracepoint grids:

```shell
# Input:  PAF file with cg:Z: CIGAR tags (e.g., scripts/toy3mode.paf)
# Output: per-alignment PDFs + combined panel (0000_*_combined.pdf)
# Deps:   matplotlib, numpy; requires cigzip binary
python3 scripts/plotting/cigar_dotplot.py \
  scripts/toy3mode.paf \
  --max-complexity 10 --fastga-spacing 100 \
  --outdir /tmp/cigar_dotplot_out \
  --cigzip $cigzip
```

Compute compression statistics reported in figure captions and text:

```shell
# Input:  data/simulated-data/benchmark.results.tsv
# Output: printed to stdout
# Deps:   tidyverse
Rscript scripts/compute_compression_stats.R
```

## ES-2bit comparison on simulated pangenomes

Compare ES-2bit (drop all `=` ops, store each `I`/`X`/`D` in 2 bits) against FL-TP (l=100), EB-TP (δ=32), and DB-TP (b=32).

```shell
zenodo=/home/guarracino/Desktop/Guarracino/tracepoints/zenodo/simulated-pangenomes
tracepoints_dir=/home/guarracino/git/tracepoints
out_dir=/home/guarracino/Desktop/Guarracino/tracepoints/es2bit-bench
mkdir -p "$out_dir"

cargo build --release --manifest-path "$tracepoints_dir/Cargo.toml" --example es2bit_bench
bench=$tracepoints_dir/target/release/examples/es2bit_bench

for l in 100 1000 10000 100000; do
  for e in 001 005 010 020; do
    prefix="set_${l}_${e}"
    tpa="$zenodo/tpas/${prefix}.standard.edit-distance.128.high.tpa"
    fa_gz="$zenodo/seqs/${prefix}.fa.gz"
    paf="$out_dir/${prefix}.paf"
    fa="$out_dir/${prefix}.fa"

    # TPA → PAF with cg:Z: tags
    cigzip decompress --input "$tpa" --decode --sequence-files "$fa_gz" > "$paf"
    zcat "$fa_gz" > "$fa"

    # Encode/decode/verify all four methods; 1000 records per  PAF
    "$bench" "$paf" "$fa" 1000 > "$out_dir/${prefix}.summary.tsv"
  done
done

# Aggregate summaries.
head -n 1 "$out_dir/set_100_001.summary.tsv" > "$out_dir/all.tsv"
for f in "$out_dir"/set_*.summary.tsv; do tail -n +2 "$f" >> "$out_dir/all.tsv"; done

# Plot and compute stats for the paper text.
Rscript scripts/plotting/plot-es2bit-bench.R "$out_dir/all.tsv" "$out_dir/plots"
Rscript scripts/compute_es2bit_stats.R "$out_dir/all.tsv"

# Copy figures into the paper (S5 = decode time, S6 = storage bits).
cp "$out_dir/plots/es2bit_decode_ns_per_edit.pdf" paper/figures/figS5_es2bit_decode_ns_per_edit.pdf
cp "$out_dir/plots/es2bit_bits_per_edit.pdf" paper/figures/figS6_es2bit_bits_per_edit.pdf
```
