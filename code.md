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
# To built these two on our HPC
# unset LIBRARY_PATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH
# CC=/usr/bin/gcc CXX=/usr/bin/g++ make BUILD_EXAMPLES=0

out_to_paf=$dir_base/scripts/out-to-paf.py
pafcheck=/moosefs/guarracino/git/pafcheck/target/release/pafcheck

cigzip=$scratch_dir/cigzip/target/release/cigzip

# Use tracepoints-dev branch from https://github.com/AndreaGuarracino/FASTGA.git
fastga_bin=$scratch_dir/FASTGA
FAtoGDB=$fastga_bin/FAtoGDB
PAFtoALN=$fastga_bin/PAFtoALN
ALNtoPAF=$fastga_bin/ALNtoPAF
ONEview=$fastga_bin/ONEview
FASTGA_THREADS=$(nproc)
```

## Simulated data

Generate:

```shell
mkdir -p $dir_base/simulated-data/alns
mkdir -p $dir_base/simulated-data/seqs

num_records=10000

for l in 100 1000 10000 100000; do
  for e in 0.001 0.01 0.05 0.10 0.20; do
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
    $pafcheck --paf $dir_base/simulated-data/alns/$prefix.paf --sequence-files $dir_base/simulated-data/seqs/$prefix.fa --threads $(nproc)
  done
done
```

Simulation - PAF CIGAR -> PAF TRACEPOINTS -> TPA -> PAF TRACEPOINTS -> PAF CIGAR:

```shell
# Run the whole simulation in a subshell so an abort exits only this subshell, NOT the terminal.
(
# Create all output directories
mkdir -p $dir_base/simulated-data/{encode,compress,decompress,decode}
mkdir -p $scratch_dir/benchmark_tmp

REPORT="$dir_base/simulated-data/benchmark.results.tsv"
TMP_DIR="$scratch_dir/benchmark_tmp"

# Strip terminal I/D from each CIGAR because FL-TP cannot represent a terminal indel due to FASTGA's ALNtoPAF requirement
# for every alignment to end in =/X. Trimming up front makes input==decoded for FL-TP and is harmless to EB-TP/DB-TP.
# Empty-after-trim alignments are dropped.
stage_inputs(){
    mkdir -p "$3"
    [ -f "$3/$4.fa" ]  || { cp "$2" "$3/$4.fa.tmp" && mv "$3/$4.fa.tmp" "$3/$4.fa"; }
    [ -f "$3/$4.paf" ] || { awk 'BEGIN{FS=OFS="\t"}
        { ci=0; for(j=13;j<=NF;j++) if($j ~ /^cg:Z:/){ci=j; break}
          if(ci==0){print; next}
          # rev: on - strand the query is reverse-complemented, so a query-consuming (I) op at the START
          # of the CIGAR consumes from the qe end (and at the END consumes from qs). Target (D) is always
          # forward. (Simulated data is all + strand, but keep this identical to the real-data trim.)
          rev=($5=="-"); c=substr($ci,6); qs=$3; qe=$4; ts=$8; te=$9; df=0; db=0; L=length(c)
          while(df<L){ w=substr(c,df+1,20); if(!match(w,/^[0-9]+[ID]/)) break
            op=substr(w,RLENGTH,1); len=substr(w,1,RLENGTH-1)+0
            if(op=="I"){ if(rev) qe-=len; else qs+=len } else ts+=len; df+=RLENGTH }
          while(df+db<L){ last=substr(c,L-db,1); if(last!="I" && last!="D") break
            s=L-db-19; if(s<df+1) s=df+1; w=substr(c,s,L-db-s+1)
            match(w,/[0-9]+[ID]$/); m=substr(w,RSTART); len=substr(m,1,length(m)-1)+0
            if(last=="I"){ if(rev) qs+=len; else qe-=len } else te-=len; db+=length(m) }
          if(df||db){ c=substr(c,df+1,L-df-db); $3=qs;$4=qe;$8=ts;$9=te; $ci="cg:Z:" c }
          if(substr($ci,6)==""){ next }                       # drop any alignment that trims to empty
          print }' "$1" > "$3/$4.paf.tmp" && mv "$3/$4.paf.tmp" "$3/$4.paf"; }
}
export -f stage_inputs

# Precompute CIGAR sizes (fast, run sequentially)
: > "$TMP_DIR/cigar_sizes.txt"
for l in 100 1000 10000 100000; do
  for e in 0.001 0.01 0.05 0.10 0.20; do
    l_nodot=$(echo $l | sed 's/\.//g')
    e_nodot=$(echo $e | sed 's/\.//g')
    prefix=set_${l_nodot}_${e_nodot}
    input_paf="$dir_base/simulated-data/alns/$prefix.paf"

    # Stage to local scratch
    stage_dir="$scratch_dir/staged/$prefix"
    stage_inputs "$input_paf" "$dir_base/simulated-data/seqs/$prefix.fa" "$stage_dir" "$prefix"
    input_paf="$stage_dir/$prefix.paf"
    seq_staged="$stage_dir/$prefix.fa"

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

    # Compress the CIGAR-only PAF with bgzip -l9
    bgzip_scratch="$stage_dir/$prefix.cigar.paf.gz"
    bgzip_file="$dir_base/simulated-data/compress/$prefix.cigar.paf.gz"
    awk -F'\t' '{
        out = $1
        for (i = 2; i <= 12; i++) out = out "\t" $i
        for (i = 13; i <= NF; i++) {
            if ($i ~ /^cg:Z:/) { out = out "\t" $i; break }
        }
        print out
    }' "$input_paf" | bgzip -l9 -c > "$bgzip_scratch"
    size_cigar_bgzip=$(wc -c < "$bgzip_scratch")
    cp -f "$bgzip_scratch" "$bgzip_file"

    # Benchmark bgzip decompression
    bgzip_decomp_log="$dir_base/simulated-data/decompress/$prefix.bgzip_decompress.log"
    \time -v bgzip -d -@ $(nproc) -c "$bgzip_scratch" > $scratch_dir/${prefix}.cigar.paf 2> "$bgzip_decomp_log"
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
            --sequence-files "$seq_staged" \
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
    [[ "$tp_type" == fastga* ]] && cmd_args="" || cmd_args="--complexity-metric $cm"
    # Compression strategy
    if [ "$tp_type" = "fastga" ]; then
        strategy_args="--strategy huffman-nocomp;huffman-nocomp"
    elif [ "$tp_type" = "fastga-no-diff" ]; then
        strategy_args="--strategy huffman-nocomp"
    elif [ "$cm" = "edit-distance" ]; then
        strategy_args="--strategy huffman-nocomp;2d-delta-nocomp"
    else  # diagonal-distance
        strategy_args="--strategy huffman-nocomp;2d-delta-nocomp"
    fi

    # Repeat 10 times for stable averages: standard with mc<=64, or any fastga variant on short sequences (l < 10k). Large fastga jobs are slow, so they run once.
    if [ "$mc" -le 64 ] || { [[ "$tp_type" == fastga* ]] && [ "$l" -lt 10000 ]; }; then
        num_repeats=10
    else
        num_repeats=1
    fi

    # Use unique scratch directory per job.
    [ -n "${scratch_dir:-}" ] || scratch_dir=$(dirname "$TMP_DIR")
    job_scratch="$scratch_dir/jobs/$full_prefix"
    mkdir -p "$job_scratch"

    # Use the local-scratch copies staged in the precompute loop
    stage_dir="$scratch_dir/staged/$prefix"
    input_paf="$stage_dir/$prefix.paf"
    seq_file="$stage_dir/$prefix.fa"

    # File paths (all in scratch)
    tp_paf_full="$job_scratch/$full_prefix.tp.full.paf"      # Full output from encode
    tp_paf="$job_scratch/$full_prefix.tp.paf"                # Subsetted: 12 cols + tp:Z:
    tpa_file="$job_scratch/$full_prefix.tpa"
    decomp_paf="$job_scratch/$full_prefix.decomp.tp.paf"
    decode_paf="$job_scratch/$full_prefix.paf"
    encode_log="$job_scratch/$full_prefix.encode.log"
    compress_log="$job_scratch/$full_prefix.compress.log"
    decompress_log="$job_scratch/$full_prefix.decompress.log"
    decode_log="$job_scratch/$full_prefix.decode.log"

    # Crash guard: if any external tool aborts, stop the whole run.
    guard(){
        echo "ERROR: $1 (exit $2) for $full_prefix; aborting run" >&2
        [ -n "${3:-}" ] && [ -s "${3:-}" ] && { echo "----- $1 log (tail) -----" >&2; tail -40 "$3" >&2; echo "-------------------------" >&2; }
        rm -rf "$job_scratch"; exit 1
    }

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

    # --- FL-TP FASTGA branch: PAFtoALN encode + ALNtoPAF -x decode (FASTGA reference impl) ---
    # Runs in the same loop as the cigzip methods. PAFtoALN hard-codes trace spacing TSPACE=100,
    # so this exists only at mc=100. Maps into the shared schema: .1aln size -> size_tpa_bytes,
    # PAFtoALN -> encode, ALNtoPAF -x -> decode; compress/decompress stages are N/A (0).
    if [ "$tp_type" = "fastga-native" ] || [ "$tp_type" = "fastga-native-nodiff" ]; then
        # 1aln no-diff: PAFtoALN -N omits the per-window diff (X) records (ALNtoPAF -x re-derives them).
        nflag=""; [ "$tp_type" = "fastga-native-nodiff" ] && nflag="-N"

        # Copy it under $prefix.paf so PAFtoALN writes the .1aln next to it.
        cp "$input_paf" "$job_scratch/$prefix.paf"
        $FAtoGDB "$seq_file" "$job_scratch/$prefix.1gdb" > "$job_scratch/gdb.log" 2>&1 || { guard "FAtoGDB" $? "$job_scratch/gdb.log"; }

        # encode (PAFtoALN); both sources are the single .fa GDB, so all names resolve
        enc_rt_sum=0; enc_mem_sum=0
        for _rep in $(seq 1 $num_repeats); do
            rm -f "$job_scratch/$prefix.1aln"
            \time -v $PAFtoALN $nflag -T$FASTGA_THREADS "$job_scratch/$prefix.paf" "$job_scratch/$prefix.1gdb" "$job_scratch/$prefix.1gdb" \
                > /dev/null 2> "$encode_log" || { guard "PAFtoALN" $? "$encode_log"; }
            enc_rt_sum=$(echo "$enc_rt_sum + $(parse_time_log "$encode_log")" | bc -l)
            enc_mem_sum=$((enc_mem_sum + $(parse_memory_log "$encode_log")))
        done
        encode_runtime=$(echo "scale=6; $enc_rt_sum / $num_repeats" | bc -l | sed 's/^\./0./')
        encode_memory=$((enc_mem_sum / num_repeats))
        size_1aln=$(wc -c < "$job_scratch/$prefix.1aln")
        # number of tracepoints in the .1aln: sum of the per-alignment trace-point list lengths
        # (each 'T' line's count). A pure-match alignment is one (0,0) tracepoint, matching cigzip.
        num_1aln_tp=$($ONEview "$job_scratch/$prefix.1aln" 2>/dev/null | awk '$1=="T"{s+=$2} END{printf "%.0f\n", s}')   # %.0f, not print: awk renders >2^31 as %.6g

        # decode (ALNtoPAF -x reconstructs the optimal X/= CIGAR)
        dec_rt_sum=0; dec_mem_sum=0
        for _rep in $(seq 1 $num_repeats); do
            ( cd "$job_scratch" && \time -v $ALNtoPAF -x -T$FASTGA_THREADS "$job_scratch/$prefix.1aln" > "$decode_paf" 2> "$decode_log" ) || { guard "ALNtoPAF" $? "$decode_log"; }
            dec_rt_sum=$(echo "$dec_rt_sum + $(parse_time_log "$decode_log")" | bc -l)
            dec_mem_sum=$((dec_mem_sum + $(parse_memory_log "$decode_log")))
        done
        decode_runtime=$(echo "scale=6; $dec_rt_sum / $num_repeats" | bc -l | sed 's/^\./0./')
        decode_memory=$((dec_mem_sum / num_repeats))

        # correctness: every record preserved AND reconstructed CIGAR consistent with sequences
        n_in=$(grep -c . "$job_scratch/$prefix.paf"); n_out=$(grep -c . "$decode_paf")
        correct="false"
        [ "$n_in" = "$n_out" ] && $pafcheck --paf "$decode_paf" --sequence-files "$seq_file" --threads $(nproc) > /dev/null 2>&1 && correct="true"

        echo -e "$l\t$e\t$tp_type\t$cm\t$mc\t$memory_mode\t$size_cigar\t$size_cigar_bgzip\t$bgzip_decomp_runtime\t$bgzip_decomp_memory\t$align_runtime\t$align_memory\t0\t$num_1aln_tp\t$encode_runtime\t$encode_memory\t$size_1aln\t0\t0\t0\t0\t$correct\t$decode_runtime\t$decode_memory\t$n_out\tNA\tNA\tNA\t$n_in" > "$TMP_DIR/$full_prefix.result.tsv"
        mv "$job_scratch/$prefix.1aln" "$dir_base/simulated-data/compress/$full_prefix.1aln"
        rm -f "$job_scratch/$prefix.paf" "$job_scratch/$prefix.1gdb" "$job_scratch/.$prefix.bps" "$decode_paf"
        echo "Completed: $full_prefix"
        return
    fi

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
            > "$tp_paf_full" 2> "$encode_log" || { guard "encode" $? "$encode_log"; }
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
    } END { printf "%.0f\n", sum }' "$tp_paf")   # %.0f, not print: awk renders >2^31 as %.6g
    
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
            2> "$compress_log" || { guard "compress" $? "$compress_log"; }
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
            --threads $(nproc) \
            2> "$decompress_log" || { guard "decompress" $? "$decompress_log"; }
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
    
    # --- DECODE (banded) ---
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
            -t $(nproc) \
            > "$decode_paf" 2> "$decode_log" || { guard "decode" $? "$decode_log"; }
        _rep_runtime=$(parse_time_log "$decode_log")
        _rep_memory=$(parse_memory_log "$decode_log")
        decode_runtime_sum=$(echo "$decode_runtime_sum + $_rep_runtime" | bc -l)
        decode_memory_sum=$((decode_memory_sum + _rep_memory))
    done
    decode_runtime=$(echo "scale=6; $decode_runtime_sum / $num_repeats" | bc -l | sed 's/^\./0./')
    decode_memory=$((decode_memory_sum / num_repeats))

    # Reconstruction score check
    sc_extract='{ s=""; for(i=13;i<=NF;i++) if($i~/^sc:i:/){s=substr($i,6)+0; break} if(s!="") print $1"|"$6"\t"s }'
    score_orig="$TMP_DIR/orig_score.${l}_${e}.sc"
    if [ ! -s "$score_orig" ]; then
        tmp_orig="$TMP_DIR/orig_score.$full_prefix.tmp.sc"
        $cigzip encode --paf "$input_paf" --type standard --complexity-metric edit-distance \
            --max-complexity 9999999 --distance edit -t $(nproc) 2>/dev/null \
            | awk -F'\t' "$sc_extract" | sort > "$tmp_orig"
        mv -f "$tmp_orig" "$score_orig"
    fi
    score_recon="$job_scratch/$full_prefix.recon.sc"
    $cigzip encode --paf "$decode_paf" --type standard --complexity-metric edit-distance \
        --max-complexity 9999999 --distance edit -t $(nproc) 2>/dev/null \
        | awk -F'\t' "$sc_extract" | sort > "$score_recon"
    read num_output_aln score_identical score_improved score_degraded <<< "$(
        join -t$'\t' "$score_orig" "$score_recon" \
        | awk -F'\t' '{ n++; if($3==$2) id++; else if($3>$2) imp++; else deg++ } END { printf "%.0f %.0f %.0f %.0f\n", n, id, imp, deg }')"
    num_input_aln=$(grep -c . "$input_paf")
    rm -f "$score_recon"

    # --- WRITE RESULT TO TEMP FILE ---
    echo -e "$l\t$e\t$tp_type\t$cm\t$mc\t$memory_mode\t$size_cigar\t$size_cigar_bgzip\t$bgzip_decomp_runtime\t$bgzip_decomp_memory\t$align_runtime\t$align_memory\t$size_tp\t$num_tps\t$encode_runtime\t$encode_memory\t$size_tpa\t$compress_runtime\t$compress_memory\t$decompress_runtime\t$decompress_memory\t$decompress_correct\t$decode_runtime\t$decode_memory\t$num_output_aln\t$score_identical\t$score_improved\t$score_degraded\t$num_input_aln" > "$TMP_DIR/$full_prefix.result.tsv"
    
    # Remove full PAF (keep only subsetted version)
    rm -f "$tp_paf_full"
    
    echo "Completed: $full_prefix"
}

export -f run_benchmark
export scratch_dir   # parallel-spawned workers need these in their environment (job_scratch uses them)
export FAtoGDB PAFtoALN ALNtoPAF ONEview pafcheck FASTGA_THREADS   # used by the fastga-native (FL-TP FASTGA) branch

# Generate all job parameters
> "$TMP_DIR/jobs.txt"
for l in 100 1000 10000 100000; do
  for e in 0.001 0.01 0.05 0.10 0.20; do
    # Look up precomputed values. Read the FIRST matching line only and let `read` split on whitespace,
    # so no embedded newline can leak into a jobs.txt line even if a stale duplicate key ever remained.
    read _ size_cigar size_cigar_bgzip bgzip_decomp_runtime bgzip_decomp_memory align_runtime align_memory \
        < <(grep -m1 "^${l}_${e} " "$TMP_DIR/cigar_sizes.txt")
    if [ -z "$size_cigar" ]; then
        echo "ERROR: no precomputed cigar size for ${l}_${e}; aborting run" >&2; exit 1
    fi

    # Per-condition bootstrap warm-up
    echo "$l $e standard edit-distance 16 high $dir_base $cigzip $size_cigar $size_cigar_bgzip $bgzip_decomp_runtime $bgzip_decomp_memory $align_runtime $align_memory $TMP_DIR" >> "$TMP_DIR/jobs.txt"

    for tp_type in standard fastga fastga-no-diff fastga-native fastga-native-nodiff; do
      if [ "$tp_type" = "fastga" ] || [ "$tp_type" = "fastga-no-diff" ]; then
          cm_list="none"
          mc_list="100 200 300 400 500 1000 1500 2000 2500 3000"
      elif [ "$tp_type" = "fastga-native" ] || [ "$tp_type" = "fastga-native-nodiff" ]; then
          cm_list="none"
          mc_list="100"
      else
          cm_list="edit-distance diagonal-distance"
          mc_list="16 32 48 64 80 96 112 128 256 512 1024"
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

# Run in parallel (adjust -j for number of concurrent jobs).
parallel -j 1 --halt now,fail=1 --colsep ' ' run_benchmark {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} < "$TMP_DIR/jobs.txt" \
    || { echo "ERROR: a benchmark job crashed; aborting run (no results merged)" >&2; exit 1; }

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
echo -e "l\te\ttp_type\tcm\tmc\tmemory_mode\tsize_cigar_bytes\tsize_cigar_bgzip_bytes\tbgzip_decompress_runtime_sec\tbgzip_decompress_memory_kb\talign_runtime_sec\talign_memory_kb\tsize_tp_bytes\tnum_tracepoints\tencode_runtime_sec\tencode_memory_kb\tsize_tpa_bytes\tcompress_runtime_sec\tcompress_memory_kb\tdecompress_runtime_sec\tdecompress_memory_kb\tdecompress_correct\tdecode_runtime_sec\tdecode_memory_kb\tnum_output_alignments\tscore_identical\tscore_improved\tscore_degraded\tnum_input_alignments" > "$REPORT"

# Merge all results into final report (sorted for reproducibility)
cat "$TMP_DIR"/*.result.tsv | sort -t$'\t' -k1,1n -k2,2n -k3,3 -k4,4 -k5,5n -k6,6 >> "$REPORT"

# Cleanup scratch
rm -rf "$TMP_DIR"
rm -rf $scratch_dir/jobs
rm -rf $scratch_dir/staged

echo "Benchmark complete. Results in: $REPORT"
)
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

```

Real indel bias, for comparison with the sweep above (full scan of all 16 primate targets,
566,038 alignments, 683 Gbp; ~30 min wall on 16 cores):

```shell
scan() { f="$1"; b=$(basename "$f" .p70.aln.paf.gz)
  zcat "$f" | awk -F'\t' -v t="$b" '
  {for(i=13;i<=NF;i++) if(substr($i,1,5)=="cg:Z:"){s=substr($i,6);num="";
    for(j=1;j<=length(s);j++){c=substr(s,j,1);
      if(c>="0"&&c<="9")num=num c;
      else{v=num+0; if(c=="I"){Ib+=v;Ie++} else if(c=="D"){Db+=v;De++} num=""}}}}
  END{printf "%s\t%d\t%d\t%.4f\t%d\t%d\t%.4f\n",t,Ib,Db,Ib/Db,Ie,De,Ie/De}'; }
export -f scan
ls $scratch_dir/t2t-ape-pangenome/*.paf.gz | xargs -P 16 -I{} bash -c 'scan "$@"' _ {}

for f in bal skew; do
  paf=$(stat -c%s $f.paf)
  id=$(awk -F'\t' '{for(i=13;i<=NF;i++) if(substr($i,1,5)=="cg:Z:"){s=substr($i,6);num="";
        for(j=1;j<=length(s);j++){c=substr(s,j,1);if(c~/[0-9]/)num=num c;
          else{v=num+0;if(c=="I")I+=v;else if(c=="D")D+=v;num=""}}}} END{printf "%.2f",I/D}' $f.paf)
  echo "=== $f (I/D=$id, PAF=$paf B) ==="
  for m in diagonal-distance edit-distance; do
    $cigzip encode --paf $f.paf --type standard --complexity-metric $m --max-complexity 32 --distance edit -o $f.$m.tp.paf
    tps=$(grep -o 'tp:Z:[^[:space:]]*' $f.$m.tp.paf | awk -F';' '{s+=NF} END{print s}')
    $cigzip compress --input $f.$m.tp.paf --type standard --complexity-metric $m --max-complexity 32 \
      --strategy "huffman-nocomp;2d-delta-nocomp" --output $f.$m.tpa
    tpa=$(stat -c%s $f.$m.tpa)
    awk -v m=$m -v tp=$tps -v t=$tpa -v p=$paf 'BEGIN{printf "  %-18s tps=%d TPA=%d B ratio=%.4f\n",m,tp,t,t/p}'
  done
done
```

### Index size and optional compression layer

Index: one varint offset per alignment record. Encode+compress each primate target, then compare
`.tpa.idx` against `.tpa` (16 targets, 566,038 alignments):

```shell
for m in diagonal-distance edit-distance; do
  for f in $scratch_dir/t2t-ape-pangenome/*.paf.gz; do
    b=$(basename "$f" .p70.aln.paf.gz)
    zcat "$f" | $cigzip encode --paf - --type standard --complexity-metric $m --max-complexity 32 \
        --distance gap-affine2p --penalties 5,8,2,24,1 -o $b.$m.tp.paf
    $cigzip compress -i $b.$m.tp.paf --type standard --complexity-metric $m --max-complexity 32 \
        --strategy "huffman-nocomp;2d-delta-nocomp" -o $b.$m.tpa
    echo "$b $m $(stat -c%s $b.$m.tpa.idx) / $(stat -c%s $b.$m.tpa)"
  done
done
```

Optional per-record compression layer: the `-nocomp` suffix pins it off (and is also the default
when no suffix is given), `-zstd` / `-bgzip` turn it on. Over the 20 simulated conditions:

```shell
for layer in nocomp zstd bgzip; do
  $cigzip compress -i $prefix.tp.paf --type standard --complexity-metric $cm --max-complexity 32 \
     --strategy "huffman-$layer;2d-delta-$layer" -o $prefix.$layer.tpa
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

Length-normalized score distributions (per-base penalty = |score| / alignment block length, PAF column 11):

```shell
# Same pipeline as above, but instead of histogramming the raw score we keep the
# alignment block length (PAF col 11) and emit score/length per alignment, then
# histogram that. Reviewer request: disentangle divergence from alignment length.

# HPRCv2
cd $scratch_dir/hprcv2-25k
find . -name "*.paf.gz" | \
  parallel -j $(nproc) \
    "$cigzip encode --paf {} \
      --type standard --complexity-metric edit-distance \
      --max-complexity 9999999 \
      --distance gap-affine2p --penalties 5,8,2,24,1 \
      -t 1 2>/dev/null \
    | awk -F'\t' '{len=\$11; for(i=13;i<=NF;i++) if(\$i~/^sc:i:/){s=substr(\$i,6)+0; if(len>0) printf \"%.5f\\n\", s/len; next}}'" \
  | sort -n | uniq -c \
  | awk 'BEGIN{OFS="\t"; print "score_per_base","frequency"}{print $2,$1}' \
  > $dir_base/hprcv2-25k/hprcv2-25k.input-scorePerBase-frequency.tsv

# T2T ape pangenome
cd $scratch_dir/t2t-ape-pangenome
find . -name "*.paf.gz" | \
  parallel -j $(nproc) \
    "$cigzip encode --paf {} \
      --type standard --complexity-metric edit-distance \
      --max-complexity 9999999 \
      --distance gap-affine2p --penalties 5,8,2,24,1 \
      -t 1 2>/dev/null \
    | awk -F'\t' '{len=\$11; for(i=13;i<=NF;i++) if(\$i~/^sc:i:/){s=substr(\$i,6)+0; if(len>0) printf \"%.5f\\n\", s/len; next}}'" \
  | sort -n | uniq -c \
  | awk 'BEGIN{OFS="\t"; print "score_per_base","frequency"}{print $2,$1}' \
  > $dir_base/t2t-ape-pangenome/t2t-ape-pangenome.input-scorePerBase-frequency.tsv

# T2T ape pangenome - length-normalized score distribution by target genome
cd $scratch_dir/t2t-ape-pangenome
find . -name "*.paf.gz" | \
  parallel -j $(nproc) \
    "$cigzip encode --paf {} \
      --type standard --complexity-metric edit-distance \
      --max-complexity 9999999 \
      --distance gap-affine2p --penalties 5,8,2,24,1 \
      -t 1 2>/dev/null \
    | awk -F'\t' '{len=\$11; target=\$6; sub(/#.*/, \"\", target); for(i=13;i<=NF;i++) if(\$i~/^sc:i:/){s=substr(\$i,6)+0; if(len>0) printf \"%.5f\\t%s\\n\", s/len, target; next}}'" \
  | sort -t$'\t' -k2,2 -k1,1n | uniq -c \
  | awk 'BEGIN{OFS="\t"; print "score_per_base","frequency","target"}{print $2,$1,$3}' \
  > $dir_base/t2t-ape-pangenome/t2t-ape-pangenome.input-scorePerBase-frequency.by-target.tsv
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
MC_LIST="32 64 128 256"         # EB-TP/DB-TP max-complexity values (mc=16 runs only as the bootstrap warm-up below, excluded from reporting)
FASTGA_MC_LIST="100 500 1000"  # FL-TP trace spacings (fastga scaling)

BENCH_METHODS="ebdb fltptpa fltp1aln"
RUN_TAG="all"

want_method(){ case " $BENCH_METHODS " in *" $1 "*) return 0;; *) return 1;; esac; }   # want_method <token>
TMP_DIR="$scratch_dir/real_benchmark_tmp"
unset _SCAN_CACHE; declare -A _SCAN_CACHE   # input PAF -> "<cigar_bytes> <alignment_count>", filled once per file
# unset first: `declare -A` on an existing array KEEPS its entries, so re-pasting this block into a
# shell that already ran would reuse stale cached values instead of rescanning.
mkdir -p $TMP_DIR
HEADER="dataset\tpaf_file\tcm\tmc\tdecode_mode\tsize_cigar_bytes\tsize_tp_bytes\tnum_tracepoints\tsize_tpa_bytes\tencode_runtime_sec\tencode_memory_kb\tcompress_runtime_sec\tcompress_memory_kb\tdecompress_runtime_sec\tdecompress_memory_kb\tdecode_runtime_sec\tdecode_memory_kb\tnum_output_alignments\tscore_identical\tscore_improved\tscore_degraded\tnum_input_alignments\tdistance"

# Resume helper: a per-iteration result row "went fine" iff it has all 23 columns, the decode ran
# (decode_runtime col16 != NA), and it has a valid output count.
iteration_ok() {
    awk -F'\t' 'END{
        if (NF!=23) exit 1
        if ($16=="" || $16=="NA") exit 1                       # decode ran
        if ($18=="" || $18=="NA" || $18+0<=0) exit 1           # has output
        if ($22=="" || $22=="NA") exit 1
        if ($3=="edit-distance" || $3=="diagonal-distance") { if ($18+0 != $22+0) exit 1 }  # EB/DB: exact
        else { if ($18+0 < $22+0) exit 1 }                     # FL-TP: gap-splits allow more outputs
        exit 0
    }' "$1"
}

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
    local DIST="${9:-gap-affine2p}"   # reconstruction distance model: edit | gap-affine2p (recorded in the output)

    paf_name=$(basename "$input_paf" .paf.gz)
    paf_name=$(basename "$paf_name" .paf)  # Handle both .paf.gz and .paf

    prefix="$paf_name.$cm.$MC.$DIST"
    job_scratch="$scratch_dir/real_jobs/$prefix"

    # --- Resume ---
    # PTe mc=16 bootstrap warm-up is never checkpointed, so every re-run repeats it first
    ckpt_dir="$TMP_DIR/$(echo "$BENCH_METHODS" | tr -s ' ' '-')${RUN_TAG:+-$RUN_TAG}"; mkdir -p "$ckpt_dir"
    ckpt="$ckpt_dir/$prefix.result.tsv"
    if [ "$MC" != "16" ] && [ -s "$ckpt" ] && iteration_ok "$ckpt"; then
        echo "SKIP (already done): $dataset/$paf_name ($cm, mc=$MC)"
        return 0
    fi
    [ "$MC" != "16" ] && rm -f "$ckpt"   # clear a stale NA/partial checkpoint before redoing
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

    # CIGAR size (12 cols + cg:Z: only, excludes other tags for fair comparison) and the input alignment count
    if [ -z "${_SCAN_CACHE[$input_paf]:-}" ]; then
        _SCAN_CACHE[$input_paf]=$(awk -F'\t' 'BEGIN{OFS="\t"}
            { out=$1; for(i=2;i<=12;i++) out=out OFS $i
              for(i=13;i<=NF;i++) if($i ~ /^cg:Z:/){ out=out OFS $i; break }
              b += length(out)+1; n++ }
            END{ printf "%.0f %.0f\n", b, n }' "$input_paf")
    fi
    read -r size_cigar num_input_aln <<< "${_SCAN_CACHE[$input_paf]}"

    # Skip empty PAF files
    if [[ "$size_cigar" -lt 10 ]]; then
        echo "Skipping empty PAF: $dataset/$paf_name ($cm, mc=$MC)"
        rm -rf "$job_scratch"
        return 0
    fi
    # PAFtoALN reuses the per-dataset FASTGA GDB staged before the main loop ($scratch_dir/<dataset>-gdb/<dataset>.1gdb).
    if [ "$cm" = "fastga-native" ] || [ "$cm" = "fastga-native-nodiff" ]; then
        local gdb="$scratch_dir/${dataset}-gdb/${dataset}.1gdb"
        # FL-TP 1aln reconstruction (ALNtoPAF) is edit-based: score orig/recon under edit-distance to match.
        local dist_args="--distance edit"
        local nflag=""; [ "$cm" = "fastga-native-nodiff" ] && nflag="-N"
        local fixed_paf="$job_scratch/$prefix.fixed.paf"

        # Preserve failure logs before the per-job dir is rm'd, so PAFtoALN/ALNtoPAF failures stay
        # debuggable (kept on moosefs, survives the scratch cleanup). Rerun them with
        # scripts/rerun_fltp1aln_failures.sh. Without this, the only trace of a failure is decode=NA.
        local fail_log_dir="$dir_base/failed_logs/$prefix"
        keep_fail_logs(){ mkdir -p "$fail_log_dir"; for L in "$encode_log" "$decode_log"; do [ -s "$L" ] && cp "$L" "$fail_log_dir/"; done; echo "  (kept failure logs in $fail_log_dir)" >&2; }

        # PAFtoALN names the output after the INPUT basename, i.e. "$prefix.fixed.1aln"
        # (not "$prefix.1aln"). Reading the wrong name silently yielded empty size/0 tps.
        local aln_file="$job_scratch/$prefix.fixed.1aln"

        if [ ! -f "$gdb" ]; then
            echo "ERROR: FASTGA GDB not found at $gdb (build it before the loop)" >&2
            rm -rf "$job_scratch"; return 1
        fi

        # The input PAF is already terminal-indel-trimmed and degenerate-filtered by the caller
        ln -sf "$input_paf" "$fixed_paf"

        # encode: PAFtoALN -> .1aln (source and target GDB are the same combined dataset GDB).
        # Retry: PAFtoALN -T96 has a rare, load-dependent heap-corruption crash (SIGABRT,
        # "corrupted double-linked list") that clears on re-run; a few attempts make it a non-event.
        pafaln_ok=0
        for _try in 1 2 3; do
            \time -v $PAFtoALN $nflag -T$FASTGA_THREADS "$fixed_paf" "$gdb" "$gdb" > /dev/null 2> "$encode_log" \
                && { pafaln_ok=1; break; }
            echo "WARN: PAFtoALN attempt $_try failed for $dataset/$paf_name (mc=$MC): $(tail -1 "$encode_log")" >&2
        done
        [ "$pafaln_ok" -eq 1 ] \
            || { echo "ERROR: PAFtoALN failed for $dataset/$paf_name (mc=$MC) after 3 attempts: $(tail -1 "$encode_log")" >&2; keep_fail_logs; rm -rf "$job_scratch"; return 1; }
        encode_runtime=$(parse_time_log "$encode_log"); encode_memory=$(parse_memory_log "$encode_log")
        size_1aln=$(wc -c < "$aln_file")
        num_tps=$($ONEview "$aln_file" 2>/dev/null | awk '$1=="T"{s+=$2} END{printf "%.0f\n", s}')   # %.0f, not print: awk renders >2^31 as %.6g

        # decode: ALNtoPAF -x reconstructs the optimal =/X CIGAR.
        # Retry: like PAFtoALN, ALNtoPAF -x -T96 can hit a rare, load-dependent transient failure
        # that clears on re-run; a few attempts make it a non-event. Any residual deterministic
        # failure (should be none after the FASTGA fork fixes) still falls through to decode=NA.
        decode_status=1
        for _try in 1 2 3; do
            decode_status=0
            ( cd "$job_scratch" && \time -v $ALNtoPAF -x -T$FASTGA_THREADS "$aln_file" > "$decode_paf" 2> "$decode_log" ) || decode_status=$?
            [ "$decode_status" -eq 0 ] && break
            echo "WARN: ALNtoPAF -x attempt $_try failed (exit $decode_status) for $dataset/$paf_name (mc=$MC): $(tail -1 "$decode_log")" >&2
        done
        if [ "$decode_status" -eq 0 ]; then
            decode_runtime=$(parse_time_log "$decode_log"); decode_memory=$(parse_memory_log "$decode_log")
        else
            echo "WARN: ALNtoPAF -x failed (exit $decode_status) for $dataset/$paf_name (mc=$MC) after 3 attempts: $(tail -1 "$decode_log"); recording size only (decode=NA)" >&2
            decode_runtime=NA; decode_memory=NA
            keep_fail_logs   # preserve encode.log/decode.log before the per-job rm below
        fi

        # Score-preservation is not evaluated for FL-TP 1aln.
        num_aln=$(grep -c . "$fixed_paf")
        score_identical=NA; score_improved=NA; score_degraded=NA

        # No separate tp stream or compression: the .1aln is the stored file.
        echo -e "$dataset\t$paf_name\t$cm\t$MC\theuristic\t$size_cigar\t$size_1aln\t$num_tps\t$size_1aln\t$encode_runtime\t$encode_memory\t0\t0\t0\t0\t$decode_runtime\t$decode_memory\t$num_aln\t$score_identical\t$score_improved\t$score_degraded\t$num_input_aln\t$DIST" > "$ckpt"

        out_dir="$dir_base/$dataset"; mkdir -p "$out_dir/encode" "$out_dir/decode"
        [ -s "$encode_log" ] && cp "$encode_log" "$out_dir/encode/$prefix.log"
        [ -s "$decode_log" ] && cp "$decode_log" "$out_dir/decode/$prefix.log"
        rm -rf "$job_scratch"
        echo "Completed (1aln$nflag): $dataset/$paf_name (mc=$MC)"
        return 0
    fi

    if [ "$DIST" = "gap-affine2p" ]; then
        dist_args="--distance gap-affine2p --penalties 5,8,2,24,1"
    else
        dist_args="--distance edit"
    fi
    fastga_contigs_args=""
    if [ "$cm" = "fastga" ]; then
        tp_type="fastga";  cm_args="";                         strategy_args="--strategy huffman-nocomp;huffman-nocomp"
        fastga_contigs_args="--fastga-contigs $scratch_dir/${dataset}-gdb/${dataset}.contigs.tsv"
    elif [ "$cm" = "fastga-no-diff" ]; then
        tp_type="fastga-no-diff"; cm_args="";                  strategy_args="--strategy huffman-nocomp"
        fastga_contigs_args="--fastga-contigs $scratch_dir/${dataset}-gdb/${dataset}.contigs.tsv"
    elif [ "$cm" = "edit-distance" ]; then
        tp_type="standard"; cm_args="--complexity-metric $cm"; strategy_args="--strategy huffman-nocomp;2d-delta-nocomp"
    else  # diagonal-distance
        tp_type="standard"; cm_args="--complexity-metric $cm"; strategy_args="--strategy huffman-nocomp;2d-delta-nocomp"
    fi

    mem_mode="adaptive"

    # --- ENCODE ---
    \time -v $cigzip encode \
        --paf "$input_paf" \
        --type $tp_type \
        $cm_args \
        $fastga_contigs_args \
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
    } END { printf "%.0f\n", sum }' "$tp_paf")   # %.0f, not print: awk renders >2^31 as %.6g

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

    # --- DECOMPRESS ---
    \time -v $cigzip decompress \
        --input "$tpa_file" \
        --output "$decomp_paf" \
        --threads $(nproc) \
        2> "$decompress_log"

    decompress_runtime=$(parse_time_log "$decompress_log")
    decompress_memory=$(parse_memory_log "$decompress_log")

    # --- DECODE (banded, default) ---
    local gdb_seq="$scratch_dir/${dataset}-gdb/${dataset}.1gdb"
    decode_status=0
    \time -v $cigzip decode \
        --paf "$decomp_paf" \
        --sequence-files "$gdb_seq" \
        --type $tp_type \
        $cm_args \
        --max-complexity $MC \
        $dist_args \
        --memory-mode $mem_mode \
        -t $(nproc) \
        > "$decode_paf" 2> "$decode_log" || decode_status=$?

    decode_runtime=$(parse_time_log "$decode_log")
    decode_memory=$(parse_memory_log "$decode_log")
    rm -f "$decomp_paf"

    if [ "$decode_status" -ne 0 ]; then
        out_dir="$dir_base/$dataset"
        mkdir -p "$out_dir/decode"
        [ -f "$decode_log" ] && cp "$decode_log" "$out_dir/decode/$prefix.FAILED.log"
        echo "ERROR: cigzip decode failed (exit $decode_status) for $dataset/$paf_name ($cm, mc=$MC); recording decode=NA (log: $out_dir/decode/$prefix.FAILED.log)" >&2
        echo -e "$dataset\t$paf_name\t$cm\t$MC\theuristic\t$size_cigar\t$size_tp\t$num_tps\t$size_tpa\t$encode_runtime\t$encode_memory\t$compress_runtime\t$compress_memory\t$decompress_runtime\t$decompress_memory\tNA\tNA\tNA\tNA\tNA\tNA\t$num_input_aln\t$DIST" > "$ckpt"
        rm -rf "$job_scratch"
        return 0
    fi

    # Score preservation is not a guaranteed property and is not reported for FL-TP.
    if [ "$cm" = "fastga" ] || [ "$cm" = "fastga-no-diff" ]; then
        num_aln=$(grep -c . "$decode_paf")
        score_identical=NA; score_improved=NA; score_degraded=NA
        rm -f "$tp_paf_full" "$decode_paf"
    else
    # Compare scores: original (from tp_paf_full) vs reconstructed (from decode_paf)
    # Sort both files by first 9 columns to ensure matching alignment order
    # (multi-threaded output may not preserve input order)
    # Lower score = better alignment (penalties are negative)
    # Also save original input alignments (with CIGAR) that lead to degraded scores
    sorted_input="$job_scratch/sorted_input.paf"
    sorted_encoded="$job_scratch/sorted_encoded.paf"
    sorted_decode="$job_scratch/sorted_decode.paf"

    # Sort all three by the first 9 cols so paste lines up the same alignment across files.
    sort -T /scratch --parallel=$(nproc) -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 -k7,7n -k8,8n -k9,9n "$input_paf" > "$sorted_input" &
    sort -T /scratch --parallel=$(nproc) -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 -k7,7n -k8,8n -k9,9n "$tp_paf_full" > "$sorted_encoded" &
    sort -T /scratch --parallel=$(nproc) -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 -k7,7n -k8,8n -k9,9n "$decode_paf" > "$sorted_decode" &
    wait
    rm -f "$tp_paf_full" "$decode_paf"

    read num_aln score_identical score_improved score_degraded <<< $(
        paste "$sorted_input" "$sorted_encoded" "$sorted_decode" | \
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
                # Output: the 12 mandatory PAF cols (identity + coords) + orig/recon score.
                # The cg:Z: CIGAR is deliberately NOT written: FL-TP degrades many records, each
                # carrying a whole-chromosome CIGAR (megabytes), which balloons this file to tens of
                # GB and dominates benchmark I/O. Identity+scores keep the diagnostic (which
                # alignments degraded, by how much); the CIGAR is recoverable from the input PAF.
                out = $1
                for (i = 2; i <= 12; i++) out = out "\t" $i
                print out "\t" orig_score "\t" recon_score >> degraded_file
            }
        }
        END { printf "%.0f %.0f %.0f %.0f\n", n, identical, improved, degraded }'
    )
    rm -f "$sorted_input" "$sorted_encoded" "$sorted_decode"
    fi

    # Write result
    echo -e "$dataset\t$paf_name\t$cm\t$MC\theuristic\t$size_cigar\t$size_tp\t$num_tps\t$size_tpa\t$encode_runtime\t$encode_memory\t$compress_runtime\t$compress_memory\t$decompress_runtime\t$decompress_memory\t$decode_runtime\t$decode_memory\t$num_aln\t$score_identical\t$score_improved\t$score_degraded\t$num_input_aln\t$DIST" > "$ckpt"

    # Determine output directory and move files immediately ($dataset is the output dir name)
    out_dir="$dir_base/$dataset"

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

trim_terminal_indels(){
    awk 'BEGIN{FS=OFS="\t"}
        { ci=0; for(j=13;j<=NF;j++) if($j ~ /^cg:Z:/){ci=j; break}
          if(ci==0){print; next}
          # rev: on - strand the query is reverse-complemented, so a query-consuming (I) op at the START
          # of the CIGAR consumes from the qe end (and at the END consumes from qs). Target (D) is always
          # forward. Getting this wrong shifts coords and breaks ~half of real alignments (pafcheck fail).
          rev=($5=="-"); c=substr($ci,6); qs=$3; qe=$4; ts=$8; te=$9; df=0; db=0; L=length(c)
          while(df<L){ w=substr(c,df+1,20); if(!match(w,/^[0-9]+[ID]/)) break
            op=substr(w,RLENGTH,1); len=substr(w,1,RLENGTH-1)+0
            if(op=="I"){ if(rev) qe-=len; else qs+=len } else ts+=len; df+=RLENGTH }
          while(df+db<L){ last=substr(c,L-db,1); if(last!="I" && last!="D") break
            s=L-db-19; if(s<df+1) s=df+1; w=substr(c,s,L-db-s+1)
            match(w,/[0-9]+[ID]$/); m=substr(w,RSTART); len=substr(m,1,length(m)-1)+0
            if(last=="I"){ if(rev) qs+=len; else qe-=len } else te-=len; db+=length(m) }
          if(df||db){ c=substr(c,df+1,L-df-db); $3=qs;$4=qe;$8=ts;$9=te; $ci="cg:Z:" c }
          if(substr($ci,6)==""){ next }                       # drop any alignment that trims to empty
          print }' "$1"
}

# Create jobs directory
mkdir -p $scratch_dir/real_jobs

# =============================================================================
# Real-data benchmark (HPRCv2 human, then T2T-ape primates)
# =============================================================================
run_real_dataset() {
    local dataset="$1"
    local pafs_dir="$dir_base/${dataset}-pafs"
    local fasta_dir="$scratch_dir/${dataset}-fasta"
    local gdb_dir="$scratch_dir/${dataset}-gdb"      # GDB: $gdb_dir/<ds>.{1gdb,1ano,contigs.tsv} + hidden .<ds>.bps
    local seqs="$TMP_DIR/${dataset}.seqlist.txt"
    mkdir -p "$dir_base/$dataset"/{encode,compress,decompress,decode} "$gdb_dir"

    ls "$fasta_dir"/*.fa > "$seqs" 2>/dev/null || : > "$seqs"

    # To build the GDB on a separate big-scratch node and ship it here, use scripts/build-fastga-gdb.sh.
    if [ ! -f "$gdb_dir/${dataset}.1gdb" ]; then
        cat $(cat "$seqs") | bgzip -@ $(nproc) -l 1 -c > "$gdb_dir/${dataset}.combined.fa.gz"
        $FAtoGDB "$gdb_dir/${dataset}.combined.fa.gz" "$gdb_dir/${dataset}.1gdb" > "$gdb_dir/${dataset}.gdb.log" 2>&1
        rm -f "$gdb_dir/${dataset}.combined.fa.gz"
    fi
    
    # Contig table for cigzip FL-TP TPA (built from the GDB.
    if [ ! -f "$gdb_dir/${dataset}.contigs.tsv" ]; then
        $ONEview "$gdb_dir/${dataset}.1gdb" | awk '$1=="S"{n=$3;p=0;next} $1=="C"{print n"\t"p"\t"p+$2;p+=$2;next} $1=="G"{p+=$2}' > "$gdb_dir/${dataset}.contigs.tsv"
    fi

    # Output table. Suffix = selection tokens joined by '-' (.ebdb / .fltp1aln / .ebdb-fltptpa-fltp1aln ...)
    local suffix=$(echo $BENCH_METHODS | tr -s ' ' '-')${RUN_TAG:+-$RUN_TAG}
    local report="$dir_base/${dataset}.benchmark.results.${suffix}.tsv"
    # Report is rebuilt atomically from the persisted checkpoints AFTER the loop (resume-safe: an
    # interrupted run leaves the previous report intact). Do not truncate it here.
    echo "Running $dataset benchmark..."
    # Process one PAF.gz at a time: copy to /scratch, decompress, run the methods, clean up.
    # Multi-node: set SHARD_N=<nodes> and SHARD_ID=0..N-1 to take every N-th target. Shards are
    # disjoint, so each node keeps its own checkpoints; merge them onto one node at the end and
    # re-run there to rebuild the full report. Unset SHARD_N = one node does everything.
    local _shard_i=0
    # LC_ALL=C sort, not the bare glob: glob order follows the locale's collation (C puts GCA_ first,
    # en_US puts chm13 first), so nodes with different locales would build different shards and
    # silently double-run some targets while skipping others. Target names contain no spaces.
    for paf in $(printf '%s\n' "$pafs_dir"/*.paf.gz | LC_ALL=C sort); do
        [ -f "$paf" ] || continue
        if [ -n "${SHARD_N:-}" ]; then
            _shard_i=$((_shard_i + 1))
            [ $(( (_shard_i - 1) % SHARD_N )) -eq "${SHARD_ID:-0}" ] || continue
        fi
        local paf_basename=$(basename "$paf")
        local scratch_paf="$scratch_dir/$paf_basename"
        local input_paf="$scratch_dir/${paf_basename%.paf.gz}.paf"
        # A finished target is marked here, so a resume skips it before the copy+zcat+prep below.
        local ckpt_dir="$TMP_DIR/$(echo "$BENCH_METHODS" | tr -s ' ' '-')${RUN_TAG:+-$RUN_TAG}"
        local target_done="$ckpt_dir/.${paf_basename%.paf.gz}.done"
        if [ -f "$target_done" ]; then
            echo "SKIP (target already complete): $dataset/${paf_basename%.paf.gz}"
            continue
        fi
        cp "$paf" "$scratch_paf"; zcat "$scratch_paf" > "$input_paf"; rm -f "$scratch_paf"
        # --- Input prep shared across methods (no genome excluded; just preprocessing) ---
        # (1) Drop degenerate (zero query/target span = pure-indel) alignments for the FL-TP input only:
        #     FL-TP cannot represent them.
        local nodegen_paf="${input_paf%.paf}.nodegen.paf"
        awk -F'\t' '{ for(i=13;i<=NF;i++) if($i ~ /^cg:Z:/){ if($3==$4 || $8==$9) next; break } print }' \
            "$input_paf" > "$nodegen_paf"
        # (2) FL-TP input: terminal-indel trim so every alignment ends in =/X (ALNtoPAF -x requires it).
        local fltp_paf="${input_paf%.paf}.fltp.paf"
        trim_terminal_indels "$nodegen_paf" > "$fltp_paf"

        # --- Bootstrap warm-up: EB-TP mc=16. Runs FIRST to warm the FS cache / stabilize timings; excluded
        #     from reporting. Needs FASTAs, so only when a FASTA-using method is selected; fltp1aln-only
        #     (GDB-only node) has no warm-up. mc=16 is NOT in MC_LIST, so this is the sole mc=16 row. ---
        if want_method ebdb || want_method fltptpa || want_method fltptpa_diff || want_method fltptpa_nodiff; then
            run_benchmark_real "$dataset" "$input_paf" "$seqs" "edit-distance" "$dir_base" "$cigzip" "$TMP_DIR" 16 "gap-affine2p"
        fi

        # ================= EB-TP / DB-TP (cigzip) -- needs FASTAs, not GDB ================================
        # Reconstruct under gap-affine2p (matching the real WFMASH scoring).
        if want_method ebdb; then
            for cm in edit-distance diagonal-distance; do
                for mc in $MC_LIST; do
                    run_benchmark_real "$dataset" "$input_paf" "$seqs" "$cm" "$dir_base" "$cigzip" "$TMP_DIR" "$mc" "gap-affine2p"
                done
            done
        fi
        # ================= end EB-TP / DB-TP =============================================================

        # ================= FL-TP TPA (cigzip fastga / fastga-no-diff) -- needs FASTAs, not GDB ===========
        # gap-affine2p at every spacing (matches WFMASH scoring, like EB/DB). edit reconstruction ONLY
        # at spacing 100, where it is compared against FL-TP 1aln (ALNtoPAF -x, edit / l=100 only).
        local fl_diff=false fl_nodiff=false
        if want_method fltptpa || want_method fltptpa_diff;   then fl_diff=true;   fi
        if want_method fltptpa || want_method fltptpa_nodiff; then fl_nodiff=true; fi
        if $fl_diff || $fl_nodiff; then
            for fl_mc in $FASTGA_MC_LIST; do
                if $fl_diff;   then run_benchmark_real "$dataset" "$fltp_paf" "$seqs" "fastga"         "$dir_base" "$cigzip" "$TMP_DIR" "$fl_mc" "gap-affine2p"; fi
                if $fl_nodiff; then run_benchmark_real "$dataset" "$fltp_paf" "$seqs" "fastga-no-diff" "$dir_base" "$cigzip" "$TMP_DIR" "$fl_mc" "gap-affine2p"; fi
            done

            if $fl_diff;   then run_benchmark_real "$dataset" "$fltp_paf" "$seqs" "fastga"         "$dir_base" "$cigzip" "$TMP_DIR" 100 "edit"; fi
            if $fl_nodiff; then run_benchmark_real "$dataset" "$fltp_paf" "$seqs" "fastga-no-diff" "$dir_base" "$cigzip" "$TMP_DIR" 100 "edit"; fi
        fi
        # ================= end FL-TP TPA =================================================================

        # ================= FL-TP 1aln (FASTGA PAFtoALN/ALNtoPAF reference) -- needs GDB, not FASTAs =======
        # l=100; ALNtoPAF -x reconstructs under edit distance.
        if want_method fltp1aln || want_method fltp1aln_diff;   then run_benchmark_real "$dataset" "$fltp_paf" "$seqs" "fastga-native"         "$dir_base" "$cigzip" "$TMP_DIR" 100 "edit"; fi
        if want_method fltp1aln || want_method fltp1aln_nodiff; then run_benchmark_real "$dataset" "$fltp_paf" "$seqs" "fastga-native-nodiff" "$dir_base" "$cigzip" "$TMP_DIR" 100 "edit"; fi
        # ================= end FL-TP 1aln ================================================================

        mkdir -p "$ckpt_dir" && touch "$target_done"
        rm -f "$input_paf" "$nodegen_paf" "$fltp_paf"
    done
    # Rebuild the report atomically from the persisted per-iteration checkpoints (this run's new rows
    # plus any carried over from prior/interrupted runs). Checkpoints are KEPT so a re-run resumes.
    # .$(hostname) in the temp name: with SHARD_N the nodes all write this same moosefs path, and a
    # shared "$report.tmp" would let two finishing nodes interleave into one another's file.
    local report_tmp="$report.tmp.$(hostname)"
    { echo -e "$HEADER"; cat "$TMP_DIR/$suffix"/*.result.tsv 2>/dev/null | grep "^$dataset" | sort; } > "$report_tmp" && mv "$report_tmp" "$report"
}

# Select methods via BENCH_METHODS
run_real_dataset "hprcv2-25k"
run_real_dataset "t2t-ape-pangenome"
# =============================================================================

# Cleanup. NOTE: $TMP_DIR holds the resume checkpoints (per-iteration result rows). Keep it if you
# might re-run to fill in failed/interrupted iterations. Delete it ONLY to force a full recompute:
#   rm -rf "$TMP_DIR"
rm -rf $scratch_dir/real_jobs

suffix=$(echo $BENCH_METHODS | tr -s ' ' '-')${RUN_TAG:+-$RUN_TAG}
echo "Real data benchmark complete (BENCH_METHODS='$BENCH_METHODS'${RUN_TAG:+, RUN_TAG='$RUN_TAG'})."
echo "  hprcv2-25k:        $dir_base/hprcv2-25k.benchmark.results.${suffix}.tsv"
echo "  t2t-ape-pangenome: $dir_base/t2t-ape-pangenome.benchmark.results.${suffix}.tsv"
```

Assemble the final table from whatever suffixed partials exist into `<ds>.benchmark.results.tsv`.

```bash
dir_base=/moosefs/guarracino/tracepoints
for ds in hprcv2-25k t2t-ape-pangenome; do
    out="$dir_base/$ds.benchmark.results.tsv"
    partials=$(ls "$dir_base/$ds".benchmark.results.*.tsv 2>/dev/null)   # excludes canonical (no middle token)
    [ -n "$partials" ] || { echo "skip $ds: no partials found" >&2; continue; }
    head -1 $(echo "$partials" | head -1) > "$out"                       # header
    tail -q -n +2 $partials | sort >> "$out"                            # bodies from all selections, sorted
    echo "merged -> $out ($(( $(wc -l < "$out") - 1 )) rows) from: $(basename -a $partials | tr '\n' ' ')"
done
```

Table 2 (`tab:comprehensive`) rows come from this merged `<ds>.benchmark.results.tsv`:
`scripts/gen_table2_rows.sh <ds>.benchmark.results.tsv <uncompressed_PAF_GiB>`.

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

Space-time tradeoff (stored size vs reconstruction time, all methods across their parameter sweeps):

```shell
# Input:  data/simulated-data/benchmark.results.tsv
#   FL-TP (fastga l 100..2000), FL-TP FASTGA (fastga-native l=100), EB-TP/DB-TP (32..1024)
# Output: space_time_tradeoff.png/.pdf
# Deps:   tidyverse, ggplot2, scales
Rscript scripts/plotting/plot-tradeoff.R
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

# WFA2-lib tools + PAF converter, used to synthesise the 0.1% set (see below).
generate_dataset=/moosefs/guarracino/git/WFA2-lib/bin/generate_dataset
align_benchmark=/moosefs/guarracino/git/WFA2-lib/bin/align_benchmark
out_to_paf=$dir_base/scripts/out-to-paf.py

for l in 100 1000 10000 100000; do
  prefix="set_${l}_0001"
  $generate_dataset -n 1000 -l "$l" -e 0.001 -o "$out_dir/${prefix}.seq"
  $align_benchmark -i "$out_dir/${prefix}.seq" -a edit-wfa --wfa-memory high --output "$out_dir/${prefix}.out" -t 1
  python3 "$out_to_paf" "$out_dir/${prefix}.out" "$out_dir/${prefix}.paf"
  awk 'BEGIN{c=1}
       /^>/{print ">target_" c; print substr($0,2)}
       /^</{print ">query_" c; print substr($0,2); c++}' \
    "$out_dir/${prefix}.seq" > "$out_dir/${prefix}.fa"
  "$bench" "$out_dir/${prefix}.paf" "$out_dir/${prefix}.fa" 1000 > "$out_dir/${prefix}.summary.tsv"
done

# 1/5/10/20% from the published zenodo TPAs.
for l in 100 1000 10000 100000; do
  for e in 001 005 010 020; do
    prefix="set_${l}_${e}"
    tpa="$zenodo/tpas/${prefix}.standard.edit-distance.128.high.tpa"
    fa_gz="$zenodo/seqs/${prefix}.fa.gz"
    paf="$out_dir/${prefix}.paf"
    fa="$out_dir/${prefix}.fa"

    # TPA → PAF with cg:Z: tags
    cigzip decompress --input "$tpa" --decode --sequence-files "$fa_gz" --threads $(nproc) > "$paf"
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
