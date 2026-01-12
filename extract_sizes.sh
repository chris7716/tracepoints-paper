#!/usr/bin/env bash

# extract_sizes.sh - Extract file sizes from benchmark reports
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <report1.md> [report2.md] ..."
    echo "Example: $0 tmp_reports/*.md"
    exit 1
fi

SIZES_BYTES="file_sizes_bytes.csv"
SIZES_MB="file_sizes_mb.csv"

# Create CSV headers
echo "Dataset,Record Count,Sequence Length,Error,Tracepoint Type,sample.seq,dataset.fa,fastga.paf,tmp.paf,cigzip.paf,tp.paf,decompressed.paf" > "$SIZES_BYTES"
echo "Dataset,Record Count,Sequence Length,Error,Tracepoint Type,sample.seq,dataset.fa,fastga.paf,tmp.paf,cigzip.paf,tp.paf,decompressed.paf" > "$SIZES_MB"

for report in "$@"; do
    if [[ ! -f "$report" ]]; then
        echo "Skipping $report (not a file)"
        continue
    fi
    
    echo "Processing: $report"
    
    # Extract parameters
    records=$(grep "| record_count (-r) |" "$report" | awk -F'|' '{print $3}' | tr -d ' ')
    length=$(grep "| length (-L) |" "$report" | awk -F'|' '{print $3}' | tr -d ' ')
    error=$(grep "| error (-E) |" "$report" | awk -F'|' '{print $3}' | tr -d ' ')
    n_id=$(grep "| N (-N) |" "$report" | awk -F'|' '{print $3}' | tr -d ' ')
    type=$(grep "| type (-T) |" "$report" | awk -F'|' '{print $3}' | tr -d ' ')
    
    dataset="N${n_id}_R${records}_L${length}_E${error}"
    echo "  Dataset: $dataset"
    
    # Extract file sizes in bytes from the "Artifacts" table
    bytes_section=$(awk '/^## Artifacts$/{flag=1;next}/^## /{flag=0}flag' "$report")
    
    # Extract the data rows (skip header and separator)
    sample_seq_bytes=""
    dataset_fa_bytes=""
    fastga_paf_bytes=""
    tmp_paf_bytes=""
    cigzip_paf_bytes=""
    tp_paf_bytes=""
    decompressed_paf_bytes=""
    
    # Parse bytes table
    while IFS= read -r line; do
        if [[ "$line" =~ sample\.dataset\.[0-9]+\.seq ]]; then
            sample_seq_bytes=$(echo "$line" | awk -F'|' '{print $3}' | tr -d ' ')
        elif [[ "$line" =~ dataset\.[0-9]+\.fa ]]; then
            dataset_fa_bytes=$(echo "$line" | awk -F'|' '{print $3}' | tr -d ' ')
        elif [[ "$line" =~ fastga-datasets.*\.paf ]]; then
            fastga_paf_bytes=$(echo "$line" | awk -F'|' '{print $3}' | tr -d ' ')
        elif [[ "$line" =~ tmp/tmp\..*\.paf ]]; then
            tmp_paf_bytes=$(echo "$line" | awk -F'|' '{print $3}' | tr -d ' ')
        elif [[ "$line" =~ cigzip-datasets.*\.paf ]] && [[ ! "$line" =~ \.tp\. ]] && [[ ! "$line" =~ decompressed ]]; then
            cigzip_paf_bytes=$(echo "$line" | awk -F'|' '{print $3}' | tr -d ' ')
        elif [[ "$line" =~ \.tp\..*\.paf ]] && [[ ! "$line" =~ decompressed ]]; then
            tp_paf_bytes=$(echo "$line" | awk -F'|' '{print $3}' | tr -d ' ')
        elif [[ "$line" =~ decompressed\.paf ]]; then
            decompressed_paf_bytes=$(echo "$line" | awk -F'|' '{print $3}' | tr -d ' ')
        fi
    done <<< "$bytes_section"
    
    # Extract file sizes in MB from the "Artifacts MB" table
    mb_section=$(awk '/^## Artifacts MB$/{flag=1;getline;getline;getline;print;exit}' "$report")
    
    # Parse MB values from the table row
    if [[ -n "$mb_section" ]]; then
        # Split the MB row by | and extract values
        IFS='|' read -ra mb_values <<< "$mb_section"
        
        # Clean up values (remove leading/trailing spaces)
        sample_seq_mb=$(echo "${mb_values[1]}" | tr -d ' ')
        dataset_fa_mb=$(echo "${mb_values[2]}" | tr -d ' ')
        fastga_paf_mb=$(echo "${mb_values[3]}" | tr -d ' ')
        tmp_paf_mb=$(echo "${mb_values[4]}" | tr -d ' ')
        cigzip_paf_mb=$(echo "${mb_values[5]}" | tr -d ' ')
        tp_paf_mb=$(echo "${mb_values[6]}" | tr -d ' ')
        decompressed_paf_mb=$(echo "${mb_values[7]}" | tr -d ' ')
    fi
    
    echo "  Bytes: $sample_seq_bytes, $dataset_fa_bytes, $fastga_paf_bytes, $tmp_paf_bytes, $cigzip_paf_bytes, $tp_paf_bytes, $decompressed_paf_bytes"
    echo "  MB: $sample_seq_mb, $dataset_fa_mb, $fastga_paf_mb, $tmp_paf_mb, $cigzip_paf_mb, $tp_paf_mb, $decompressed_paf_mb"
    
    # Write to CSV files
    echo "$dataset,$records,$length,$error,$type,$sample_seq_bytes,$dataset_fa_bytes,$fastga_paf_bytes,$tmp_paf_bytes,$cigzip_paf_bytes,$tp_paf_bytes,$decompressed_paf_bytes" >> "$SIZES_BYTES"
    echo "$dataset,$records,$length,$error,$type,$sample_seq_mb,$dataset_fa_mb,$fastga_paf_mb,$tmp_paf_mb,$cigzip_paf_mb,$tp_paf_mb,$decompressed_paf_mb" >> "$SIZES_MB"
    
    echo ""
done

echo "=== FILES CREATED ==="
echo "Bytes: $SIZES_BYTES ($(wc -l < "$SIZES_BYTES") records)"
echo "MB: $SIZES_MB ($(wc -l < "$SIZES_MB") records)"

echo ""
echo "=== PREVIEW - BYTES TABLE ==="
head -3 "$SIZES_BYTES" | column -t -s,

echo ""
echo "=== PREVIEW - MB TABLE ==="
head -3 "$SIZES_MB" | column -t -s,
