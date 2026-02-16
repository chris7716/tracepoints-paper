#!/usr/bin/env bash
# filepath: /Users/hasitha/Documents/biology/tracepoints-paper/run_paf_benchmark.sh
set -euo pipefail

# Runs cigzip encode/decode benchmarks on an existing PAF file
# Writes a Markdown report that includes:
#  - Input params
#  - Each command (verbatim)
#  - Output of each sub-command (stdout + stderr), including /usr/bin/time -v
#
# Usage:
#   ./run_paf_benchmark.sh -p <paf_file> -f <fasta_file> -N <N> [-C <complexity_metric>] [-M <max_complexity>] [-T <type>] [-t <threads>] [-o <output_dir>]

# ---------- Defaults / paths ----------
CM="diagonal-distance"  # Default complexity metric
MC="100"               # Default max complexity
TYPE="standard"        # Default tracepoint type
THREADS="8"            # Default threads
OUTPUT_DIR="data/benchmarks"  # Default output directory

# Tool paths (override via env if needed)
CIGZIP_BIN="${CIGZIP_BIN:-../cigzip/target/debug/cigzip}"

# Prefer /usr/bin/time -v if available
TIME_BIN="$(command -v /usr/bin/time >/dev/null 2>&1 && echo '/usr/bin/time -v' || echo time)"

# ---------- Args ----------
PAF_FILE="" FASTA_FILE="" N_ID=""

usage() {
  echo "Usage: $0 -p <paf_file> -f <fasta_file> -N <N> [-C <complexity_metric>] [-M <max_complexity>] [-T <type>] [-t <threads>] [-o <output_dir>]"
  echo ""
  echo "Required parameters:"
  echo "  -p <paf_file>          Input PAF file to benchmark"
  echo "  -f <fasta_file>        Corresponding FASTA file (for sequence data)"
  echo "  -N <N>                 Benchmark ID number (for file naming)"
  echo ""
  echo "Optional parameters:"
  echo "  -C <complexity_metric>  Complexity metric (default: diagonal-distance)"
  echo "                         Options: edit-distance, diagonal-distance"
  echo "  -M <max_complexity>     Max complexity threshold (default: 100)"
  echo "  -T <type>              Tracepoint type (default: standard)"
  echo "                         Options: standard, mixed, variable, fastga"
  echo "  -t <threads>           Number of threads (default: 8)"
  echo "  -o <output_dir>        Output directory (default: data/benchmarks)"
  echo ""
  echo "Examples:"
  echo "  $0 -p data/alignments/sample.paf -f data/sequences/sample.fa -N 1001"
  echo "  $0 -p data/alignments/sample.paf -f data/sequences/sample.fa -N 1001 -C edit-distance -M 32 -T variable -t 16"
  exit 1
}

while getopts ":p:f:N:C:M:T:t:o:h" opt; do
  case "$opt" in
    p) PAF_FILE="$OPTARG" ;;
    f) FASTA_FILE="$OPTARG" ;;
    N) N_ID="$OPTARG" ;;
    C) CM="$OPTARG" ;;
    M) MC="$OPTARG" ;;
    T) TYPE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    h) usage ;;
    *) echo "Invalid option: -$OPTARG" >&2; usage ;;
  esac
done

[[ -z "$PAF_FILE" || -z "$FASTA_FILE" || -z "$N_ID" ]] && usage

# ---------- Validation ----------
if [[ ! -f "$PAF_FILE" ]]; then
  echo "Error: PAF file not found: $PAF_FILE" >&2
  exit 1
fi

if [[ ! -f "$FASTA_FILE" ]]; then
  echo "Error: FASTA file not found: $FASTA_FILE" >&2
  exit 1
fi

if [[ ! -x "$CIGZIP_BIN" ]]; then
  echo "Error: cigzip binary not found or not executable: $CIGZIP_BIN" >&2
  exit 1
fi

# Validate complexity metric options
case "$CM" in
  edit-distance|diagonal-distance)
    ;;
  *)
    echo "Warning: Unknown complexity metric '$CM'. Proceeding anyway..." >&2
    ;;
esac

# Validate type options
case "$TYPE" in
  standard|mixed|variable|fastga)
    ;;
  *)
    echo "Warning: Unknown tracepoint type '$TYPE'. Proceeding anyway..." >&2
    ;;
esac

# Validate max complexity is a number
if ! [[ "$MC" =~ ^[0-9]+$ ]]; then
  echo "Error: Max complexity must be a positive integer, got: '$MC'" >&2
  exit 1
fi

# Validate threads is a number
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
  echo "Error: Threads must be a positive integer, got: '$THREADS'" >&2
  exit 1
fi

# ---------- Setup ----------
mkdir -p "$OUTPUT_DIR" reports

# Report filename - include all relevant parameters
TS="$(date +%Y%m%d-%H%M%S)"
REPORT="reports/paf_benchmark_${TYPE}_N${N_ID}_${CM}_mc${MC}_${TS}.md"

# Get input file info
PAF_BASENAME=$(basename "$PAF_FILE" .paf)
FASTA_BASENAME=$(basename "$FASTA_FILE")

# Output filenames
ENCODED_PAF="${OUTPUT_DIR}/${PAF_BASENAME}.${TYPE}.${CM}.mc${MC}.encoded.paf"
DECODED_PAF="${OUTPUT_DIR}/${PAF_BASENAME}.${TYPE}.${CM}.mc${MC}.decoded.paf"
MINIMAL_PAF="${OUTPUT_DIR}/${PAF_BASENAME}.${TYPE}.${CM}.mc${MC}.minimal.paf"

# ---------- Helper functions ----------
section()     { echo -e "\n## $1\n" >> "$REPORT"; }
subsection()  { echo -e "\n### $1\n" >> "$REPORT"; }
kv()          { printf "| %s | %s |\n" "$1" "$2" >> "$REPORT"; }
code_bash()   { echo '```bash' >> "$REPORT"; echo "$1" >> "$REPORT"; echo '```' >> "$REPORT"; }
code_text()   { echo '```text' >> "$REPORT"; }
code_end()    { echo '```' >> "$REPORT"; }

filesize() { 
  [[ -f "$1" ]] && stat -c%s "$1" 2>/dev/null || wc -c <"$1" 2>/dev/null || echo "0"
}

filesize_mb() { 
  local size_bytes=$(filesize "$1")
  awk "BEGIN {printf \"%.2f\", $size_bytes/1024/1024}"
}

# Get PAF record count
paf_records() {
  [[ -f "$1" ]] && wc -l < "$1" 2>/dev/null || echo "0"
}

# Capture stdout+stderr of a command into the report
run_and_capture() {
  local desc="$1"; shift
  local cmd="$*"
  subsection "$desc"
  echo "**Command**:" >> "$REPORT"
  code_bash "$cmd"
  echo "**Output**:" >> "$REPORT"
  {
    code_text
    eval "$cmd" 2>&1
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
      echo "Command failed with exit code: $exit_code"
    fi
    code_end
  } >> "$REPORT"
  return ${exit_code:-0}
}

# ---------- Report header ----------
{
  echo "# cigzip PAF Benchmark Report"
  echo
  echo "- Timestamp: \`$(date -Is)\`"
  echo "- Input PAF: \`$(realpath "$PAF_FILE")\`"
  echo "- Input FASTA: \`$(realpath "$FASTA_FILE")\`"
  echo
  echo "## Input Parameters"
  echo
  echo "| Parameter | Value |"
  echo "|---|---|"
} > "$REPORT"

kv "PAF file (-p)" "$PAF_FILE"
kv "FASTA file (-f)" "$FASTA_FILE"
kv "Benchmark ID (-N)" "$N_ID"
kv "Complexity metric (-C)" "$CM"
kv "Max complexity (-M)" "$MC"
kv "Tracepoint type (-T)" "$TYPE"
kv "Threads (-t)" "$THREADS"
kv "Output directory (-o)" "$OUTPUT_DIR"

section "Input File Information"
{
  echo "| File | Size (bytes) | Size (MB) | Records/Lines |"
  echo "|---|---|---|---|"
  echo "| $PAF_FILE | $(filesize "$PAF_FILE") | $(filesize_mb "$PAF_FILE") | $(paf_records "$PAF_FILE") |"
  echo "| $FASTA_FILE | $(filesize "$FASTA_FILE") | $(filesize_mb "$FASTA_FILE") | $(grep -c '^>' "$FASTA_FILE" 2>/dev/null || echo "n/a") |"
} >> "$REPORT"

section "Tool Information"
{
  echo '```text'
  echo "CIGZIP_BIN : $CIGZIP_BIN"
  echo "TIME_BIN   : $TIME_BIN"
  echo
  echo "# cigzip version/help:"
  ($CIGZIP_BIN --version 2>&1 || $CIGZIP_BIN --help 2>&1 | head -n 10 || echo "No version/help output")
  echo '```'
} >> "$REPORT"

section "Benchmark Pipeline"

# ---------- 1) Basic encode (no max-complexity limit) ----------
run_and_capture "1) cigzip encode (baseline, no complexity limit)" \
  "$TIME_BIN $CIGZIP_BIN encode -p \"$PAF_FILE\" --type \"$TYPE\" --complexity-metric \"$CM\" --max-complexity 999999999 -t \"$THREADS\" > \"$ENCODED_PAF\""

# ---------- 2) Basic decode (verification) ----------
run_and_capture "2) cigzip decode (verification)" \
  "$TIME_BIN $CIGZIP_BIN decode -p \"$ENCODED_PAF\" --type \"$TYPE\" --complexity-metric \"$CM\" --sequence-files \"$FASTA_FILE\" -t \"$THREADS\" > \"$DECODED_PAF\""

# ---------- 4) Encode with max-complexity (main benchmark) ----------
run_and_capture "4) cigzip encode with max-complexity=$MC (main benchmark)" \
  "$TIME_BIN $CIGZIP_BIN encode -p \"$DECODED_PAF\" --type \"$TYPE\" --complexity-metric \"$CM\" --minimal --max-complexity \"$MC\" -t \"$THREADS\" > \"$MINIMAL_PAF\""

# ---------- 5) Decode with max-complexity ----------
DECODED_MINIMAL_PAF="${OUTPUT_DIR}/${PAF_BASENAME}.${TYPE}.${CM}.mc${MC}.minimal.decoded.paf"
run_and_capture "5) cigzip decode with max-complexity=$MC" \
  "$TIME_BIN $CIGZIP_BIN decode -p \"$MINIMAL_PAF\" --type \"$TYPE\" --complexity-metric \"$CM\" --sequence-files \"$FASTA_FILE\" --max-complexity \"$MC\" -t \"$THREADS\" > \"$DECODED_MINIMAL_PAF\""

# ---------- 6) Verify max-complexity decode ----------
run_and_capture "6) Verification: diff original vs max-complexity decoded" \
  "diff <(sort \"$DECODED_PAF\") <(sort \"$DECODED_MINIMAL_PAF\") | wc -l"

# ---------- 7) Compression analysis ----------
subsection "7) Compression Analysis"
{
  echo "**Analysis**:" >> "$REPORT"
  code_text
  
  orig_size=$(filesize "$PAF_FILE")
  encoded_size=$(filesize "$ENCODED_PAF")
  minimal_size=$(filesize "$MINIMAL_PAF")
  
  if [[ "$orig_size" -gt 0 ]]; then
    encoded_ratio=$(awk "BEGIN {printf \"%.4f\", $encoded_size/$orig_size}")
    minimal_ratio=$(awk "BEGIN {printf \"%.4f\", $minimal_size/$orig_size}")
    encoded_reduction=$(awk "BEGIN {printf \"%.2f\", (1-$encoded_size/$orig_size)*100}")
    minimal_reduction=$(awk "BEGIN {printf \"%.2f\", (1-$minimal_size/$orig_size)*100}")
    
    echo "Original PAF size: $orig_size bytes ($(filesize_mb "$PAF_FILE") MB)"
    echo "Encoded size (no limit): $encoded_size bytes ($(filesize_mb "$ENCODED_PAF") MB)"
    echo "Encoded size (mc=$MC): $minimal_size bytes ($(filesize_mb "$MINIMAL_PAF") MB)"
    echo
    echo "Compression ratios:"
    echo "- No limit: $encoded_ratio (${encoded_reduction}% reduction)"
    echo "- Max complexity $MC: $minimal_ratio (${minimal_reduction}% reduction)"
    echo
    echo "Space savings with max complexity: $(awk "BEGIN {printf \"%.2f\", ($encoded_size-$minimal_size)/1024/1024}") MB"
  else
    echo "Error: Could not determine original file size"
  fi
  
  code_end
} >> "$REPORT"

# ---------- Performance summary ----------
section "Performance Summary"
{
  echo "| Metric | Value |"
  echo "|---|---|"
  echo "| Original PAF Records | $(paf_records "$PAF_FILE") |"
  echo "| Original Size (MB) | $(filesize_mb "$PAF_FILE") |"
  echo "| Encoded Size (no limit, MB) | $(filesize_mb "$ENCODED_PAF") |"
  echo "| Encoded Size (mc=$MC, MB) | $(filesize_mb "$MINIMAL_PAF") |"
  echo "| Compression Ratio (no limit) | $(awk "BEGIN {orig=$(filesize "$PAF_FILE"); enc=$(filesize "$ENCODED_PAF"); if(orig>0) printf \"%.4f\", enc/orig; else print \"n/a\"}") |"
  echo "| Compression Ratio (mc=$MC) | $(awk "BEGIN {orig=$(filesize "$PAF_FILE"); min=$(filesize "$MINIMAL_PAF"); if(orig>0) printf \"%.4f\", min/orig; else print \"n/a\"}") |"
  echo "| Complexity Metric | $CM |"
  echo "| Max Complexity Threshold | $MC |"
  echo "| Tracepoint Type | $TYPE |"
  echo "| Threads Used | $THREADS |"
} >> "$REPORT"

section "Output Files"
{
  echo "| File | Description | Size (bytes) | Size (MB) |"
  echo "|---|---|---|---|"
  echo "| $ENCODED_PAF | Encoded (no complexity limit) | $(filesize "$ENCODED_PAF") | $(filesize_mb "$ENCODED_PAF") |"
  echo "| $DECODED_PAF | Decoded (verification) | $(filesize "$DECODED_PAF") | $(filesize_mb "$DECODED_PAF") |"
  echo "| $MINIMAL_PAF | Encoded (mc=$MC) | $(filesize "$MINIMAL_PAF") | $(filesize_mb "$MINIMAL_PAF") |"
  echo "| $DECODED_MINIMAL_PAF | Decoded (mc=$MC) | $(filesize "$DECODED_MINIMAL_PAF") | $(filesize_mb "$DECODED_MINIMAL_PAF") |"
} >> "$REPORT"

# ---------- Clean up temporary verification files (optional) ----------
if [[ -f "$DECODED_PAF" && -f "$DECODED_MINIMAL_PAF" ]]; then
  echo ""
  echo "Note: Verification files created. You may want to delete them to save space:"
  echo "  rm \"$DECODED_PAF\" \"$DECODED_MINIMAL_PAF\""
fi

echo ""
echo "✓ Benchmark completed successfully!"
echo "📊 Report written: $REPORT"
echo "📁 Output files in: $OUTPUT_DIR"

# Optional: Show quick summary
echo ""
echo "=== Quick Summary ==="
echo "Original size:    $(filesize_mb "$PAF_FILE") MB"
echo "Encoded (no lim): $(filesize_mb "$ENCODED_PAF") MB"
echo "Encoded (mc=$MC):  $(filesize_mb "$MINIMAL_PAF") MB"
echo "Compression:      $(awk "BEGIN {orig=$(filesize "$PAF_FILE"); min=$(filesize "$MINIMAL_PAF"); if(orig>0) printf \"%.1f%%\", (1-min/orig)*100; else print \"n/a\"}")"
