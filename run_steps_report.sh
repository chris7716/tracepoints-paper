#!/usr/bin/env bash
set -euo pipefail

# Runs exactly your listed steps and writes a Markdown report that
# includes:
#  - Input params (without FastGA tmp)
#  - Each command (verbatim)
#  - **Output** of each sub-command (stdout + stderr), including /usr/bin/time -v
#
# Usage:
#   ./run_steps_report.sh -r <record_count> -L <length> -E <error> -N <N> -T <type>

# ---------- Defaults / paths ----------
FASTGA_TMP="/home/hasitha/data/projects/fastga-tmp"
CM="edit-distance"
MC="32"
THREADS="8"

# Tool paths (override via env if needed)
WFA2_BIN="${WFA2_BIN:-../WFA2-lib/bin/generate_dataset}"
FASTGA_BIN="${FASTGA_BIN:-../FASTGA/FastGA}"
CIGZIP_BIN="${CIGZIP_BIN:-../cigzip/target/debug/cigzip}"

# Prefer /usr/bin/time -v if available
TIME_BIN="$(command -v /usr/bin/time >/dev/null 2>&1 && echo /usr/bin/time || echo time)"

# ---------- Args ----------
RECORDS="" LENGTH="" ERROR="" N_ID="" TYPE=""
usage() {
  echo "Usage: $0 -r <record_count> -L <length> -E <error> -N <N> -T <type>"
  exit 1
}

while getopts ":r:L:E:N:T:" opt; do
  case "$opt" in
    r) RECORDS="$OPTARG" ;;
    L) LENGTH="$OPTARG" ;;
    E) ERROR="$OPTARG" ;;
    N) N_ID="$OPTARG" ;;
    T) TYPE="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "$RECORDS" || -z "$LENGTH" || -z "$ERROR" || -z "$N_ID" || -z "$TYPE" ]] && usage

# ---------- Checks & prep ----------
for exe in "$WFA2_BIN" "$FASTGA_BIN" "$CIGZIP_BIN"; do
  [[ -x "$exe" ]] || { echo "Missing executable: $exe" >&2; exit 2; }
done

export PATH="$PATH:$(pwd)"
mkdir -p wf-dataset fastga-datasets cigzip-datasets tmp reports

# Report filename
TS="$(date +%Y%m%d-%H%M%S)"
REPORT="reports/run_${TYPE}_N${N_ID}_${TS}.md"

# ---------- Filenames (exactly as in your steps) ----------
SAMPLE_SEQ="data/simulated/wf-dataset/sample.dataset.${N_ID}.seq"
FASTA="data/simulated/wf-dataset/dataset.${N_ID}.fa"
FASTGA_PAF="data/simulated/fastga-datasets/dataset.${N_ID}.paf"

TMP_TYPED="data/simulated/tmp/tmp.${TYPE}.${N_ID}.paf"
CIGZIP_OUT_TYPED="data/simulated/cigzip-datasets/dataset.${TYPE}.${N_ID}.paf"

name="dataset.${TYPE}.${N_ID}"
cm="${CM}"
mc="${MC}"

# ---------- Helpers ----------
section()     { echo -e "\n## $1\n" >> "$REPORT"; }
subsection()  { echo -e "\n### $1\n" >> "$REPORT"; }
kv()          { printf "| %s | %s |\n" "$1" "$2" >> "$REPORT"; }
code_bash()   { echo '```bash' >> "$REPORT"; echo "$1" >> "$REPORT"; echo '```' >> "$REPORT"; }
code_text()   { echo '```text' >> "$REPORT"; }
code_end()    { echo '```' >> "$REPORT"; }
filesize()    { [[ -f "$1" ]] && stat -c%s "$1" 2>/dev/null || wc -c <"$1" 2>/dev/null || echo "n/a"; }
filesize_mb() { 
  local size_bytes=$(stat -c%s "$1" 2>/dev/null || wc -c <"$1" 2>/dev/null || echo "0")
  awk "BEGIN {printf \"%.2f\", $size_bytes/1024/1024}"
}

# Capture stdout+stderr of a command into the report under **Output**
run_and_capture() {
  local desc="$1"; shift
  local cmd="$*"
  subsection "$desc"
  echo "**Command**:" >> "$REPORT"
  code_bash "$cmd"
  echo "**Output**:" >> "$REPORT"
  {
    code_text
    bash -c "$cmd" 2>&1
    code_end
  } >> "$REPORT"
}

# ---------- Report header ----------
{
  echo "# cigzip / FastGA run"
  echo
  echo "- Timestamp: \`$(date -Is)\`"
  echo
  echo "## Inputs"
  echo
  echo "| Parameter | Value |"
  echo "|---|---|"
} > "$REPORT"

kv "record_count (-r)" "$RECORDS"
kv "length (-L)" "$LENGTH"
kv "error (-E)" "$ERROR"
kv "N (-N)" "$N_ID"
kv "type (-T)" "$TYPE"
kv "complexity metric (cm)" "$cm"
kv "max complexity (mc)" "$mc"
kv "threads" "$THREADS"

section "Tool Binaries"
{
  echo '```text'
  echo "WFA2_BIN   : $WFA2_BIN"
  echo "FASTGA_BIN : $FASTGA_BIN"
  echo "CIGZIP_BIN : $CIGZIP_BIN"
  echo
  echo "# WFA2 help/version (if any):"
  ($WFA2_BIN --help 2>&1 | head -n 5) || echo "(no output)"
  echo
  echo "# FastGA banner (no -h):"
  ($FASTGA_BIN 2>&1 | head -n 5) || echo "(no output)"
  echo
  echo "# cigzip help (if any):"
  ($CIGZIP_BIN --help 2>&1 | head -n 5) || echo "(no output)"
  echo '```'
} >> "$REPORT"

section "Pipeline"

# ---------- 1) Generate dataset ----------
run_and_capture "1) WFA2 generate_dataset" \
  "$WFA2_BIN -n ${RECORDS} -l ${LENGTH} -e ${ERROR} -o ${SAMPLE_SEQ}"

# ---------- 2) Convert to FASTA (exact awk) ----------
subsection "2) Convert to FASTA (awk)"
echo "**Command**:" >> "$REPORT"
code_bash "awk '
BEGIN { counter=1 }
/^>/ { 
    print \">pattern_\" counter
    print substr(\$0, 2)
}
/^</ { 
    print \">text_\" counter
    print substr(\$0, 2)
    counter++
}' ${SAMPLE_SEQ} > ${FASTA}"
echo "**Output**:" >> "$REPORT"
{
  code_text
  awk '
  BEGIN { counter=1 }
  /^>/ {
      print ">pattern_" counter
      print substr($0, 2)
  }
  /^</ {
      print ">text_" counter
      print substr($0, 2)
      counter++
  }' "${SAMPLE_SEQ}" > "${FASTA}" 2>&1
  echo "FASTA written: ${FASTA}"
  code_end
} >> "$REPORT"

# ---------- 3) FastGA self-alignment ----------
run_and_capture "3) FastGA self-alignment → ${FASTGA_PAF}" \
  "$FASTGA_BIN -P${FASTGA_TMP} -pafx ${FASTA} ${FASTA} > ${FASTGA_PAF}"

# ---------- 4) cigzip encode (type=<type>) ----------
run_and_capture "4) cigzip encode (type=${TYPE}) → ${TMP_TYPED}" \
  "$CIGZIP_BIN encode -p ${FASTGA_PAF} --type standard --complexity-metric edit-distance --max-complexity 99999999 > ${TMP_TYPED}"

# ---------- 5) cigzip decode (type=<type>) ----------
run_and_capture "5) cigzip decode (type=${TYPE}) → ${CIGZIP_OUT_TYPED}" \
  "$CIGZIP_BIN decode -p ${TMP_TYPED} --type standard --complexity-metric edit-distance --sequence-files ${FASTA} > ${CIGZIP_OUT_TYPED}"

# ---------- 6) cigzip encode with max-complexity (tracepoints) ----------
run_and_capture "6) cigzip encode (max-complexity=${mc}, threads=${THREADS}) --minimal → cigzip-datasets/${name}.tp.mc${mc}.paf" \
  "${TIME_BIN} ${CIGZIP_BIN} encode -p data/simulated/cigzip-datasets/${name}.paf --type ${TYPE} --complexity-metric ${cm} --max-complexity ${mc} -t ${THREADS} --minimal > data/simulated/cigzip-datasets/${name}.tp.mc${mc}.paf"

# ---------- 7) cigzip decode + diff verification ----------
#subsection "7) cigzip decode + diff verification"
#echo "**Command**:" >> "$REPORT"
#code_bash "${TIME_BIN} -v ${CIGZIP_BIN} decode -p data/simulated/cigzip-datasets/${name}.tp.mc${mc}.paf --type ${TYPE} --complexity-metric ${cm} --sequence-files ${FASTA} --max-complexity ${mc} -t ${THREADS} > data/simulated/cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf && diff <(sort data/simulated/cigzip-datasets/${name}.paf) <(sort data/simulated/cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf) | wc -l"
#echo "**Output**:" >> "$REPORT"
#{
#  code_text
#  ${TIME_BIN} -v "${CIGZIP_BIN}" decode -p "data/simulated/cigzip-datasets/${name}.tp.mc${mc}.paf" --type "${TYPE}" --complexity-metric "${cm}" --sequence-files "${FASTA}" --max-complexity "${mc}" -t "${THREADS}" > "data/simulated/cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf"
#  diff <(sort "data/simulated/cigzip-datasets/${name}.paf") <(sort "data/simulated/cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf") | wc -l
#  code_end
#} >> "$REPORT"

# ---------- 7) cigzip decode + diff verification ----------
run_and_capture "7) cigzip decode + diff verification" \
  "${TIME_BIN} -v ${CIGZIP_BIN} decode -p data/simulated/cigzip-datasets/${name}.tp.mc${mc}.paf --type ${TYPE} --complexity-metric ${cm} --sequence-files ${FASTA} --max-complexity ${mc} -t ${THREADS} > data/simulated/cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf && diff <(sort data/simulated/cigzip-datasets/${name}.paf) <(sort data/simulated/cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf) | wc -l"

# ---------- Artifacts summary ----------
section "Artifacts"
{
  echo
  echo "| File | Size (bytes) |"
  echo "|---|---|"
  echo "| ${SAMPLE_SEQ} | $(filesize "${SAMPLE_SEQ}") |"
  echo "| ${FASTA} | $(filesize "${FASTA}") |"
  echo "| ${FASTGA_PAF} | $(filesize "${FASTGA_PAF}") |"
  echo "| ${TMP_TYPED} | $(filesize "${TMP_TYPED}") |"
  echo "| ${CIGZIP_OUT_TYPED} | $(filesize "${CIGZIP_OUT_TYPED}") |"
  echo "| cigzip-datasets/${name}.tp.mc${mc}.paf | $(filesize "data/simulated/cigzip-datasets/${name}.tp.mc${mc}.paf") |"
  echo "| cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf | $(filesize "data/simulated/cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf") |"
} >> "$REPORT"

echo "Report written: ${REPORT}"

section "Artifacts MB"
{
  echo
  echo "| ${SAMPLE_SEQ} | ${FASTA} | ${FASTGA_PAF} | ${TMP_TYPED} | ${CIGZIP_OUT_TYPED} | cigzip-datasets/${name}.tp.mc${mc}.paf | cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf |"
  echo "|---|---|---|---|---|---|---|"
  echo "| $(filesize_mb "${SAMPLE_SEQ}") | $(filesize_mb "${FASTA}") | $(filesize_mb "${FASTGA_PAF}") | $(filesize_mb "${TMP_TYPED}") | $(filesize_mb "${CIGZIP_OUT_TYPED}") | $(filesize_mb "data/simulated/cigzip-datasets/${name}.tp.mc${mc}.paf") | $(filesize_mb "data/simulated/cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf") |"
} >> "$REPORT"

echo "Report written: ${REPORT}"
