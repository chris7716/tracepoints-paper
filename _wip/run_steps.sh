#!/usr/bin/env bash
set -euo pipefail

# Runs:
# 1) WFA2 generate_dataset
# 2) awk -> FASTA conversion
# 3) FastGA self-alignment to PAF
# 4) cigzip encode (type=<type>)
# 5) cigzip decode (type=<type>)
# 6) cigzip encode with max-complexity (tracepoints)
# 7) cigzip decode + diff verification
#
# Matches your filenames and flags exactly.

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
mkdir -p wf-dataset fastga-datasets cigzip-datasets tmp

# ---------- Filenames (exactly as in your steps) ----------
SAMPLE_SEQ="wf-dataset/sample.dataset.${N_ID}.seq"
FASTA="wf-dataset/dataset.${N_ID}.fa"
FASTGA_PAF="fastga-datasets/dataset.${N_ID}.paf"

TMP_TYPED="tmp/tmp.${TYPE}.${N_ID}.paf"
CIGZIP_OUT_TYPED="cigzip-datasets/dataset.${TYPE}.${N_ID}.paf"

name="dataset.${TYPE}.${N_ID}"
cm="${CM}"
mc="${MC}"

# ---------- 1) Generate dataset ----------
"$WFA2_BIN" -n "${RECORDS}" -l "${LENGTH}" -e "${ERROR}" -o "${SAMPLE_SEQ}"

# ---------- 2) Convert to FASTA (exact awk) ----------
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
}' "${SAMPLE_SEQ}" > "${FASTA}"

# ---------- 3) FastGA self-alignment ----------
"$FASTGA_BIN" -P"${FASTGA_TMP}" -pafx "${FASTA}" "${FASTA}" > "${FASTGA_PAF}"

# ---------- 4) cigzip encode (type=<type>) ----------
"$CIGZIP_BIN" encode -p "${FASTGA_PAF}" --type "${TYPE}" --complexity-metric edit-distance --max-complexity 99999999 > "${TMP_TYPED}"

# ---------- 5) cigzip decode (type=<type>) ----------
"$CIGZIP_BIN" decode -p "${TMP_TYPED}" --type "${TYPE}" --complexity-metric edit-distance --sequence-files "${FASTA}" > "${CIGZIP_OUT_TYPED}"

# ---------- 6) cigzip encode with max-complexity (tracepoints) ----------
"$CIGZIP_BIN" encode -p "cigzip-datasets/${name}.paf" --type "${TYPE}" --complexity-metric "${cm}" --max-complexity "${mc}" -t "${THREADS}" > "cigzip-datasets/${name}.tp.mc${mc}.paf"

# ---------- 7) cigzip decode + diff verification ----------
${TIME_BIN} -v "$CIGZIP_BIN" decode -p "cigzip-datasets/${name}.tp.mc${mc}.paf" --type "${TYPE}" --complexity-metric "${cm}" --sequence-files "${FASTA}" --max-complexity "${mc}" -t "${THREADS}" > "cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf" && \
diff <(sort "cigzip-datasets/${name}.paf") <(sort "cigzip-datasets/${name}.tp.mc${mc}.decompressed.paf") | wc -l
