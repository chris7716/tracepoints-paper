#!/usr/bin/env bash
# Build one self-contained FASTGA GDB from a set of (gzipped) FASTAs, decoupled from the
# benchmark so it can run on a big-scratch machine (e.g. tux09) and the GDB shipped to the
# experiment node. The GDB (<dataset>.1gdb + hidden .<dataset>.bps) holds a 2-bit copy of
# every base; PAFtoALN/ALNtoPAF read from it, so the FASTAs are NOT needed on the exp node.
#
# Usage:
#   [FASTGA=/scratch/FASTGA] [DEST=user@host:/path] \
#     ./scripts/build-fastga-gdb.sh <dataset> <out_dir> <fasta.gz...>
#
# Examples (run on tux09):
#   ./scripts/build-fastga-gdb.sh hprcv2-25k /scratch/gdb /moosefs/pangenomes/HPRCv2/*.fa.gz
#   ./scripts/build-fastga-gdb.sh primates  /scratch/gdb /moosefs/guarracino/t2t-ape-pangenome/*.fa.gz
#   # stream the finished GDB straight to the experiment node (no intermediate tar on scratch):
#   DEST=guarracino@tuxNN:/scratch ./scripts/build-fastga-gdb.sh hprcv2-25k /scratch/gdb /moosefs/pangenomes/HPRCv2/*.fa.gz
#
# NOTE: the GDB scaffold names come from the FASTA headers, so these FASTAs must use the
# same sequence names as the PAF query/target columns (PanSN), or PAFtoALN won't resolve.

set -euo pipefail

FASTGA="${FASTGA:-/scratch/FASTGA}"
FAtoGDB="$FASTGA/FAtoGDB"
DEST="${DEST:-}"   # optional user@host:/dir — if set, stream the GDB there with no local tar

[ $# -ge 3 ] || { echo "usage: $0 <dataset> <out_dir> <fasta.gz...>" >&2; exit 1; }
dataset="$1"; out_dir="$2"; shift 2
fastas=("$@")
[ -x "$FAtoGDB" ] || { echo "ERROR: FAtoGDB not executable at $FAtoGDB (set FASTGA=...)" >&2; exit 1; }
mkdir -p "$out_dir"

gdb="$out_dir/$dataset.1gdb"
bps="$out_dir/.$dataset.bps"
combined="$out_dir/$dataset.combined.fa.gz"

echo "[build-gdb] dataset=$dataset  inputs=${#fastas[@]}  out_dir=$out_dir  FASTGA=$FASTGA"

# Inputs are already gzipped; concatenating gzip members gives one valid multi-member gzip
# stream that FAtoGDB (zlib gzopen) reads transparently. No decompress/recompress.
echo "[build-gdb] concatenating ${#fastas[@]} gzipped FASTAs -> $combined"
cat "${fastas[@]}" > "$combined"

echo "[build-gdb] FAtoGDB -> $gdb (this is the heavy step)"
"$FAtoGDB" "$combined" "$gdb" > "$out_dir/$dataset.gdb.log" 2>&1
rm -f "$combined"   # GDB is self-contained; drop the temporary combined FASTA

[ -f "$gdb" ] || { echo "ERROR: $gdb not produced (see $out_dir/$dataset.gdb.log)" >&2; exit 1; }
[ -f "$bps" ] || { echo "ERROR: expected bases file $bps not found; check FAtoGDB output naming" >&2; exit 1; }

# GDB component files to ship (always .1gdb + hidden .bps; .1ano only if masking was detected)
comps=("$dataset.1gdb" ".$dataset.bps")
[ -f "$out_dir/$dataset.1ano" ]  && comps+=("$dataset.1ano")
[ -f "$out_dir/.$dataset.1ano" ] && comps+=(".$dataset.1ano")

echo "[build-gdb] GDB ready:"; ls -lah "$gdb" "$bps" 2>/dev/null

if [ -n "$DEST" ]; then
    # Stream the GDB to the experiment node without writing an intermediate tar to scratch.
    echo "[build-gdb] streaming GDB -> $DEST"
    ( cd "$out_dir" && tar -cf - "${comps[@]}" ) | ssh "${DEST%%:*}" "mkdir -p '${DEST#*:}' && tar -C '${DEST#*:}' -xf -"
    echo "[build-gdb] transferred. On the experiment node the GDB is ${DEST#*:}/$dataset.1gdb"
else
    echo
    echo "[build-gdb] To transfer WITHOUT an intermediate tar (avoids doubling scratch):"
    echo "  ( cd $out_dir && tar -cf - ${comps[*]} ) | ssh user@expnode \"mkdir -p /dest && tar -C /dest -xf -\""
    echo "Then run_benchmark_real's fastga-native branch must point at /dest/$dataset.1gdb"
    echo "(it reads \$scratch_dir/<dataset>.1gdb, so place it there or symlink)."
fi

echo "[build-gdb] done."
