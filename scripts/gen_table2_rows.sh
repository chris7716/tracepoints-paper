#!/bin/bash
# Generate Table 2 (tab:comprehensive) rows from a real-pangenome benchmark.results.tsv.
# Usage: gen_table2_rows.sh <benchmark.results.tsv> <uncompressed_PAF_GiB>
# Per (cm, mc, distance): storage_GiB, tracepoints, compr_ratio, recon_h, mem_GiB, improved%.
#   storage = sum(size_tpa_bytes)/2^30
#   ratio   = storage_GiB / uncompressed_PAF_GiB
#   recon_h = sum(decompress_runtime_sec + decode_runtime_sec)/3600
#   mem_GiB = max(decode_memory_kb)/2^20                (peak across files)
#   improved% = 100*sum(score_improved)/sum(num_output_alignments)
set -euo pipefail
awk -F'\t' -v U="$2" '
NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
{
  k = $(h["cm"]) "|" $(h["mc"]) "|" $(h["distance"])
  if ($(h["decode_runtime_sec"]) == "NA") next
  tpa[k]  += $(h["size_tpa_bytes"]);  tp[k] += $(h["num_tracepoints"])
  rt[k]   += $(h["decompress_runtime_sec"]) + $(h["decode_runtime_sec"])
  m = $(h["decode_memory_kb"]) + 0;   if (m > mx[k]) mx[k] = m
  imp[k]  += $(h["score_improved"]);  nout[k] += $(h["num_output_alignments"]);  n[k]++
}
END {
  printf "%-34s %6s %8s %10s %8s %8s %8s %8s\n","cm|mc|dist","files","storGiB","tracepts","ratio","recon_h","memGiB","impr%"
  for (k in n) { g = tpa[k]/1073741824
    printf "%-34s %6d %8.2f %9.1fM %8.3f %8.2f %8.2f %8.2f\n",
           k, n[k], g, tp[k]/1e6, g/U, rt[k]/3600, mx[k]/1048576, (nout[k]?100*imp[k]/nout[k]:0)
  }
}' "$1" | sort
