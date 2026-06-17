#!/usr/bin/env Rscript

# Compute numbers reported in the ES-2bit discussion paragraph and
# supplementary section. Run after es2bit_bench produces all.tsv.
#
# Usage:
#   Rscript scripts/compute_es2bit_stats.R <all.tsv>

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: compute_es2bit_stats.R <all.tsv>")

raw <- read_tsv(args[1], show_col_types = FALSE)
raw <- raw %>%
  mutate(
    n   = as.integer(str_extract(basename(file), "(?<=set_)\\d+")),
    eps_code = str_extract(basename(file), "(?<=_)\\d+(?=\\.paf)"),
    eps = as.numeric(paste0("0.", str_sub(eps_code, 2)))
  )

cat("eps values in data:", sort(unique(raw$eps)), "\n\n")

# Highest error rate and longest sequences
max_eps <- max(raw$eps)
max_n   <- max(raw$n)
r <- raw %>% filter(n == max_n, eps == max_eps)

cat("=== Discussion paragraph numbers ===\n\n")

cat(sprintf("Condition: n=%s, eps=%s (mean_e=%.1f)\n",
    scales::comma(max_n), max_eps, r$mean_e))
cat(sprintf("DB-TP disk bits/edit: %.2f\n", r$dbtp_disk_bits_per_e))
cat(sprintf("ES-2bit disk bits/edit: %.2f\n", r$es2bit_disk_bits_per_e))
cat(sprintf("Fold difference (ES-2bit / DB-TP): %.0fx\n",
    r$es2bit_disk_bits_per_e / r$dbtp_disk_bits_per_e))
cat(sprintf("ES-2bit decode ns/edit (all conditions): %.0f -- %.0f\n",
    min(raw$es2bit_dec_ns_per_e), max(raw$es2bit_dec_ns_per_e)))

# Decode speed ratios vs ES-2bit
cat(sprintf("FL-TP / ES-2bit decode ratio: %.1fx -- %.1fx\n",
    min(raw$fltp_dec_ns_per_e / raw$es2bit_dec_ns_per_e),
    max(raw$fltp_dec_ns_per_e / raw$es2bit_dec_ns_per_e)))
cat(sprintf("EB-TP / ES-2bit decode ratio: %.1fx -- %.1fx\n",
    min(raw$ebtp_dec_ns_per_e / raw$es2bit_dec_ns_per_e),
    max(raw$ebtp_dec_ns_per_e / raw$es2bit_dec_ns_per_e)))
cat(sprintf("DB-TP / ES-2bit decode ratio: %.1fx -- %.1fx\n",
    min(raw$dbtp_dec_ns_per_e / raw$es2bit_dec_ns_per_e),
    max(raw$dbtp_dec_ns_per_e / raw$es2bit_dec_ns_per_e)))

# EB-TP bits/edit at n=100k
r100k <- raw %>% filter(n == max_n)
cat(sprintf("EB-TP disk bits/edit at n=%s: %.2f -- %.2f\n",
    scales::comma(max_n), min(r100k$ebtp_disk_bits_per_e), max(r100k$ebtp_disk_bits_per_e)))

cat("\n=== Supplementary section numbers ===\n\n")

r100k <- raw %>% filter(n == max_n)
cat(sprintf("DB-TP decode ns/edit at n=%s: %.0f -- %.0f\n",
    scales::comma(max_n),
    min(r100k$dbtp_dec_ns_per_e), max(r100k$dbtp_dec_ns_per_e)))
cat(sprintf("CIGAR bits/edit (all conditions): %.0f -- %.0f\n",
    min(raw$cigar_disk_bits_per_e), max(raw$cigar_disk_bits_per_e)))

# Gzip gap
cat(sprintf("After gzip at n=%s, eps=%s:\n", scales::comma(max_n), max_eps))
cat(sprintf("  DB-TP gz B/edit: %.4f\n", r$dbtp_gz_B_per_e))
cat(sprintf("  ES-2bit gz B/edit: %.4f\n", r$es2bit_gz_B_per_e))
gz_ratio <- r$es2bit_gz_B_per_e / r$dbtp_gz_B_per_e
cat(sprintf("  Gzip ratio (ES-2bit / DB-TP): %.0fx\n", gz_ratio))

cat("\nDone.\n")
