#!/usr/bin/env Rscript

# Script to compute all numbers reported in Figure 1 caption and text
# for the tracepoints paper

library(tidyverse)

# Configuration: set data directory via environment variable or use default relative path
data_dir <- Sys.getenv("TRACEPOINTS_DATA_DIR", unset = "data")
sim_data_dir <- file.path(data_dir, "simulated-data")
real_data_dir <- file.path(data_dir, "real-data")

# Read data
data <- read_tsv(file.path(sim_data_dir, "benchmark.results.tsv"), show_col_types = FALSE)

# Filter for b=32 / delta=32 (mc=32) and high memory mode
df <- data %>%
  filter(memory_mode == "high") %>%
  filter((tp_type == "fastga" & mc == 100) | (tp_type == "standard" & mc == 32))

# Calculate compression ratios (smaller = better compression)
df <- df %>%
  mutate(
    ratio_bgzip = size_cigar_bgzip_bytes / size_cigar_bytes,
    ratio_tpa = size_tpa_bytes / size_cigar_bytes,
    method = case_when(
      tp_type == "fastga" ~ "FL-TP",
      cm == "edit-distance" ~ "EB-TP",
      cm == "diagonal-distance" ~ "DB-TP"
    )
  )

cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat("COMPRESSION STATISTICS FOR FIGURE 1\n")
cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# =============================================================================
# Claim 1: "For short sequences (100 bp), BGZIP achieves the best compression
#           (0.14--0.23 ratio across error rates)"
# =============================================================================

cat("CLAIM 1: BGZIP compression ratios at 100 bp across error rates\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

bgzip_100bp <- df %>%
  filter(l == 100, method == "FL-TP") %>%  # Use FL-TP rows which have BGZIP data
  select(e, ratio_bgzip) %>%
  arrange(e)

print(bgzip_100bp)

cat(sprintf("\nBGZIP at 100 bp: range %.2f -- %.2f\n",
            min(bgzip_100bp$ratio_bgzip),
            max(bgzip_100bp$ratio_bgzip)))

# Compare with other methods at 100 bp
cat("\nAll methods at 100 bp:\n")
methods_100bp <- df %>%
  filter(l == 100) %>%
  select(e, method, ratio_tpa) %>%
  pivot_wider(names_from = method, values_from = ratio_tpa) %>%
  left_join(bgzip_100bp, by = "e") %>%
  rename(BGZIP = ratio_bgzip)

print(methods_100bp)

cat("\n")

# =============================================================================
# Claim 2: "At 100 Kb, DB-TP TPA achieves 10--13× better compression than FL-TP
#           and 30--146× better than BGZIP across error rates"
# =============================================================================

cat("CLAIM 2: DB-TP improvement over FL-TP and BGZIP at 100 Kb\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

data_100kb <- df %>%
  filter(l == 100000) %>%
  select(e, method, ratio_tpa, ratio_bgzip)

# Get DB-TP ratios
dbtp_100kb <- data_100kb %>%
  filter(method == "DB-TP") %>%
  select(e, dbtp_ratio = ratio_tpa)

# Get FL-TP ratios
fltp_100kb <- data_100kb %>%
  filter(method == "FL-TP") %>%
  select(e, fltp_ratio = ratio_tpa, bgzip_ratio = ratio_bgzip)

# Compute improvement factors
comparison_100kb <- dbtp_100kb %>%
  left_join(fltp_100kb, by = "e") %>%
  mutate(
    dbtp_vs_fltp = fltp_ratio / dbtp_ratio,
    dbtp_vs_bgzip = bgzip_ratio / dbtp_ratio
  )

print(comparison_100kb)

cat(sprintf("\nDB-TP vs FL-TP at 100 Kb: %.1f -- %.1f× improvement\n",
            min(comparison_100kb$dbtp_vs_fltp),
            max(comparison_100kb$dbtp_vs_fltp)))

cat(sprintf("DB-TP vs BGZIP at 100 Kb: %.0f -- %.0f× improvement\n",
            min(comparison_100kb$dbtp_vs_bgzip),
            max(comparison_100kb$dbtp_vs_bgzip)))

cat("\n")

# =============================================================================
# Claim 3: "EB-TP TPA achieves 1.3--5.2× improvement over FL-TP at 100 Kb"
# =============================================================================

cat("CLAIM 3: EB-TP improvement over FL-TP at 100 Kb\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

# Get EB-TP ratios
ebtp_100kb <- data_100kb %>%
  filter(method == "EB-TP") %>%
  select(e, ebtp_ratio = ratio_tpa)

# Compute improvement factors
comparison_ebtp_100kb <- ebtp_100kb %>%
  left_join(fltp_100kb, by = "e") %>%
  mutate(
    ebtp_vs_fltp = fltp_ratio / ebtp_ratio,
    ebtp_vs_bgzip = bgzip_ratio / ebtp_ratio
  )

print(comparison_ebtp_100kb)

cat(sprintf("\nEB-TP vs FL-TP at 100 Kb: %.1f -- %.1f× improvement\n",
            min(comparison_ebtp_100kb$ebtp_vs_fltp),
            max(comparison_ebtp_100kb$ebtp_vs_fltp)))

cat("\n")

# =============================================================================
# Additional: Full table of compression ratios for reference
# =============================================================================

cat("FULL COMPRESSION RATIO TABLE (all lengths, all error rates)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

full_table <- df %>%
  select(l, e, method, ratio_tpa) %>%
  pivot_wider(names_from = method, values_from = ratio_tpa) %>%
  left_join(
    df %>%
      filter(method == "FL-TP") %>%
      select(l, e, BGZIP = ratio_bgzip),
    by = c("l", "e")
  ) %>%
  arrange(l, e) %>%
  select(l, e, BGZIP, `FL-TP`, `EB-TP`, `DB-TP`)

print(full_table, n = Inf)

cat("\n")

# =============================================================================
# Summary statistics
# =============================================================================

cat("SUMMARY: Compression improvement by sequence length\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

summary_by_length <- df %>%
  select(l, e, method, ratio_tpa, ratio_bgzip) %>%
  pivot_wider(names_from = method, values_from = ratio_tpa) %>%
  mutate(
    dbtp_vs_fltp = `FL-TP` / `DB-TP`,
    dbtp_vs_bgzip = ratio_bgzip / `DB-TP`,
    ebtp_vs_fltp = `FL-TP` / `EB-TP`
  ) %>%
  group_by(l) %>%
  summarise(
    dbtp_vs_fltp_min = min(dbtp_vs_fltp),
    dbtp_vs_fltp_max = max(dbtp_vs_fltp),
    dbtp_vs_bgzip_min = min(dbtp_vs_bgzip),
    dbtp_vs_bgzip_max = max(dbtp_vs_bgzip),
    ebtp_vs_fltp_min = min(ebtp_vs_fltp),
    ebtp_vs_fltp_max = max(ebtp_vs_fltp)
  )

print(summary_by_length)

cat("\n")

# =============================================================================
# DECODING COST STATISTICS (TPA -> PAF)
# =============================================================================

cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat("DECODING COST STATISTICS (TPA -> PAF)\n")
cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Compute TPA -> PAF time = decompress + decode
df_decode <- df %>%
  mutate(
    # TPA -> PAF runtime
    tpa_to_paf_runtime = case_when(
      method == "FL-TP" ~ decompress_runtime_sec + decode_runtime_sec,
      TRUE ~ decompress_runtime_sec + decode_heuristic_runtime_sec
    ),
    # TPA -> PAF peak memory (max of decompress and decode)
    tpa_to_paf_memory_kb = case_when(
      method == "FL-TP" ~ pmax(decompress_memory_kb, decode_memory_kb),
      TRUE ~ pmax(decompress_memory_kb, decode_heuristic_memory_kb)
    ),
    tpa_to_paf_memory_mb = tpa_to_paf_memory_kb / 1024
  )

cat("Decoding costs at 100 Kb, 10% error rate:\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

decode_100kb_10pct <- df_decode %>%
  filter(l == 100000, e == 0.1) %>%
  select(method, align_runtime_sec, tpa_to_paf_runtime, tpa_to_paf_memory_mb)

print(decode_100kb_10pct)

# Get alignment time (ORIGINAL)
align_time <- decode_100kb_10pct %>% filter(method == "FL-TP") %>% pull(align_runtime_sec)

cat(sprintf("\nORIGINAL alignment time: %.2f sec\n", align_time))
cat(sprintf("FL-TP TPA->PAF time: %.2f sec (%.1f× faster)\n",
            decode_100kb_10pct %>% filter(method == "FL-TP") %>% pull(tpa_to_paf_runtime),
            align_time / (decode_100kb_10pct %>% filter(method == "FL-TP") %>% pull(tpa_to_paf_runtime))))
cat(sprintf("EB-TP TPA->PAF time: %.2f sec (%.1f× faster)\n",
            decode_100kb_10pct %>% filter(method == "EB-TP") %>% pull(tpa_to_paf_runtime),
            align_time / (decode_100kb_10pct %>% filter(method == "EB-TP") %>% pull(tpa_to_paf_runtime))))
cat(sprintf("DB-TP TPA->PAF time: %.2f sec (%.1f× faster)\n",
            decode_100kb_10pct %>% filter(method == "DB-TP") %>% pull(tpa_to_paf_runtime),
            align_time / (decode_100kb_10pct %>% filter(method == "DB-TP") %>% pull(tpa_to_paf_runtime))))

cat("\n")
cat("Full decoding costs table (100 Kb, all error rates):\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

decode_100kb_all <- df_decode %>%
  filter(l == 100000) %>%
  select(e, method, align_runtime_sec, tpa_to_paf_runtime, tpa_to_paf_memory_mb) %>%
  arrange(e, method)

print(decode_100kb_all)

cat("\n")
cat("Speedup over re-alignment (100 Kb):\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

speedup_table <- df_decode %>%
  filter(l == 100000) %>%
  select(e, method, align_runtime_sec, tpa_to_paf_runtime) %>%
  mutate(speedup = align_runtime_sec / tpa_to_paf_runtime) %>%
  select(e, method, speedup) %>%
  pivot_wider(names_from = method, values_from = speedup)

print(speedup_table)

cat(sprintf("\nSpeedup range: %.0f--%.0f× faster than re-alignment\n",
            min(speedup_table$`FL-TP`, speedup_table$`EB-TP`, speedup_table$`DB-TP`),
            max(speedup_table$`FL-TP`, speedup_table$`EB-TP`, speedup_table$`DB-TP`)))

cat("\n")

# =============================================================================
# DECOMPRESSION PERCENTAGE OF TOTAL RECONSTRUCTION TIME
# =============================================================================

cat("DECOMPRESSION AS PERCENTAGE OF TOTAL RECONSTRUCTION TIME\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

df_decomp_frac <- data %>%
  filter(memory_mode == "high") %>%
  filter((tp_type == "fastga" & mc == 100) | (tp_type == "standard" & mc == 32)) %>%
  filter(e != 0.05) %>%
  mutate(
    method = case_when(
      tp_type == "fastga" ~ "FL-TP",
      cm == "edit-distance" ~ "EB-TP",
      cm == "diagonal-distance" ~ "DB-TP"
    ),
    total_runtime = case_when(
      method == "FL-TP" ~ decompress_runtime_sec + decode_runtime_sec,
      TRUE ~ decompress_runtime_sec + decode_heuristic_runtime_sec
    ),
    decomp_pct = ifelse(total_runtime > 0,
                        decompress_runtime_sec / total_runtime * 100, 0),
    length_label = factor(
      case_when(
        l == 100 ~ "100 bp", l == 1000 ~ "1 Kb",
        l == 10000 ~ "10 Kb", l == 100000 ~ "100 Kb"
      ),
      levels = c("100 bp", "1 Kb", "10 Kb", "100 Kb")
    ),
    error_label = paste0(e * 100, "%")
  ) %>%
  select(length_label, error_label, method, decompress_runtime_sec, total_runtime, decomp_pct) %>%
  arrange(length_label, error_label, method)

cat("Full table:\n")
print(df_decomp_frac, n = Inf)

cat("\nPivoted (decompression %):\n")
decomp_pivot <- df_decomp_frac %>%
  mutate(decomp_pct = round(decomp_pct, 2)) %>%
  select(length_label, error_label, method, decomp_pct) %>%
  pivot_wider(names_from = method, values_from = decomp_pct)
print(decomp_pivot, n = Inf)

cat(sprintf("\nMax decompression fraction: %.1f%% (%s)\n",
            max(df_decomp_frac$decomp_pct),
            df_decomp_frac %>% filter(decomp_pct == max(decomp_pct)) %>%
              mutate(label = paste(method, length_label, error_label)) %>% pull(label)))
cat(sprintf("Mean decompression fraction: %.1f%%\n", mean(df_decomp_frac$decomp_pct)))
cat(sprintf("Conditions with decomp < 3%%: %d / %d\n",
            sum(df_decomp_frac$decomp_pct < 3), nrow(df_decomp_frac)))

cat("\n")

# =============================================================================
# REAL DATA TABLE STATISTICS (Table 1)
# =============================================================================

cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat("REAL DATA TABLE STATISTICS (Table 1)\n")
cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Helper function to compute method stats from combined benchmark file
# Input: raw data with columns cm, mc, size_cigar_bytes, size_tpa_bytes, etc.
# Output: one row per (cm, mc) combination with aggregated stats
compute_method_stats <- function(data) {
  data %>%
    mutate(
      method = case_when(
        cm == "diagonal-distance" ~ "DB-TP",
        cm == "edit-distance" ~ "EB-TP"
      ),
      method_label = case_when(
        cm == "edit-distance" ~ paste0("EB-TP (\u03b4=", mc, ")"),
        cm == "diagonal-distance" ~ paste0("DB-TP (b=", mc, ")")
      ),
      decode_total_runtime = decompress_runtime_sec + decode_runtime_sec
    ) %>%
    group_by(cm, mc, method, method_label) %>%
    summarise(
      total_cigar_bytes = sum(size_cigar_bytes),
      total_tpa_bytes = sum(size_tpa_bytes),
      compression_ratio = sum(size_tpa_bytes) / sum(size_cigar_bytes),
      total_decode_runtime_sec = sum(decode_total_runtime),
      total_decode_runtime_hrs = sum(decode_total_runtime) / 3600,
      peak_decode_memory_kb = max(decode_memory_kb),
      peak_decode_memory_gib = max(decode_memory_kb) / 1024 / 1024,
      total_alignments = sum(num_alignments),
      total_identical = sum(score_identical),
      total_improved = sum(score_improved),
      total_degraded = sum(score_degraded),
      pct_identical = sum(score_identical) / sum(num_alignments) * 100,
      pct_improved = sum(score_improved) / sum(num_alignments) * 100,
      pct_degraded = sum(score_degraded) / sum(num_alignments) * 100,
      pct_preserved = (sum(score_identical) + sum(score_improved)) / sum(num_alignments) * 100,
      n_files = n(),
      .groups = "drop"
    )
}

# -----------------------------------------------------------------------------
# Human-Human (HPRCv2) dataset
# -----------------------------------------------------------------------------

cat("HUMAN-HUMAN (HPRCv2) DATASET\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

human_data <- read_tsv(file.path(real_data_dir, "hprcv2-25k.benchmark.results.tsv"), show_col_types = FALSE)
human_data <- human_data %>% filter(mc != 16)
human_stats <- compute_method_stats(human_data)

cat("Human summary by method (all thresholds):\n")
print(human_stats %>% select(method_label, total_tpa_bytes, compression_ratio,
                              total_decode_runtime_hrs, peak_decode_memory_gib,
                              pct_identical, pct_improved, pct_preserved))

cat("\n--- TABLE 1 VALUES (Human-Human) ---\n")
cat(sprintf("Total CIGAR size: %.1f GiB\n", human_stats$total_cigar_bytes[1] / 1024^3))

# PAF size from manuscript Table 1 (external measurement)
human_paf_size_gib <- 940.1

cat(sprintf("PAF size (from Table 1): %.1f GiB\n", human_paf_size_gib))

for (i in seq_len(nrow(human_stats))) {
  row <- human_stats[i, ]
  storage_gib <- row$total_tpa_bytes / 1024^3
  ratio_vs_paf <- storage_gib / human_paf_size_gib
  cat(sprintf("\n%s:\n", row$method_label))
  cat(sprintf("  Storage: %.1f GiB\n", storage_gib))
  cat(sprintf("  Compression ratio vs CIGAR: %.3f\u00d7\n", row$compression_ratio))
  cat(sprintf("  Compression ratio vs PAF: %.3f\u00d7 (Table 1 value)\n", ratio_vs_paf))
  cat(sprintf("  TPA->PAF time: %.1f hours (= sum of all files)\n", row$total_decode_runtime_hrs))
  cat(sprintf("  TPA->PAF peak memory: %.0f GiB (= max across files)\n", row$peak_decode_memory_gib))
  cat(sprintf("  Score preserved: %.2f%% (identical: %.2f%%, improved: %.2f%%, degraded: %.2f%%)\n",
              row$pct_preserved, row$pct_identical, row$pct_improved, row$pct_degraded))
}

cat("\n")

# -----------------------------------------------------------------------------
# Primate-Primate (T2T Apes) dataset
# -----------------------------------------------------------------------------

cat("PRIMATE-PRIMATE (T2T APES) DATASET\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

primate_data <- read_tsv(file.path(real_data_dir, "t2t-ape-pangenome.benchmark.results.tsv"), show_col_types = FALSE)
primate_stats <- compute_method_stats(primate_data)

cat("Primate summary by method (all thresholds):\n")
print(primate_stats %>% select(method_label, total_tpa_bytes, compression_ratio,
                                total_decode_runtime_hrs, peak_decode_memory_gib,
                                pct_identical, pct_improved, pct_preserved))

cat("\n--- TABLE 1 VALUES (Primate-Primate) ---\n")
cat(sprintf("Total CIGAR size: %.1f GiB\n", primate_stats$total_cigar_bytes[1] / 1024^3))

# PAF size from manuscript Table 1 (external measurement)
primate_paf_size_gib <- 78.9

cat(sprintf("PAF size (from Table 1): %.1f GiB\n", primate_paf_size_gib))

for (i in seq_len(nrow(primate_stats))) {
  row <- primate_stats[i, ]
  storage_gib <- row$total_tpa_bytes / 1024^3
  ratio_vs_paf <- storage_gib / primate_paf_size_gib
  cat(sprintf("\n%s:\n", row$method_label))
  cat(sprintf("  Storage: %.2f GiB\n", storage_gib))
  cat(sprintf("  Compression ratio vs CIGAR: %.3f\u00d7\n", row$compression_ratio))
  cat(sprintf("  Compression ratio vs PAF: %.3f\u00d7 (Table 1 value)\n", ratio_vs_paf))
  cat(sprintf("  TPA->PAF time: %.1f hours (= sum of all files)\n", row$total_decode_runtime_hrs))
  cat(sprintf("  TPA->PAF peak memory: %.0f GiB (= max across files)\n", row$peak_decode_memory_gib))
  cat(sprintf("  Score preserved: %.2f%% (identical: %.2f%%, improved: %.2f%%, degraded: %.2f%%)\n",
              row$pct_preserved, row$pct_identical, row$pct_improved, row$pct_degraded))
}

cat("\n")

# -----------------------------------------------------------------------------
# Tracepoint counts (for text)
# -----------------------------------------------------------------------------

cat("TRACEPOINT COUNTS AND SEGMENT SIZES (for text)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Tracepoint counts from benchmark files — all EB-TP δ values and DB-TP
# Use total alignment length from aln-stats TOTAL row for approximate segment sizes
human_aln_stats_tp <- read_tsv(file.path(real_data_dir, "hprcv2-25k.aln-stats.tsv"), show_col_types = FALSE)
primate_aln_stats_tp <- read_tsv(file.path(real_data_dir, "t2t-ape-pangenome.aln-stats.tsv"), show_col_types = FALSE)
human_total_aln_length_bp <- (human_aln_stats_tp %>% filter(file == "TOTAL"))$total_length_GB * 1e9
primate_total_aln_length_bp <- (primate_aln_stats_tp %>% filter(file == "TOTAL"))$total_length_GB * 1e9

cat("Human tracepoint counts (all methods/thresholds):\n")
human_tp_counts <- human_data %>%
  mutate(
    method_label = case_when(
      cm == "diagonal-distance" ~ paste0("DB-TP (b=", mc, ")"),
      cm == "edit-distance" ~ paste0("EB-TP (\u03b4=", mc, ")")
    )
  ) %>%
  group_by(cm, mc, method_label) %>%
  summarise(total_tracepoints = sum(num_tracepoints), .groups = "drop") %>%
  mutate(approx_avg_seg_bp = human_total_aln_length_bp / total_tracepoints)
print(human_tp_counts)

# Segment sizes from tp-stats files (weighted average) — only available for mc=32
human_dbtp_stats <- read_tsv(file.path(real_data_dir, "hprcv2-25k.diagonal-distance-32.tp-stats.tsv"), show_col_types = FALSE)
human_ebtp_stats <- read_tsv(file.path(real_data_dir, "hprcv2-25k.edit-distance-32.tp-stats.tsv"), show_col_types = FALSE)

human_dbtp_avg_seg <- sum(human_dbtp_stats$num_pairs * human_dbtp_stats$val1_avg) / sum(human_dbtp_stats$num_pairs)
human_ebtp_avg_seg <- sum(human_ebtp_stats$num_pairs * human_ebtp_stats$val1_avg) / sum(human_ebtp_stats$num_pairs)

cat(sprintf("\nHuman segment sizes (weighted avg from tp-stats, mc=32 only):\n"))
cat(sprintf("  DB-TP: %.1f Kb\n", human_dbtp_avg_seg / 1000))
cat(sprintf("  EB-TP: %.1f Kb\n", human_ebtp_avg_seg / 1000))

cat(sprintf("\nHuman summary (all EB-TP thresholds):\n"))
for (i in seq_len(nrow(human_tp_counts))) {
  row <- human_tp_counts[i, ]
  if (row$approx_avg_seg_bp >= 1000) {
    cat(sprintf("  %s: %.1f billion tracepoints, approx avg segment %.1f Kb\n",
                row$method_label, row$total_tracepoints / 1e9, row$approx_avg_seg_bp / 1000))
  } else {
    cat(sprintf("  %s: %.0f million tracepoints, approx avg segment %.0f bp\n",
                row$method_label, row$total_tracepoints / 1e6, row$approx_avg_seg_bp))
  }
}

# DB-TP vs each EB-TP threshold
human_dbtp_tp <- human_tp_counts %>% filter(cm == "diagonal-distance") %>% pull(total_tracepoints)
human_ebtp_tp <- human_tp_counts %>% filter(cm == "edit-distance")
cat("\nHuman tracepoint ratios (EB-TP / DB-TP):\n")
for (i in seq_len(nrow(human_ebtp_tp))) {
  row <- human_ebtp_tp[i, ]
  cat(sprintf("  %s / DB-TP: %.1f\u00d7\n", row$method_label, row$total_tracepoints / human_dbtp_tp))
}

cat("\nPrimate tracepoint counts (all methods/thresholds):\n")
primate_tp_counts <- primate_data %>%
  mutate(
    method_label = case_when(
      cm == "diagonal-distance" ~ paste0("DB-TP (b=", mc, ")"),
      cm == "edit-distance" ~ paste0("EB-TP (\u03b4=", mc, ")")
    )
  ) %>%
  group_by(cm, mc, method_label) %>%
  summarise(total_tracepoints = sum(num_tracepoints), .groups = "drop") %>%
  mutate(approx_avg_seg_bp = primate_total_aln_length_bp / total_tracepoints)
print(primate_tp_counts)

# Segment sizes from tp-stats files — only available for mc=32
primate_dbtp_stats <- read_tsv(file.path(real_data_dir, "t2t-ape-pangenome.diagonal-distance-32.tp-stats.tsv"), show_col_types = FALSE)
primate_ebtp_stats <- read_tsv(file.path(real_data_dir, "t2t-ape-pangenome.edit-distance-32.tp-stats.tsv"), show_col_types = FALSE)

primate_dbtp_avg_seg <- sum(primate_dbtp_stats$num_pairs * primate_dbtp_stats$val1_avg) / sum(primate_dbtp_stats$num_pairs)
primate_ebtp_avg_seg <- sum(primate_ebtp_stats$num_pairs * primate_ebtp_stats$val1_avg) / sum(primate_ebtp_stats$num_pairs)

cat(sprintf("\nPrimate segment sizes (weighted avg from tp-stats, mc=32 only):\n"))
cat(sprintf("  DB-TP: %.1f Kb\n", primate_dbtp_avg_seg / 1000))
cat(sprintf("  EB-TP: %.0f bp\n", primate_ebtp_avg_seg))

cat(sprintf("\nPrimate summary (all EB-TP thresholds):\n"))
for (i in seq_len(nrow(primate_tp_counts))) {
  row <- primate_tp_counts[i, ]
  if (row$approx_avg_seg_bp >= 1000) {
    cat(sprintf("  %s: %.0f million tracepoints, approx avg segment %.1f Kb\n",
                row$method_label, row$total_tracepoints / 1e6, row$approx_avg_seg_bp / 1000))
  } else {
    cat(sprintf("  %s: %.0f million tracepoints, approx avg segment %.0f bp\n",
                row$method_label, row$total_tracepoints / 1e6, row$approx_avg_seg_bp))
  }
}

# DB-TP vs each EB-TP threshold
primate_dbtp_tp <- primate_tp_counts %>% filter(cm == "diagonal-distance") %>% pull(total_tracepoints)
primate_ebtp_tp <- primate_tp_counts %>% filter(cm == "edit-distance")
cat("\nPrimate tracepoint ratios (EB-TP / DB-TP):\n")
for (i in seq_len(nrow(primate_ebtp_tp))) {
  row <- primate_ebtp_tp[i, ]
  cat(sprintf("  %s / DB-TP: %.1f\u00d7\n", row$method_label, row$total_tracepoints / primate_dbtp_tp))
}

cat("\n")

# -----------------------------------------------------------------------------
# Bits per base (for text)
# -----------------------------------------------------------------------------

cat("BITS PER BASE (for text)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Human alignment stats - use TOTAL row only
human_aln_stats <- read_tsv(file.path(real_data_dir, "hprcv2-25k.aln-stats.tsv"), show_col_types = FALSE)
human_total_row <- human_aln_stats %>% filter(file == "TOTAL")
human_total_bases <- human_total_row$total_length_GB * 1e9  # Convert GB to bases

cat(sprintf("Human total aligned bases: %.1f Tb (from TOTAL row)\n", human_total_bases / 1e12))

human_bpb <- human_stats %>%
  mutate(bits_per_base = (total_tpa_bytes * 8) / human_total_bases)
cat("Human bits per base:\n")
for (i in seq_len(nrow(human_bpb))) {
  row <- human_bpb[i, ]
  cat(sprintf("  %s: %.3f bits/bp\n", row$method_label, row$bits_per_base))
}

# Primate alignment stats - use TOTAL row only
primate_aln_stats <- read_tsv(file.path(real_data_dir, "t2t-ape-pangenome.aln-stats.tsv"), show_col_types = FALSE)
primate_total_row <- primate_aln_stats %>% filter(file == "TOTAL")
primate_total_bases <- primate_total_row$total_length_GB * 1e9

cat(sprintf("\nPrimate total aligned bases: %.1f Gb (from TOTAL row)\n", primate_total_bases / 1e9))

primate_bpb <- primate_stats %>%
  mutate(bits_per_base = (total_tpa_bytes * 8) / primate_total_bases)
cat("Primate bits per base:\n")
for (i in seq_len(nrow(primate_bpb))) {
  row <- primate_bpb[i, ]
  cat(sprintf("  %s: %.3f bits/bp\n", row$method_label, row$bits_per_base))
}

cat("\n")

# =============================================================================
# ADDITIONAL VERIFICATIONS FOR MANUSCRIPT
# =============================================================================

cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat("ADDITIONAL MANUSCRIPT NUMBER VERIFICATIONS\n")
cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# SIMULATED DATA: Heuristic vs Standard WFA (Supp S2/S3)
# Line 34: "runtime by up to 3.7×" and "memory by up to 5.4×"
# -----------------------------------------------------------------------------

cat("SIMULATED: Heuristic vs Standard WFA speedup (Supp S2/S3)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Filter for DB-TP at 100 Kb with mc=32 and high memory mode (where heuristic matters)
heuristic_comparison <- data %>%
  filter(memory_mode == "high") %>%
  filter(tp_type == "standard" & cm == "diagonal-distance" & mc == 32) %>%
  filter(l == 100000) %>%
  mutate(
    runtime_speedup = decode_runtime_sec / decode_heuristic_runtime_sec,
    memory_reduction = decode_memory_kb / decode_heuristic_memory_kb
  ) %>%
  select(l, e, decode_runtime_sec, decode_heuristic_runtime_sec, runtime_speedup,
         decode_memory_kb, decode_heuristic_memory_kb, memory_reduction)

cat("DB-TP at 100 Kb, standard WFA vs heuristic WFA (high memory mode):\n")
print(heuristic_comparison)

cat(sprintf("\nMax runtime speedup (standard / heuristic): %.1f×\n", max(heuristic_comparison$runtime_speedup)))
cat(sprintf("Max memory reduction (standard / heuristic): %.1f×\n", max(heuristic_comparison$memory_reduction)))

# Also check across all lengths and error rates
heuristic_all <- data %>%
  filter(memory_mode == "high") %>%
  filter(tp_type == "standard" & cm == "diagonal-distance" & mc == 32) %>%
  mutate(
    runtime_speedup = decode_runtime_sec / decode_heuristic_runtime_sec,
    memory_reduction = decode_memory_kb / decode_heuristic_memory_kb
  )

cat(sprintf("\nAcross all lengths (high memory mode):\n"))
cat(sprintf("  Max runtime speedup: %.1f×\n", max(heuristic_all$runtime_speedup)))
cat(sprintf("  Max memory reduction: %.1f×\n", max(heuristic_all$memory_reduction)))

cat("\n")

# -----------------------------------------------------------------------------
# SIMULATED DATA: Tracepoint counts at 100 Kb (Supp S4)
# Line 35: FL-TP ~1M, DB-TP 13K at 10% (77× fewer), 2K at 1% (500× fewer), EB-TP 289K at 10%
# -----------------------------------------------------------------------------

cat("SIMULATED: Tracepoint counts at 100 Kb (Supp S4, Line 35)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

tp_counts_100kb <- data %>%
  filter(l == 100000, memory_mode == "high") %>%
  filter((tp_type == "fastga" & mc == 100) | (tp_type == "standard" & mc == 32)) %>%
  mutate(
    method = case_when(
      tp_type == "fastga" ~ "FL-TP",
      cm == "edit-distance" ~ "EB-TP",
      cm == "diagonal-distance" ~ "DB-TP"
    )
  ) %>%
  select(e, method, num_tracepoints) %>%
  pivot_wider(names_from = method, values_from = num_tracepoints)

cat("Tracepoint counts at 100 Kb:\n")
print(tp_counts_100kb)

# FL-TP produces ~1M tracepoints regardless of divergence
cat(sprintf("\nFL-TP tracepoints at 100 Kb: %.0fK -- %.0fK (manuscript: ~1M)\n",
            min(tp_counts_100kb$`FL-TP`) / 1000,
            max(tp_counts_100kb$`FL-TP`) / 1000))

# DB-TP at 10% divergence
dbtp_10pct <- tp_counts_100kb %>% filter(e == 0.1) %>% pull(`DB-TP`)
fltp_10pct <- tp_counts_100kb %>% filter(e == 0.1) %>% pull(`FL-TP`)
cat(sprintf("DB-TP at 10%%: %.0fK tracepoints (manuscript: 13K)\n", dbtp_10pct / 1000))
cat(sprintf("  FL-TP / DB-TP = %.0f× fewer (manuscript: 77×)\n", fltp_10pct / dbtp_10pct))

# DB-TP at 1% divergence
dbtp_1pct <- tp_counts_100kb %>% filter(e == 0.01) %>% pull(`DB-TP`)
fltp_1pct <- tp_counts_100kb %>% filter(e == 0.01) %>% pull(`FL-TP`)
cat(sprintf("DB-TP at 1%%: %.0fK tracepoints (manuscript: 2K)\n", dbtp_1pct / 1000))
cat(sprintf("  FL-TP / DB-TP = %.0f× fewer (manuscript: 500×)\n", fltp_1pct / dbtp_1pct))

# EB-TP at 10% divergence
ebtp_10pct <- tp_counts_100kb %>% filter(e == 0.1) %>% pull(`EB-TP`)
cat(sprintf("EB-TP at 10%%: %.0fK tracepoints (manuscript: 289K)\n", ebtp_10pct / 1000))

cat("\n")

# -----------------------------------------------------------------------------
# REAL DATA: Dataset description from TOTAL rows (Lines 45-46)
# -----------------------------------------------------------------------------

cat("REAL DATA: Dataset description verification (Lines 45-46)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Human dataset from TOTAL row
human_total <- human_aln_stats %>% filter(file == "TOTAL")

cat("Human-Human (HPRCv2) - from TOTAL row:\n")
cat(sprintf("  Alignments: %s (manuscript: 389,234,668)\n", format(human_total$num_alignments, big.mark = ",")))
cat(sprintf("  Total length: %.1f Tb (manuscript: 76.7 Tb)\n", human_total$total_length_GB / 1000))
cat(sprintf("  Mean length: %.0f Kb (manuscript: 197 Kb)\n", human_total$len_avg / 1000))
cat(sprintf("  Max length: %.1f Mb (manuscript: 30.1 Mb)\n", human_total$len_max / 1e6))
cat(sprintf("  Mean identity: %.1f%% ± %.1f%% (manuscript: 98.3%% ± 2.2%%)\n",
            human_total$id_avg * 100, human_total$id_std * 100))

cat("\n")

# Primate dataset from TOTAL row
primate_total <- primate_aln_stats %>% filter(file == "TOTAL")

cat("Primate-Primate (T2T apes) - from TOTAL row:\n")
cat(sprintf("  Alignments: %s (manuscript: 566,038)\n", format(primate_total$num_alignments, big.mark = ",")))
cat(sprintf("  Total length: %.0f Gb (manuscript: 683 Gb)\n", primate_total$total_length_GB))
cat(sprintf("  Mean length: %.2f Mb (manuscript: 1.21 Mb)\n", primate_total$len_avg / 1e6))
cat(sprintf("  Max length: %.0f Mb (manuscript: 170 Mb)\n", primate_total$len_max / 1e6))
cat(sprintf("  Mean identity: %.1f%% ± %.1f%% (manuscript: 91.2%% ± 8.5%%)\n",
            primate_total$id_avg * 100, primate_total$id_std * 100))

cat("\n")

# -----------------------------------------------------------------------------
# REAL DATA: DB-TP vs EB-TP ratio comparison (Line 51)
# "DB-TP ... outperforming EB-TP ... by 1.7× and 3.3× respectively"
# -----------------------------------------------------------------------------

cat("REAL DATA: DB-TP vs EB-TP compression ratio comparison (Line 51)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

for (mc_val in c(32, 64, 128)) {
  cat(sprintf("--- mc=%d ---\n", mc_val))

  human_ebtp_ratio <- human_stats %>% filter(method == "EB-TP", mc == mc_val) %>% pull(compression_ratio)
  human_dbtp_ratio <- human_stats %>% filter(method == "DB-TP", mc == mc_val) %>% pull(compression_ratio)

  primate_ebtp_ratio <- primate_stats %>% filter(method == "EB-TP", mc == mc_val) %>% pull(compression_ratio)
  primate_dbtp_ratio <- primate_stats %>% filter(method == "DB-TP", mc == mc_val) %>% pull(compression_ratio)

  if (length(human_ebtp_ratio) > 0 && length(human_dbtp_ratio) > 0) {
    human_improvement <- human_ebtp_ratio / human_dbtp_ratio
    cat(sprintf("Human: EB-TP ratio=%.3f, DB-TP ratio=%.3f\n", human_ebtp_ratio, human_dbtp_ratio))
    cat(sprintf("  DB-TP outperforms EB-TP by: %.1f\u00d7\n", human_improvement))
  } else {
    cat("  Human: data not yet available for this mc\n")
  }

  if (length(primate_ebtp_ratio) > 0 && length(primate_dbtp_ratio) > 0) {
    primate_improvement <- primate_ebtp_ratio / primate_dbtp_ratio
    cat(sprintf("Primate: EB-TP ratio=%.3f, DB-TP ratio=%.3f\n", primate_ebtp_ratio, primate_dbtp_ratio))
    cat(sprintf("  DB-TP outperforms EB-TP by: %.1f\u00d7\n", primate_improvement))
  } else {
    cat("  Primate: data not yet available for this mc\n")
  }

  cat("\n")
}

cat("\n")

# -----------------------------------------------------------------------------
# REAL DATA: Bits per base verification (Line 72)
# -----------------------------------------------------------------------------

cat("REAL DATA: Bits per base verification (Line 72)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Human total bases from TOTAL row
human_total_bases_bp <- human_total$total_length_GB * 1e9

cat("Human bits per base:\n")
for (i in seq_len(nrow(human_stats))) {
  row <- human_stats[i, ]
  bpb <- (row$total_tpa_bytes * 8) / human_total_bases_bp
  cat(sprintf("  %s: %.3f bits/bp (formula: %.1f GiB \u00d7 8 / %.1f Tb)\n",
              row$method_label, bpb, row$total_tpa_bytes / 1024^3, human_total_bases_bp / 1e12))
}

# Also compute for uncompressed PAF (from Table 1: 940.1 GiB)
paf_size_human <- 940.1 * 1024^3  # from manuscript Table 1
cat(sprintf("  PAF: %.2f bits/bp (from Table 1: 940.1 GiB / %.1f Tb)\n",
            (paf_size_human * 8) / human_total_bases_bp, human_total_bases_bp / 1e12))

cat("\nManuscript claims (Line 72, mc=32): DB-TP 0.003, EB-TP 0.005, PAF 0.10\n")

# Primate total bases from TOTAL row
primate_total_bases_bp <- primate_total$total_length_GB * 1e9

cat("\nPrimate bits per base:\n")
for (i in seq_len(nrow(primate_stats))) {
  row <- primate_stats[i, ]
  bpb <- (row$total_tpa_bytes * 8) / primate_total_bases_bp
  cat(sprintf("  %s: %.3f bits/bp (formula: %.2f GiB \u00d7 8 / %.1f Gb)\n",
              row$method_label, bpb, row$total_tpa_bytes / 1024^3, primate_total_bases_bp / 1e9))
}

cat("\nManuscript claims (Line 72, mc=32): DB-TP 0.007, EB-TP 0.022\n")

cat("\n")

# -----------------------------------------------------------------------------
# REAL DATA: Per-file runtime and speed/memory differences (Line 82)
# -----------------------------------------------------------------------------

cat("REAL DATA: Per-file runtime and speed/memory differences (Line 82)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Per-file runtimes (mc=32 for manuscript claims)
human_perfile <- human_stats %>%
  mutate(avg_runtime_per_file = total_decode_runtime_sec / n_files)

cat("Human per-file statistics:\n")
for (i in seq_len(nrow(human_perfile))) {
  row <- human_perfile[i, ]
  cat(sprintf("  %s: %.0f sec/file (%.1f hours total, %d files)\n",
              row$method_label, row$avg_runtime_per_file, row$total_decode_runtime_hrs, row$n_files))
}
cat("Manuscript claims (mc=32): DB-TP 118 sec/file, EB-TP 54 sec/file\n")

primate_perfile <- primate_stats %>%
  mutate(avg_runtime_per_file = total_decode_runtime_sec / n_files)

cat("\nPrimate per-file statistics:\n")
for (i in seq_len(nrow(primate_perfile))) {
  row <- primate_perfile[i, ]
  cat(sprintf("  %s: %.0f sec/file (%.1f hours total, %d files)\n",
              row$method_label, row$avg_runtime_per_file, row$total_decode_runtime_hrs, row$n_files))
}
cat("Manuscript claims (mc=32): DB-TP 1068 sec/file, EB-TP 48 sec/file\n")

# Speed and memory differences (all EB-TP thresholds vs DB-TP)
cat("\nSpeed and memory differences (DB-TP b=32 vs EB-TP, all thresholds):\n")
primate_dbtp_perfile <- primate_perfile %>% filter(method == "DB-TP", mc == 32)
primate_dbtp_mem <- primate_stats %>% filter(method == "DB-TP", mc == 32) %>% pull(peak_decode_memory_gib)

for (mc_val in c(32, 64, 128)) {
  ebtp_pf <- primate_perfile %>% filter(method == "EB-TP", mc == mc_val)
  ebtp_mem <- primate_stats %>% filter(method == "EB-TP", mc == mc_val) %>% pull(peak_decode_memory_gib)
  if (nrow(ebtp_pf) > 0) {
    speed_diff <- primate_dbtp_perfile$avg_runtime_per_file / ebtp_pf$avg_runtime_per_file
    memory_diff <- primate_dbtp_mem / ebtp_mem
    cat(sprintf("  Primate DB-TP / EB-TP (\u03b4=%d): speed %.1f\u00d7, memory %.1f\u00d7\n",
                mc_val, speed_diff, memory_diff))
  }
}

cat("\n")

# -----------------------------------------------------------------------------
# REAL DATA: Decompression fractions (for text, line 84)
# -----------------------------------------------------------------------------

cat("REAL DATA: Decompression fractions (for text, line 84)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Helper function: compute decompression fraction per method/threshold
compute_decomp_fractions <- function(data, dataset_name) {
  data %>%
    mutate(
      method_label = case_when(
        cm == "diagonal-distance" ~ paste0("DB-TP (b=", mc, ")"),
        cm == "edit-distance" ~ paste0("EB-TP (\u03b4=", mc, ")")
      ),
      total_runtime = decompress_runtime_sec + decode_runtime_sec,
      decomp_frac = ifelse(total_runtime > 0, decompress_runtime_sec / total_runtime * 100, 0)
    ) %>%
    group_by(cm, mc, method_label) %>%
    summarise(
      min_decomp_pct = min(decomp_frac),
      max_decomp_pct = max(decomp_frac),
      mean_decomp_pct = mean(decomp_frac),
      .groups = "drop"
    ) %>%
    mutate(dataset = dataset_name)
}

human_decomp_fracs <- compute_decomp_fractions(human_data, "Human")
primate_decomp_fracs <- compute_decomp_fractions(primate_data, "Primate")

cat("Human decompression fractions:\n")
for (i in seq_len(nrow(human_decomp_fracs))) {
  row <- human_decomp_fracs[i, ]
  cat(sprintf("  %s: %.1f--%.1f%% (mean %.1f%%)\n",
              row$method_label, row$min_decomp_pct, row$max_decomp_pct, row$mean_decomp_pct))
}

cat("\nPrimate decompression fractions:\n")
for (i in seq_len(nrow(primate_decomp_fracs))) {
  row <- primate_decomp_fracs[i, ]
  cat(sprintf("  %s: %.1f--%.1f%% (mean %.1f%%)\n",
              row$method_label, row$min_decomp_pct, row$max_decomp_pct, row$mean_decomp_pct))
}

cat("\n")

# -----------------------------------------------------------------------------
# BGZIP DECOMPRESSION BENCHMARKS (for Table)
# -----------------------------------------------------------------------------

cat("BGZIP DECOMPRESSION BENCHMARKS (for Table)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Human BGZIP decompression
human_bgzip <- read_tsv(file.path(real_data_dir, "hprcv2-25k.bgzip-decompression.benchmark.tsv"), show_col_types = FALSE)

human_bgzip_total_runtime_h <- sum(human_bgzip$decompress_runtime_sec, na.rm = TRUE) / 3600
human_bgzip_peak_memory_gib <- max(human_bgzip$decompress_memory_kb, na.rm = TRUE) / 1024 / 1024
human_bgzip_total_size_gib <- sum(human_bgzip$size_gz_bytes, na.rm = TRUE) / 1024^3

cat(sprintf("Human BGZIP decompression:\n"))
cat(sprintf("  Files: %d (with runtime data: %d)\n",
            nrow(human_bgzip), sum(!is.na(human_bgzip$decompress_runtime_sec))))
cat(sprintf("  Total BGZIP size: %.1f GiB\n", human_bgzip_total_size_gib))
cat(sprintf("  Total decompress runtime: %.2f hours (%.1f sec)\n",
            human_bgzip_total_runtime_h, sum(human_bgzip$decompress_runtime_sec, na.rm = TRUE)))
cat(sprintf("  Peak decompress memory: %.3f GiB (%.1f MB)\n",
            human_bgzip_peak_memory_gib, max(human_bgzip$decompress_memory_kb, na.rm = TRUE) / 1024))
cat(sprintf("  Avg runtime per file: %.1f sec\n",
            mean(human_bgzip$decompress_runtime_sec, na.rm = TRUE)))

# Primate BGZIP decompression
primate_bgzip <- read_tsv(file.path(real_data_dir, "t2t-ape-pangenome.bgzip-decompression.benchmark.tsv"), show_col_types = FALSE)

primate_bgzip_total_runtime_h <- sum(primate_bgzip$decompress_runtime_sec, na.rm = TRUE) / 3600
primate_bgzip_peak_memory_gib <- max(primate_bgzip$decompress_memory_kb, na.rm = TRUE) / 1024 / 1024
primate_bgzip_total_size_gib <- sum(primate_bgzip$size_gz_bytes, na.rm = TRUE) / 1024^3

cat(sprintf("\nPrimate BGZIP decompression:\n"))
cat(sprintf("  Files: %d (with runtime data: %d)\n",
            nrow(primate_bgzip), sum(!is.na(primate_bgzip$decompress_runtime_sec))))
cat(sprintf("  Total BGZIP size: %.1f GiB\n", primate_bgzip_total_size_gib))
cat(sprintf("  Total decompress runtime: %.2f hours (%.1f sec)\n",
            primate_bgzip_total_runtime_h, sum(primate_bgzip$decompress_runtime_sec, na.rm = TRUE)))
cat(sprintf("  Peak decompress memory: %.3f GiB (%.1f MB)\n",
            primate_bgzip_peak_memory_gib, max(primate_bgzip$decompress_memory_kb, na.rm = TRUE) / 1024))
cat(sprintf("  Avg runtime per file: %.1f sec\n",
            mean(primate_bgzip$decompress_runtime_sec, na.rm = TRUE)))

cat(sprintf("\n--- TABLE VALUES (BGZIP decompression) ---\n"))
cat(sprintf("Human:   time = %.1f h,  peak mem = %.3f GiB (< 1 GiB)\n",
            human_bgzip_total_runtime_h, human_bgzip_peak_memory_gib))
cat(sprintf("Primate: time = %.1f h,  peak mem = %.3f GiB (< 1 GiB)\n",
            primate_bgzip_total_runtime_h, primate_bgzip_peak_memory_gib))

# Storage comparison: BGZIP vs adaptive TPA (computed from stats)
human_bgzip_storage <- human_bgzip_total_size_gib
primate_bgzip_storage <- primate_bgzip_total_size_gib

cat(sprintf("\nBGZIP storage vs TPA (all methods):\n"))
for (i in seq_len(nrow(human_stats))) {
  row <- human_stats[i, ]
  storage_gib <- row$total_tpa_bytes / 1024^3
  cat(sprintf("  Human:   BGZIP / %s = %.1f\u00d7\n", row$method_label, human_bgzip_storage / storage_gib))
}
for (i in seq_len(nrow(primate_stats))) {
  row <- primate_stats[i, ]
  storage_gib <- row$total_tpa_bytes / 1024^3
  cat(sprintf("  Primate: BGZIP / %s = %.1f\u00d7\n", row$method_label, primate_bgzip_storage / storage_gib))
}

cat("\n")

# -----------------------------------------------------------------------------
# SIMULATED DATA: BGZIP decompression stats (for text)
# -----------------------------------------------------------------------------

cat("SIMULATED: BGZIP decompression stats (for text)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

bgzip_sim <- data %>%
  filter(memory_mode == "high", tp_type == "fastga", mc == 100) %>%
  select(l, e, bgzip_decompress_runtime_sec, bgzip_decompress_memory_kb)

cat("BGZIP decompression at 100 Kb:\n")
bgzip_100kb <- bgzip_sim %>% filter(l == 100000)
print(bgzip_100kb)
cat(sprintf("\n  Runtime range: %.2f -- %.2f sec\n",
            min(bgzip_100kb$bgzip_decompress_runtime_sec),
            max(bgzip_100kb$bgzip_decompress_runtime_sec)))
cat(sprintf("  Peak memory range: %.0f -- %.0f MB\n",
            min(bgzip_100kb$bgzip_decompress_memory_kb) / 1024,
            max(bgzip_100kb$bgzip_decompress_memory_kb) / 1024))

cat("\n")

# =============================================================================
# ABSTRACT PLACEHOLDER NUMBERS
# =============================================================================

cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat("ABSTRACT PLACEHOLDER NUMBERS\n")
cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

cat("Template: Adaptive tracepoints achieve \\red{X}x better compression\n")
cat("  than fixed-length encodings (\\red{Y}x vs \\red{Z}x total compression)\n")
cat("  with minimal score deviation \\red{(<W% in P% of cases)}\n\n")

# --- X, Y, Z: Compression factors from simulated data (100 Kb) ---
# FL-TP is only benchmarked on simulated data, so comparison comes from there.
# Compression factor = CIGAR_size / TPA_size (larger = better compression).

abstract_comp <- df %>%
  filter(l == 100000) %>%
  mutate(compression_factor = size_cigar_bytes / size_tpa_bytes) %>%
  select(e, method, compression_factor) %>%
  pivot_wider(names_from = method, values_from = compression_factor)

abstract_comp <- abstract_comp %>%
  mutate(
    dbtp_vs_fltp = `DB-TP` / `FL-TP`,
    ebtp_vs_fltp = `EB-TP` / `FL-TP`
  )

cat("Compression FACTORS at 100 Kb (CIGAR / TPA, larger = better):\n")
print(abstract_comp)

cat(sprintf("\n  Z (FL-TP factor): %.0f--%.0fx\n",
            min(abstract_comp$`FL-TP`), max(abstract_comp$`FL-TP`)))
cat(sprintf("  Y (DB-TP factor): %.0f--%.0fx\n",
            min(abstract_comp$`DB-TP`), max(abstract_comp$`DB-TP`)))
cat(sprintf("  X (DB-TP / FL-TP improvement): %.0f--%.0fx\n",
            min(abstract_comp$dbtp_vs_fltp), max(abstract_comp$dbtp_vs_fltp)))
cat(sprintf("  EB-TP factor: %.0f--%.0fx\n",
            min(abstract_comp$`EB-TP`), max(abstract_comp$`EB-TP`)))

cat("\nSuggested abstract values (pick representative error rate or range):\n")
cat("  At 1%% error:  DB-TP %.0fx vs FL-TP %.0fx => %.0fx better\n" %>%
      sprintf(abstract_comp %>% filter(e == 0.01) %>% pull(`DB-TP`),
              abstract_comp %>% filter(e == 0.01) %>% pull(`FL-TP`),
              abstract_comp %>% filter(e == 0.01) %>% pull(dbtp_vs_fltp)))
cat("  At 5%% error:  DB-TP %.0fx vs FL-TP %.0fx => %.0fx better\n" %>%
      sprintf(abstract_comp %>% filter(e == 0.05) %>% pull(`DB-TP`),
              abstract_comp %>% filter(e == 0.05) %>% pull(`FL-TP`),
              abstract_comp %>% filter(e == 0.05) %>% pull(dbtp_vs_fltp)))
cat("  At 10%% error: DB-TP %.0fx vs FL-TP %.0fx => %.0fx better\n" %>%
      sprintf(abstract_comp %>% filter(e == 0.1) %>% pull(`DB-TP`),
              abstract_comp %>% filter(e == 0.1) %>% pull(`FL-TP`),
              abstract_comp %>% filter(e == 0.1) %>% pull(dbtp_vs_fltp)))

cat("\n")

# --- W, P: Score preservation from real data ---
# Combine human and primate datasets for overall stats.

cat("Score preservation across REAL datasets:\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Per-dataset, per-method
score_summary <- bind_rows(
  human_stats %>% mutate(dataset = "Human"),
  primate_stats %>% mutate(dataset = "Primate")
) %>%
  select(dataset, method, method_label, mc, total_alignments, total_identical, total_improved,
         total_degraded, pct_identical, pct_improved, pct_degraded, pct_preserved)

cat("Per dataset and method:\n")
print(score_summary %>% select(dataset, method_label, total_alignments,
                                pct_identical, pct_improved, pct_degraded, pct_preserved))

# Overall across both datasets and both methods
overall_score <- score_summary %>%
  summarise(
    total_alignments = sum(total_alignments),
    total_identical = sum(total_identical),
    total_improved = sum(total_improved),
    total_degraded = sum(total_degraded),
    pct_degraded = sum(total_degraded) / sum(total_alignments) * 100,
    pct_preserved = (sum(total_identical) + sum(total_improved)) / sum(total_alignments) * 100
  )

cat(sprintf("\nOverall (all datasets, all methods):\n"))
cat(sprintf("  Total alignments: %s\n", format(overall_score$total_alignments, big.mark = ",")))
cat(sprintf("  Degraded: %s (%.4f%%)\n",
            format(overall_score$total_degraded, big.mark = ","), overall_score$pct_degraded))
cat(sprintf("  Preserved (identical + improved): %.2f%%\n", overall_score$pct_preserved))

# Per-method across datasets
per_method_score <- score_summary %>%
  group_by(method, mc, method_label) %>%
  summarise(
    total_alignments = sum(total_alignments),
    total_degraded = sum(total_degraded),
    pct_degraded = sum(total_degraded) / sum(total_alignments) * 100,
    pct_preserved = (sum(total_identical) + sum(total_improved)) / sum(total_alignments) * 100,
    .groups = "drop"
  )

cat("\nPer method (across datasets):\n")
print(per_method_score)

cat(sprintf("\nSuggested abstract values for score deviation:\n"))
cat(sprintf("  Max degradation %%: %.2f%%\n", max(score_summary$pct_degraded)))
cat(sprintf("  Min preserved %%:   %.2f%%\n", min(score_summary$pct_preserved)))
cat(sprintf("  => (<%.1f%% degradation in %.1f%% of cases)\n",
            max(score_summary$pct_degraded),
            min(score_summary$pct_preserved)))

cat("\n")

# =============================================================================
# ABSTRACT NUMBERS VERIFICATION
# =============================================================================

cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat("ABSTRACT NUMBERS VERIFICATION\n")
cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# --- Simulated data (100 Kb): DB-TP vs FL-TP and EB-TP vs FL-TP ---
cat("SIMULATED (100 Kb): DB-TP and EB-TP vs FL-TP\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat(sprintf("  DB-TP vs FL-TP: %.1f--%.1f× (abstract: 10--13×)\n",
            min(abstract_comp$dbtp_vs_fltp), max(abstract_comp$dbtp_vs_fltp)))
cat(sprintf("  EB-TP vs FL-TP: %.1f--%.1f× (abstract: 1.3--5×)\n",
            min(abstract_comp$ebtp_vs_fltp), max(abstract_comp$ebtp_vs_fltp)))

# --- Real data: compression factor vs uncompressed PAF ---
cat("\nREAL DATA: compression factor vs uncompressed PAF\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n")

real_data_factors <- bind_rows(
  human_stats %>%
    mutate(dataset = "Human",
           compression_factor_vs_paf = 1 / (total_tpa_bytes / 1024^3 / human_paf_size_gib)),
  primate_stats %>%
    mutate(dataset = "Primate",
           compression_factor_vs_paf = 1 / (total_tpa_bytes / 1024^3 / primate_paf_size_gib))
) %>%
  select(dataset, method, method_label, compression_factor_vs_paf)

cat("Per dataset and method:\n")
for (i in seq_len(nrow(real_data_factors))) {
  row <- real_data_factors[i, ]
  cat(sprintf("  %s %s: %.1f\u00d7\n", row$dataset, row$method_label, row$compression_factor_vs_paf))
}

cat(sprintf("\n  Combined range: %.1f--%.1f\u00d7 (abstract: 23--143\u00d7)\n",
            min(real_data_factors$compression_factor_vs_paf),
            max(real_data_factors$compression_factor_vs_paf)))

# --- Total real-data alignments ---
total_real_alignments <- sum(human_stats$total_alignments[1], primate_stats$total_alignments[1])
cat(sprintf("\n  Total real-data alignments: %s (abstract: 780M)\n",
            format(total_real_alignments, big.mark = ",")))

cat("\n")

# =============================================================================
# CROSS-REFERENCE: Line 34 manuscript claims (simulated data, 100 Kb)
# =============================================================================

cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat("CROSS-REFERENCE: LINE 34 MANUSCRIPT CLAIMS\n")
cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

# Claim: "requires 13--144× more storage than adaptive tracepoint formats at 100 Kb"
bgzip_vs_adaptive_100kb <- comparison_100kb %>%
  left_join(comparison_ebtp_100kb %>% select(e, ebtp_ratio), by = "e") %>%
  mutate(
    bgzip_vs_dbtp = bgzip_ratio / dbtp_ratio,
    bgzip_vs_ebtp = bgzip_ratio / ebtp_ratio
  )

cat("Verify: BGZIP vs adaptive TPA at 100 Kb (manuscript: 13--144×):\n")
cat(sprintf("  BGZIP/DB-TP: %.0f--%.0f×\n",
            min(bgzip_vs_adaptive_100kb$bgzip_vs_dbtp),
            max(bgzip_vs_adaptive_100kb$bgzip_vs_dbtp)))
cat(sprintf("  BGZIP/EB-TP: %.0f--%.0f×\n",
            min(bgzip_vs_adaptive_100kb$bgzip_vs_ebtp),
            max(bgzip_vs_adaptive_100kb$bgzip_vs_ebtp)))
cat(sprintf("  Combined range: %.0f--%.0f×\n",
            min(c(bgzip_vs_adaptive_100kb$bgzip_vs_dbtp, bgzip_vs_adaptive_100kb$bgzip_vs_ebtp)),
            max(c(bgzip_vs_adaptive_100kb$bgzip_vs_dbtp, bgzip_vs_adaptive_100kb$bgzip_vs_ebtp))))

# Claim: speedup vs re-alignment at 100 Kb
speedup_all <- c(speedup_table$`FL-TP`, speedup_table$`EB-TP`, speedup_table$`DB-TP`)
cat(sprintf("\nVerify: speedup vs re-alignment at 100 Kb (manuscript: up to 117×):\n"))
cat(sprintf("  Range (all methods, all error rates): %.1f--%.1f×\n",
            min(speedup_all), max(speedup_all)))
cat(sprintf("  Range (5--20%% error only): %.1f--%.1f×\n",
            min(speedup_all[rep(speedup_table$e, 3) >= 0.05]),
            max(speedup_all[rep(speedup_table$e, 3) >= 0.05])))

# Claim: "under 2 seconds at 100 Kb" (BGZIP)
cat(sprintf("\nVerify: BGZIP decompression at 100 Kb (manuscript: under 2 sec):\n"))
cat(sprintf("  Range: %.2f--%.2f sec\n",
            min(bgzip_100kb$bgzip_decompress_runtime_sec),
            max(bgzip_100kb$bgzip_decompress_runtime_sec)))

# Claim: "memory by up to 5.4×" (heuristic WFA)
cat(sprintf("\nVerify: DB-TP memory reduction via heuristic WFA (manuscript: 5.4×):\n"))
cat(sprintf("  At 100 Kb: max %.1f×\n", max(heuristic_comparison$memory_reduction)))
cat(sprintf("  All lengths: max %.1f×\n", max(heuristic_all$memory_reduction)))

cat("\n")

# =============================================================================
# CROSS-REFERENCE: Line 35 decompression fraction
# =============================================================================

cat("CROSS-REFERENCE: LINE 35 DECOMPRESSION FRACTION\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

cat(sprintf("Worst case: %.1f%% (manuscript: less than 20%%)\n", max(df_decomp_frac$decomp_pct)))
cat(sprintf("Conditions < 3%%: %d/%d = %.0f%% (manuscript: 'at most 3%%')\n",
            sum(df_decomp_frac$decomp_pct < 3), nrow(df_decomp_frac),
            sum(df_decomp_frac$decomp_pct < 3) / nrow(df_decomp_frac) * 100))

cat("\n")

# =============================================================================
# CROSS-REFERENCE: Line 85 BGZIP vs adaptive TPA (real data)
# =============================================================================

cat("CROSS-REFERENCE: LINE 85 BGZIP VS ADAPTIVE TPA (REAL DATA)\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

real_bgzip_vs_tpa <- c(
  human_bgzip_storage / (human_stats$total_tpa_bytes / 1024^3),
  primate_bgzip_storage / (primate_stats$total_tpa_bytes / 1024^3)
)

cat("Per method:\n")
for (i in seq_len(nrow(human_stats))) {
  row <- human_stats[i, ]
  storage_gib <- row$total_tpa_bytes / 1024^3
  cat(sprintf("  Human   BGZIP/%s = %.1f\u00d7\n", row$method_label, human_bgzip_storage / storage_gib))
}
for (i in seq_len(nrow(primate_stats))) {
  row <- primate_stats[i, ]
  storage_gib <- row$total_tpa_bytes / 1024^3
  cat(sprintf("  Primate BGZIP/%s = %.1f\u00d7\n", row$method_label, primate_bgzip_storage / storage_gib))
}
cat(sprintf("Range: %.0f--%.0f\u00d7 (manuscript: 7--XXX\u00d7)\n",
            min(real_bgzip_vs_tpa), max(real_bgzip_vs_tpa)))

cat("\n")

# =============================================================================
# CROSS-REFERENCE: Discussion placeholder values
# =============================================================================

cat("CROSS-REFERENCE: DISCUSSION PLACEHOLDER VALUES\n")
cat("-" %>% rep(70) %>% paste0(collapse = ""), "\n\n")

cat("Line 5: 'achieves X× better compression (Y× vs Z× total compression)'\n")
cat(sprintf("  X (DB-TP vs FL-TP improvement at 100 Kb): %.1f--%.1f×\n",
            min(abstract_comp$dbtp_vs_fltp), max(abstract_comp$dbtp_vs_fltp)))
cat(sprintf("  Y (DB-TP compression factor at 100 Kb): %.0f--%.0f×\n",
            min(abstract_comp$`DB-TP`), max(abstract_comp$`DB-TP`)))
cat(sprintf("  Z (FL-TP compression factor at 100 Kb): %.0f--%.0f×\n",
            min(abstract_comp$`FL-TP`), max(abstract_comp$`FL-TP`)))

cat(sprintf("\nLine 13: 'Real-world compression ratios of X× ...'\n"))
cat(sprintf("  Range (all methods, both datasets vs PAF): %.0f--%.0f×\n",
            min(real_data_factors$compression_factor_vs_paf),
            max(real_data_factors$compression_factor_vs_paf)))

cat("\n")
cat("=" %>% rep(70) %>% paste0(collapse = ""), "\n")
cat("Script completed.\n")
