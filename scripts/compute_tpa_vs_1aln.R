#!/usr/bin/env Rscript
# TPA+index vs FASTGA .1aln.
#
# .1aln is a binary ONEcode file whose footer carries the byte index for random access
# (ONElib '&' lines), so its size already includes its index. TPA keeps the index in a
# separate .tpa.idx, which the reported storage excludes. Comparing a bare .tpa against
# a .1aln is therefore not like for like: this script adds the measured index back.
#
# Inputs (all committed):
#   data/real-data/fltp-tpa-index-sizes.tsv        idx per target, real pangenomes
#   data/simulated-data/fltp-tpa-index-sizes.tsv   idx per condition, simulated
#   data/real-data/*.benchmark.results.tsv         1aln sizes (fastga-native*)
#   data/simulated-data/benchmark.results.tsv      1aln sizes (fastga-native*)
suppressPackageStartupMessages(library(tidyverse))
d <- Sys.getenv("TRACEPOINTS_DATA_DIR", unset = "data")

pair <- tibble(tpa_cm = c("fastga", "fastga-no-diff"),
               aln_cm = c("fastga-native", "fastga-native-nodiff"),
               method = c("FL-TP", "FL-TP nd"))

# ---- real pangenomes -------------------------------------------------------------
idx_real <- read_tsv(file.path(d, "real-data/fltp-tpa-index-sizes.tsv"), show_col_types = FALSE) %>%
  group_by(dataset, tpa_cm = method) %>%
  summarise(tpa = sum(tpa_bytes), idx = sum(idx_bytes), .groups = "drop")

aln_real <- bind_rows(
  read_tsv(file.path(d, "real-data/hprcv2-25k.benchmark.results.tsv"), show_col_types = FALSE) %>%
    mutate(dataset = "hprcv2-25k"),
  read_tsv(file.path(d, "real-data/t2t-ape-pangenome.benchmark.results.tsv"), show_col_types = FALSE) %>%
    mutate(dataset = "t2t-ape-pangenome")) %>%
  filter(cm %in% pair$aln_cm, mc == 100) %>%
  group_by(dataset, aln_cm = cm) %>% summarise(aln = sum(as.numeric(size_tpa_bytes)), .groups = "drop")

real <- idx_real %>% left_join(pair, by = "tpa_cm") %>% left_join(aln_real, by = c("dataset", "aln_cm"))

# ---- simulated -------------------------------------------------------------------
idx_sim <- read_tsv(file.path(d, "simulated-data/fltp-tpa-index-sizes.tsv"), show_col_types = FALSE) %>%
  select(l, e, tpa_cm = tp_type, tpa = tpa_bytes, idx = idx_bytes)

aln_sim <- read_tsv(file.path(d, "simulated-data/benchmark.results.tsv"), show_col_types = FALSE) %>%
  filter(tp_type %in% pair$aln_cm, mc == 100) %>%
  select(l, e, aln_cm = tp_type, aln = size_tpa_bytes)

sim <- idx_sim %>% left_join(pair, by = "tpa_cm") %>% left_join(aln_sim, by = c("l", "e", "aln_cm"))

fmt <- function(x) sprintf("%.2f", x)
cat("=== TPA vs 1aln, with TPA's index counted ===\n\n")
cat("REAL PANGENOMES (totals, GiB)\n")
real %>% mutate(across(c(tpa, idx, aln), ~ .x / 1024^3),
                `TPA` = fmt(tpa), `+idx` = fmt(tpa + idx), `1aln` = fmt(aln),
                `1aln/TPA` = fmt(aln / tpa), `1aln/(TPA+idx)` = fmt(aln / (tpa + idx)),
                `idx %` = sprintf("%.1f%%", 100 * idx / tpa)) %>%
  select(dataset, method, TPA, `+idx`, `1aln`, `1aln/TPA`, `1aln/(TPA+idx)`, `idx %`) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nSIMULATED (per condition, ratio range across the 20 conditions)\n")
sim %>% mutate(r_bare = aln / tpa, r_fair = aln / (tpa + idx)) %>%
  group_by(method) %>%
  summarise(n = n(),
            `1aln/TPA` = sprintf("%.2f-%.2f", min(r_bare), max(r_bare)),
            `1aln/(TPA+idx)` = sprintf("%.2f-%.2f", min(r_fair), max(r_fair)),
            `idx %` = sprintf("%.1f-%.1f%%", 100 * min(idx / tpa), 100 * max(idx / tpa)),
            .groups = "drop") %>% as.data.frame() %>% print(row.names = FALSE)

cat(sprintf("\nSentence values -- simulated: 1aln stores %.2f-%.2fx more than TPA+idx\n",
            min(sim$aln / (sim$tpa + sim$idx)), max(sim$aln / (sim$tpa + sim$idx))))
cat(sprintf("                  pangenomes: 1aln stores %.2f-%.2fx more than TPA+idx\n",
            min(real$aln / (real$tpa + real$idx)), max(real$aln / (real$tpa + real$idx))))
