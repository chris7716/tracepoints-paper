#!/usr/bin/env Rscript

# Space-time tradeoff figure (reviewer request):
#   x = stored file size (bytes, log), y = reconstruction time (stored format -> PAF),
#   one connected line per tracepoint method across its parameter sweep
#   (FL-TP spacing l, EB-TP delta, DB-TP b), plus the two CIGAR baselines as single
#   points per panel (plain CIGAR and BGZIP-compressed CIGAR), faceted by sequence
#   length x error rate.
#
# The y-axis uses a pseudo-log (symlog) scale so the baselines, whose reconstruction
# time is ~0 (plain CIGAR is already the alignment; BGZIP only decompresses), can be
# shown at 0 while the tracepoint times keep their log-like spread.
#
# Input:  data/simulated-data/benchmark.results.tsv
# Output: paper/figures/figS8_space_time_tradeoff.{png,pdf}
# Deps:   tidyverse, ggplot2, scales, ggrepel, ggh4x

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(ggrepel)
  library(ggh4x)
})

data_dir <- Sys.getenv("TRACEPOINTS_DATA_DIR", unset = "data")
fig_dir  <- Sys.getenv("TRACEPOINTS_FIG_DIR", unset = "paper/figures")

tsv <- file.path(data_dir, "simulated-data", "benchmark.results.tsv")
if (!file.exists(tsv)) {
  stop("Simulated TSV not found: ", tsv,
       "\nRun the simulated-data benchmark block in code.md first.")
}
message("Reading: ", tsv)
data <- read_tsv(tsv, show_col_types = FALSE)

# Colors/shapes match the other paper figures.
lev  <- c("Plain CIGAR", "CG BGZIP", "FL-TP TPA", "FL-TP TPA nd", "FL-TP 1aln", "EB-TP TPA", "DB-TP TPA")
cols <- c("FL-TP TPA" = "#7570b3", "FL-TP TPA nd" = "#54278f", "FL-TP 1aln" = "#9e9ac8",
          "EB-TP TPA" = "#d95f02", "DB-TP TPA" = "#377eb8", "CG BGZIP" = "#1b9e77", "Plain CIGAR" = "#555555")
shp  <- c("FL-TP TPA" = 16, "FL-TP TPA nd" = 15, "FL-TP 1aln" = 17, "EB-TP TPA" = 16,
          "DB-TP TPA" = 16, "CG BGZIP" = 18, "Plain CIGAR" = 4)

mk_labels <- function(d) d %>% mutate(
  length_label = factor(case_when(
    l == 100 ~ "100 bp", l == 1000 ~ "1 Kbp",
    l == 10000 ~ "10 Kbp", l == 100000 ~ "100 Kbp"),
    levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp")),
  error_label = factor(paste0(e * 100, "%"),
                       levels = c("0.1%", "1%", "5%", "10%", "20%")))

# Tracepoint parameter sweeps. recon_time = stored-format -> PAF (decompress + decode).
# Drop the standard edit-distance mc=16 bootstrap warm-up rows and the single-point
# FASTGA .1aln reference (l=100 only).
tp <- data %>%
  filter(decompress_correct == TRUE | decompress_correct == "true") %>%
  filter(!(tp_type == "standard" & mc == 16)) %>%
  filter(tp_type != "fastga-native") %>%
  mutate(method = case_when(
           tp_type == "fastga"         ~ "FL-TP TPA",
           tp_type == "fastga-no-diff" ~ "FL-TP TPA nd",
           cm == "edit-distance"       ~ "EB-TP TPA",
           cm == "diagonal-distance"   ~ "DB-TP TPA"),
         size = size_tpa_bytes,
         time = decompress_runtime_sec + decode_runtime_sec) %>%
  mk_labels() %>% arrange(method, mc) %>%
  transmute(l, e, length_label, error_label, method, mc, size, time, kind = "sweep")

# CIGAR baselines: one point per (length, error). Plain CIGAR needs no reconstruction
# (time = 0); BGZIP only decompresses (bgzip_decompress_runtime_sec).
base <- data %>% filter(memory_mode == "high") %>% distinct(l, e, .keep_all = TRUE)
cig <- base %>% transmute(l, e, method = "Plain CIGAR", mc = NA_real_,
                          size = size_cigar_bytes, time = 0, kind = "base") %>% mk_labels()
bgz <- base %>% transmute(l, e, method = "CG BGZIP", mc = NA_real_,
                          size = size_cigar_bgzip_bytes,
                          time = bgzip_decompress_runtime_sec, kind = "base") %>% mk_labels()

# FASTGA reference: a single FL-TP 1aln operating point per panel.
# Keep its trace spacing in mc (=100) so the point is labelled like the FL-TP TPA sweep.
aln <- data %>% filter(tp_type == "fastga-native") %>%
  mutate(method = "FL-TP 1aln", size = size_tpa_bytes,
         time = decompress_runtime_sec + decode_runtime_sec, kind = "base") %>%
  mk_labels() %>% transmute(l, e, length_label, error_label, method, mc, size, time, kind)

allpts  <- bind_rows(tp, aln, cig, bgz) %>% mutate(method = factor(method, levels = lev))
sweep   <- filter(allpts, kind == "sweep")
labeled <- filter(allpts, !is.na(mc))   # sweeps + the FL-TP 1aln point

p <- ggplot(allpts, aes(x = size, y = time, color = method, shape = method)) +
  geom_line(data = sweep, aes(group = method), linewidth = 0.5, alpha = 0.8) +
  geom_point(size = 1.8) +
  ggrepel::geom_text_repel(data = labeled, aes(label = mc), size = 2.0, show.legend = FALSE,
                           min.segment.length = 0, segment.size = 0.2, segment.alpha = 0.4,
                           box.padding = 0.45, point.padding = 0.25,
                           force = 4, force_pull = 0.4,
                           max.overlaps = Inf, max.iter = 100000, max.time = 3, seed = 42) +
  ggh4x::facet_grid2(error_label ~ length_label, scales = "free", independent = "all") +
  scale_x_log10(breaks = scales::breaks_log(n = 11), labels = label_bytes()) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.05),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100)) +
  scale_color_manual(values = cols, name = "Method", limits = lev) +
  scale_shape_manual(values = shp, name = "Method", limits = lev) +
  labs(x = "File size", y = "Reconstruction time (sec)") +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "figS8_space_time_tradeoff.png"), p, width = 16, height = 15, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS8_space_time_tradeoff.pdf"), p, width = 16, height = 15, bg = "white")
message("Figure saved: ", file.path(fig_dir, "figS8_space_time_tradeoff.{png,pdf}"))
