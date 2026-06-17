#!/usr/bin/env Rscript

# Space-time tradeoff figure (reviewer request):
#   x = TPA size (bytes), y = reconstruction time (decompress + decode, exact/--no-banded),
#   one connected line per method across its parameter sweep (FL-TP spacing, EB-TP delta, DB-TP b),
#   faceted by sequence length x error rate.
#
# Input:  data/simulated-data/benchmark.results.tsv
#   methods (tp_type): FL-TP (fastga, l 100..2000), FL-TP FASTGA (fastga-native, l=100 only),
#   EB-TP (edit-distance, 32..1024), DB-TP (diagonal-distance, 32..1024)
# Output: paper/figures/space_time_tradeoff.{png,pdf}
# Deps:   tidyverse, ggplot2, scales

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

# Colors and labels match the other paper figures.
# FL-TP 1aln is a single operating point (l=100), not a parameter sweep, and at short
# lengths it sits ~2x to the right of the TPA cluster, squeezing the panels; we omit it
# here so the three swept methods fill each panel.
method_colors <- c("FL-TP TPA" = "#7570b3", "EB-TP TPA" = "#d95f02", "DB-TP TPA" = "#377eb8")
method_levels <- c("FL-TP TPA", "EB-TP TPA", "DB-TP TPA")

mk_labels <- function(d) {
  d %>% mutate(
    length_label = factor(case_when(
      l == 100    ~ "100 bp",
      l == 1000   ~ "1 Kbp",
      l == 10000  ~ "10 Kbp",
      l == 100000 ~ "100 Kbp"),
      levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp")),
    error_label = factor(paste0(e * 100, "%"),
                         levels = c("0.1%", "1%", "5%", "10%", "20%"))
  )
}

df <- data %>%
  filter(decompress_correct == TRUE | decompress_correct == "true") %>%
  # drop per-condition bootstrap warm-up rows (standard edit-distance mc=16, excluded from reporting)
  filter(!(tp_type == "standard" & mc == 16)) %>%
  # omit the single-point FASTGA reference (l=100) so the swept methods fill each panel
  filter(tp_type != "fastga-native") %>%
  mutate(
    method = factor(case_when(
      tp_type == "fastga"            ~ "FL-TP TPA",
      tp_type == "fastga-native"     ~ "FL-TP 1aln",
      cm == "edit-distance"          ~ "EB-TP TPA",
      cm == "diagonal-distance"      ~ "DB-TP TPA"
    ), levels = method_levels),
    # reconstruction time = the stored-format -> alignment step. For cigzip that is
    # decompress + decode; for FL-TP FASTGA decompress is 0 (ALNtoPAF goes straight from .1aln).
    recon_time = decompress_runtime_sec + decode_runtime_sec
  ) %>%
  mk_labels() %>%
  arrange(method, mc)

# FL-TP FASTGA has a single operating point (l=100), so it shows as a point with no connecting line.
p <- ggplot(df, aes(x = size_tpa_bytes, y = recon_time, color = method, group = method)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_point(aes(shape = method), size = 1.8) +
  ggrepel::geom_text_repel(aes(label = mc), size = 2.3, show.legend = FALSE,
                           min.segment.length = 0, segment.size = 0.2, segment.alpha = 0.4,
                           box.padding = 0.15, point.padding = 0.1, max.overlaps = Inf) +
  scale_shape_manual(values = c("FL-TP TPA" = 16, "EB-TP TPA" = 16, "DB-TP TPA" = 16),
                     name = "Method", drop = FALSE) +
  ggh4x::facet_grid2(error_label ~ length_label, scales = "free", independent = "all") +
  scale_x_log10(breaks = scales::breaks_log(n = 6), labels = label_bytes()) +
  scale_y_log10(breaks = scales::breaks_log(n = 7)) +
  annotation_logticks(sides = "bl", colour = "grey55", alpha = 0.6) +
  scale_color_manual(values = method_colors, name = "Method", drop = FALSE) +
  labs(
    x = "File size",
    y = "Reconstruction time (sec)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "figS8_space_time_tradeoff.png"), p, width = 11, height = 12, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS8_space_time_tradeoff.pdf"), p, width = 11, height = 12, bg = "white")
message("Figure saved: ", file.path(fig_dir, "figS8_space_time_tradeoff.{png,pdf}"))
