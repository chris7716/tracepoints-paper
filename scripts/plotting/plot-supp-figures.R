#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggh4x)

# Configuration: set directories via environment variables or use default relative paths
data_dir <- Sys.getenv("TRACEPOINTS_DATA_DIR", unset = "data")
fig_dir <- Sys.getenv("TRACEPOINTS_FIG_DIR", unset = "paper/figures")

# Read data
data <- read_tsv(file.path(data_dir, "simulated-data", "benchmark.results.tsv"), show_col_types = FALSE)

# Common theme (matching compression_ratios figure)
common_theme <- theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 15, face = "bold"),
    strip.text = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(t = -5, unit = "pt")
  )

# =============================================================================
# Figure S1: Effect of δ parameter on compression ratio
# =============================================================================

df_delta <- data %>%
  filter(memory_mode == "high") %>%
  filter(tp_type == "standard") %>%
  filter(mc %in% c(32, 64, 128, 256, 512, 1024)) %>%
  filter(e != 0.05) %>%
  mutate(
    ratio_tpa = size_tpa_bytes / size_cigar_bytes,
    method = factor(ifelse(cm == "edit-distance", "EB-TP", "DB-TP"), levels = c("EB-TP", "DB-TP")),
    length_label = factor(
      case_when(
        l == 100 ~ "100 bp",
        l == 1000 ~ "1 Kbp",
        l == 10000 ~ "10 Kbp",
        l == 100000 ~ "100 Kbp"
      ),
      levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp")
    ),
    error_label = factor(paste0(e * 100, "%"), levels = c("0.1%", "1%", "10%", "20%")),
    delta_label = factor(paste0("t=", mc), levels = c("t=32", "t=64", "t=128", "t=256", "t=512", "t=1024"))
  )

p_delta <- ggplot(df_delta, aes(x = error_label, y = ratio_tpa, fill = delta_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.78) +
  facet_grid2(method ~ length_label, scales = "free_y", independent = "y") +
  scale_fill_manual(values = c("t=32" = "#9ecae1", "t=64" = "#4292c6", "t=128" = "#08519c",
                               "t=256" = "#253494", "t=512" = "#6a017a", "t=1024" = "#2d004b"),
                    name = "Threshold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Error rate",
    y = "Compression ratio"
  ) +
  common_theme +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(fig_dir, "figS1_delta_effect.png"), p_delta, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS1_delta_effect.pdf"), p_delta, width = 10, height = 6, bg = "white")

message("Figure S1 saved: figS1_delta_effect.png")

# =============================================================================
# Figure S2: Number of tracepoints comparison
# =============================================================================

df_tp_counts <- data %>%
  filter(memory_mode == "high") %>%
  filter((tp_type == "fastga" & mc == 100) | (tp_type == "standard" & mc == 32)) %>%
  filter(e != 0.05) %>%
  mutate(
    method = factor(
      case_when(
        tp_type == "fastga" ~ "FL-TP",
        cm == "edit-distance" ~ "EB-TP",
        cm == "diagonal-distance" ~ "DB-TP"
      ),
      levels = c("FL-TP", "EB-TP", "DB-TP")
    ),
    length_label = factor(
      case_when(
        l == 100 ~ "100 bp",
        l == 1000 ~ "1 Kbp",
        l == 10000 ~ "10 Kbp",
        l == 100000 ~ "100 Kbp"
      ),
      levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp")
    ),
    error_label = factor(paste0(e * 100, "%"), levels = c("0.1%", "1%", "10%", "20%"))
  )

p_tp_counts <- ggplot(df_tp_counts, aes(x = error_label, y = num_tracepoints, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ length_label, nrow = 2, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = c("FL-TP" = "#7570b3", "EB-TP" = "#d95f02", "DB-TP" = "#377eb8"), name = "Method") +
  labs(
    x = "Error rate",
    y = "Number of tracepoints"
  ) +
  common_theme +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(fig_dir, "figS2_tracepoint_counts.png"), p_tp_counts, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS2_tracepoint_counts.pdf"), p_tp_counts, width = 10, height = 6, bg = "white")

message("Figure S2 saved: figS2_tracepoint_counts.pdf")

# =============================================================================
# Figure S3: Alignment score distributions (HPRCv2 and T2T-Primates)
# =============================================================================

score_dir <- file.path(data_dir, "real-data")

df_hprc <- read_tsv(file.path(score_dir, "hprcv2-25k.input-score-frequency.tsv"),
                    show_col_types = FALSE) %>%
  mutate(dataset = "HPRCv2")

df_t2t <- read_tsv(file.path(score_dir, "t2t-ape-pangenome.input-score-frequency.by-target.tsv"),
                   show_col_types = FALSE) %>%
  mutate(dataset = "T2T-Primates")

# Map targets to species (grouped as in plot-supp-figure-primates.R)
species_levels <- c("Human", "Chimpanzee", "Bonobo", "Gorilla", "B. Orangutan", "S. Orangutan", "Siamang")

df_hprc_plot <- df_hprc %>%
  mutate(neg_score = -score, species = NA_character_) %>%
  filter(neg_score > 0)

df_t2t_plot <- df_t2t %>%
  mutate(
    neg_score = -score,
    species = case_when(
      target %in% c("chm13", "grch38", "hg002") ~ "Human",
      target == "mPanTro3"  ~ "Chimpanzee",
      target == "mPanPan1"  ~ "Bonobo",
      target == "mGorGor1"  ~ "Gorilla",
      target == "mPonAbe1"  ~ "B. Orangutan",
      target == "mPonPyg2"  ~ "S. Orangutan",
      target == "mSymSyn1"  ~ "Siamang"
    )
  ) %>%
  filter(neg_score > 0) %>%
  group_by(neg_score, dataset, species) %>%
  summarise(frequency = sum(frequency), .groups = "drop")

df_scores <- bind_rows(df_hprc_plot, df_t2t_plot) %>%
  mutate(
    dataset = factor(dataset, levels = c("HPRCv2", "T2T-Primates")),
    species = factor(species, levels = species_levels)
  )

# Log-spaced bins shared across both datasets (kills integer-score spikes and
# makes bar heights reflect the total alignment count)
log_bin_edges <- 10 ^ seq(
  log10(max(1, min(df_scores$neg_score))),
  log10(max(df_scores$neg_score)),
  length.out = 81  # 80 bins
)
df_scores_binned <- df_scores %>%
  mutate(bin_idx = findInterval(neg_score, log_bin_edges, all.inside = TRUE)) %>%
  group_by(dataset, bin_idx) %>%
  summarise(frequency = sum(frequency), .groups = "drop") %>%
  mutate(
    bin_low  = log_bin_edges[bin_idx],
    bin_high = log_bin_edges[bin_idx + 1],
    bin_center = 10 ^ ((log10(bin_low) + log10(bin_high)) / 2)
  )

species_colors <- c(
  "Human"        = "#E41A1C",
  "Chimpanzee"   = "#377EB8",
  "Bonobo"       = "#4DAF4A",
  "Gorilla"      = "#984EA3",
  "B. Orangutan"  = "#FF7F00",
  "S. Orangutan"  = "#FFFF33",
  "Siamang"      = "#A65628"
)

# Compact y-axis labels: 100k, 200k, etc.
compact_label <- function(x) {
  ifelse(x >= 1e6, paste0(x / 1e6, "M"),
  ifelse(x >= 1e3, paste0(x / 1e3, "k"),
  as.character(x)))
}

p_scores <- ggplot(df_scores_binned,
                   aes(xmin = bin_low, xmax = bin_high, ymin = 0, ymax = frequency)) +
  geom_rect() +
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") +
  scale_x_log10(breaks = 10^(1:7), labels = scales::comma) +
  scale_y_continuous(labels = compact_label, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "|Alignment score|", y = "Frequency") +
  common_theme

ggsave(file.path(fig_dir, "figS3_score_distributions.png"), p_scores, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS3_score_distributions.pdf"), p_scores, width = 10, height = 6, bg = "white")

message("Figure S3 saved: figS3_score_distributions.pdf")

# =============================================================================
# Figure S7: Absolute file sizes (simulated data) — complements main Fig 2
# =============================================================================
library(scales)

df_sizes <- data %>%
  filter(memory_mode == "high") %>%
  filter((tp_type == "fastga" & mc == 100) | (tp_type == "standard" & mc == 32)) %>%
  filter(e != 0.05) %>%
  mutate(
    bgzip_bytes = size_cigar_bgzip_bytes,
    tpa_bytes   = size_tpa_bytes,
    method = case_when(tp_type == "fastga" ~ "FL-TP",
                       cm == "edit-distance" ~ "EB-TP",
                       cm == "diagonal-distance" ~ "DB-TP"),
    length_label = factor(
      case_when(l == 100 ~ "100 bp", l == 1000 ~ "1 Kbp",
                l == 10000 ~ "10 Kbp", l == 100000 ~ "100 Kbp"),
      levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp")),
    error_label = factor(paste0(e * 100, "%"), levels = c("0.1%", "1%", "10%", "20%"))
  )

fltp_s <- df_sizes %>% filter(method == "FL-TP") %>%
  select(l, e, length_label, error_label, bgzip_bytes, fltp_tpa = tpa_bytes)
ebtp_s <- df_sizes %>% filter(method == "EB-TP") %>%
  select(l, e, length_label, error_label, ebtp_tpa = tpa_bytes)
dbtp_s <- df_sizes %>% filter(method == "DB-TP") %>%
  select(l, e, length_label, error_label, dbtp_tpa = tpa_bytes)

sizes_wide <- fltp_s %>%
  left_join(ebtp_s, by = c("l", "e", "length_label", "error_label")) %>%
  left_join(dbtp_s, by = c("l", "e", "length_label", "error_label"))

sizes_long <- sizes_wide %>%
  pivot_longer(cols = c(bgzip_bytes, fltp_tpa, ebtp_tpa, dbtp_tpa),
               names_to = "format", values_to = "bytes") %>%
  mutate(format_label = factor(
    case_when(format == "bgzip_bytes" ~ "CG BGZIP",
              format == "fltp_tpa" ~ "FL-TP TPA",
              format == "ebtp_tpa" ~ "EB-TP TPA",
              format == "dbtp_tpa" ~ "DB-TP TPA"),
    levels = c("CG BGZIP", "FL-TP TPA", "EB-TP TPA", "DB-TP TPA")))

size_colors <- scale_fill_manual(
  values = c("CG BGZIP" = "#1b9e77", "FL-TP TPA" = "#7570b3",
             "EB-TP TPA" = "#d95f02", "DB-TP TPA" = "#377eb8"),
  name = "Method"
)

# Row 1 (100 bp, 1 Kb): linear free_y, no break needed
p_sizes_row1 <- ggplot(
    filter(sizes_long, length_label %in% c("100 bp", "1 Kbp")),
    aes(x = error_label, y = bytes, fill = format_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ length_label, nrow = 1, scales = "free_y") +
  size_colors +
  scale_y_continuous(labels = label_bytes(units = "auto_si"), n.breaks = 8,
                     expand = expansion(mult = c(0, 0.1)), minor_breaks = NULL) +
  labs(x = NULL, y = "File size") +
  common_theme +
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Row 2 (10 Kb, 100 Kb): broken y-axis, per-panel break value
make_broken_size_panel <- function(df_panel, low_max, high_max,
                                   low_visual_end = 0.6, gap_visual = 0.03) {
  high_visual_start <- low_visual_end + gap_visual
  high_visual_end <- 1.0
  low_scale <- low_visual_end / low_max
  high_scale <- (high_visual_end - high_visual_start) / (high_max - low_max)
  tf <- function(y) ifelse(y <= low_max, y * low_scale,
                           high_visual_start + (y - low_max) * high_scale)

  df_plot <- df_panel %>% mutate(bytes_plot = tf(pmin(bytes, high_max)))

  low_ticks  <- pretty(c(0, low_max), n = 4); low_ticks <- low_ticks[low_ticks <= low_max]
  high_ticks <- pretty(c(low_max, high_max), n = 3); high_ticks <- high_ticks[high_ticks > low_max & high_ticks <= high_max]
  brks   <- c(tf(low_ticks), tf(high_ticks))
  labs_y <- label_bytes(units = "auto_si")(c(low_ticks, high_ticks))

  ggplot(df_plot, aes(x = error_label, y = bytes_plot, fill = format_label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = low_visual_end, ymax = high_visual_start, fill = "white") +
    geom_hline(yintercept = low_visual_end,   linewidth = 0.4) +
    geom_hline(yintercept = high_visual_start, linewidth = 0.4) +
    facet_wrap(~ length_label, nrow = 1) +
    size_colors +
    scale_y_continuous(breaks = brks, labels = labs_y,
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "File size") +
    common_theme
}

# Peak values (bytes): 10 Kb ≈ 16 MB; 100 Kb ≈ 150 MB
p_sizes_10k  <- make_broken_size_panel(filter(sizes_long, length_label == "10 Kbp"),
                                       low_max = 2.0e6,  high_max = 1.7e7)
p_sizes_100k <- make_broken_size_panel(filter(sizes_long, length_label == "100 Kbp"),
                                       low_max = 1.5e7, high_max = 1.7e8) +
  theme(axis.title.y = element_blank())

library(patchwork)
p_sizes_row2 <- (p_sizes_10k + p_sizes_100k) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_sizes_row2 <- p_sizes_row2 & labs(x = "Error rate")

p_sizes <- p_sizes_row1 / p_sizes_row2 +
  plot_layout(heights = c(0.45, 0.55), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "figS7_absolute_sizes.png"), p_sizes, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS7_absolute_sizes.pdf"), p_sizes, width = 10, height = 6, bg = "white")
message("Figure S7 saved: figS7_absolute_sizes.pdf")

message("\nAll supplementary figures (S1, S2, S3, S7) generated successfully!")
