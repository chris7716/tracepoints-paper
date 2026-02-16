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
        l == 1000 ~ "1 Kb",
        l == 10000 ~ "10 Kb",
        l == 100000 ~ "100 Kb"
      ),
      levels = c("100 bp", "1 Kb", "10 Kb", "100 Kb")
    ),
    error_label = factor(paste0(e * 100, "%"), levels = c("1%", "10%", "20%")),
    delta_label = factor(paste0("t=", mc), levels = c("t=32", "t=64", "t=128", "t=256", "t=512", "t=1024"))
  )

p_delta <- ggplot(df_delta, aes(x = error_label, y = ratio_tpa, fill = delta_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.8) +
  facet_grid2(method ~ length_label, scales = "free_y", independent = "y") +
  scale_fill_manual(values = c("t=32" = "#4292c6", "t=64" = "#2171b5", "t=128" = "#08519c",
                               "t=256" = "#08306b", "t=512" = "#6a51a3", "t=1024" = "#4a1486"),
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
        l == 1000 ~ "1 Kb",
        l == 10000 ~ "10 Kb",
        l == 100000 ~ "100 Kb"
      ),
      levels = c("100 bp", "1 Kb", "10 Kb", "100 Kb")
    ),
    error_label = factor(paste0(e * 100, "%"), levels = c("1%", "10%", "20%"))
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

p_scores <- ggplot(df_scores, aes(x = neg_score, y = frequency)) +
  geom_col(width = 0.10) +
  facet_wrap(~ dataset, nrow = 2, scales = "free_y") +
  scale_x_log10(breaks = 10^(1:7), labels = scales::comma) +
  scale_y_continuous(labels = compact_label, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "|Alignment score|", y = "Frequency") +
  common_theme

ggsave(file.path(fig_dir, "figS3_score_distributions.png"), p_scores, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS3_score_distributions.pdf"), p_scores, width = 10, height = 6, bg = "white")

message("Figure S3 saved: figS3_score_distributions.pdf")

message("\nAll supplementary figures (S1, S2, S3) generated successfully!")
