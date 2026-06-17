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
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
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

# Parameter sweeps: EB-TP/DB-TP vary the threshold (delta, b); FL-TP varies the
# trace spacing l. The value sets differ, so the two are given SEPARATE fill
# scales (and legends) via ggnewscale
library(ggnewscale)

delta_levels  <- c(32, 64, 128, 256, 512, 1024)     # EB-TP (delta), DB-TP (b)
fastga_levels <- c(100, 200, 300, 500, 1000, 2000)  # FL-TP (trace spacing l)
threshold_cols <- setNames(c("#9ecae1", "#4292c6", "#08519c",
                             "#253494", "#6a017a", "#2d004b"), as.character(delta_levels))
length_cols    <- setNames(c("#fdd0a2", "#fdae6b", "#fd8d3c",
                             "#f16913", "#d94801", "#7f2704"), as.character(fastga_levels))

df_delta <- data %>%
  filter(memory_mode == "high") %>%
  filter((tp_type == "standard" & mc %in% delta_levels) |
         (tp_type == "fastga"   & mc %in% fastga_levels)) %>%
  filter(e %in% c(0.001, 0.01, 0.05, 0.10, 0.20)) %>%
  mutate(
    ratio_tpa = size_tpa_bytes / size_cigar_bytes,
    method = factor(
      case_when(
        tp_type == "fastga" ~ "FL-TP TPA",
        cm == "edit-distance" ~ "EB-TP TPA",
        cm == "diagonal-distance" ~ "DB-TP TPA"
      ),
      levels = c("FL-TP TPA", "EB-TP TPA", "DB-TP TPA")
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
    error_label = factor(paste0(e * 100, "%"), levels = c("0.1%", "1%", "5%", "10%", "20%")),
    threshold_label = factor(as.character(mc), levels = as.character(delta_levels)),
    spacing_label   = factor(as.character(mc), levels = as.character(fastga_levels))
  )

df_fltp <- df_delta %>% filter(method == "FL-TP TPA")
df_ebdb <- df_delta %>% filter(method != "FL-TP TPA")

p_delta <- ggplot() +
  # FL-TP: trace spacing l -> "Length" legend (orange palette)
  geom_bar(data = df_fltp, aes(x = error_label, y = ratio_tpa, fill = spacing_label),
           stat = "identity", position = position_dodge(width = 0.9), width = 0.78) +
  scale_fill_manual(values = length_cols, name = "Length (l)",
                    guide = guide_legend(order = 1, nrow = 1)) +
  new_scale_fill() +
  # EB-TP / DB-TP: threshold delta, b -> "Threshold" legend (blue-purple palette)
  geom_bar(data = df_ebdb, aes(x = error_label, y = ratio_tpa, fill = threshold_label),
           stat = "identity", position = position_dodge(width = 0.9), width = 0.78) +
  scale_fill_manual(values = threshold_cols, name = "Threshold (t)",
                    guide = guide_legend(order = 2, nrow = 1)) +
  facet_grid2(method ~ length_label, scales = "free_y", independent = "y") +
  # Denser y ticks: ~0.1 on the large-range panels, finer on the small ones
  # (free_y panels range from ~0.08 to ~0.7, so a fixed 0.1 step would leave the
  # smallest panels with a single tick).
  scale_y_continuous(n.breaks = 8, expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Error rate",
    y = "Compression ratio"
  ) +
  common_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.box = "vertical", legend.spacing.y = unit(12, "pt"))

ggsave(file.path(fig_dir, "figS1_delta_effect.png"), p_delta, width = 10, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS1_delta_effect.pdf"), p_delta, width = 10, height = 7.5, bg = "white")

message("Figure S1 saved: figS1_delta_effect.png")

# =============================================================================
# Figure S2: Number of tracepoints comparison
# =============================================================================

df_tp_counts <- data %>%
  filter(memory_mode == "high") %>%
  filter((tp_type == "fastga" & mc == 100) | (tp_type == "standard" & mc == 32)) %>%
  filter(e %in% c(0.001, 0.01, 0.05, 0.10, 0.20)) %>%
  mutate(
    method = factor(
      case_when(
        tp_type == "fastga" ~ "FL-TP TPA",
        cm == "edit-distance" ~ "EB-TP TPA",
        cm == "diagonal-distance" ~ "DB-TP TPA"
      ),
      levels = c("FL-TP TPA", "EB-TP TPA", "DB-TP TPA")
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
    error_label = factor(paste0(e * 100, "%"), levels = c("0.1%", "1%", "5%", "10%", "20%"))
  )

p_tp_counts <- ggplot(df_tp_counts, aes(x = error_label, y = num_tracepoints, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ length_label, nrow = 2, scales = "free_y") +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = c("FL-TP TPA" = "#7570b3", "EB-TP TPA" = "#d95f02", "DB-TP TPA" = "#377eb8"), name = "Method") +
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
# Figure S3: Alignment score distributions, raw and length-normalized, by target
# =============================================================================

score_dir <- file.path(data_dir, "real-data")
RAW_LAB <- "|Alignment score|"
PB_LAB  <- "|Alignment score / alignment length|"

# Map PanSN target prefix to species. HPRCv2 is all human-vs-human, coloured "Human".
tgt_species <- function(t) case_when(
  t %in% c("chm13", "grch38", "hg002") ~ "Human",
  t == "mPanTro3" ~ "Chimpanzee",
  t == "mPanPan1" ~ "Bonobo",
  t == "mGorGor1" ~ "Gorilla",
  t == "mPonAbe1" ~ "B. Orangutan",
  t == "mPonPyg2" ~ "S. Orangutan",
  t == "mSymSyn1" ~ "Siamang",
  TRUE ~ NA_character_
)
species_levels <- c("Human", "Chimpanzee", "Bonobo", "Gorilla", "B. Orangutan", "S. Orangutan", "Siamang")
species_colors <- c("Human" = "#E41A1C", "Chimpanzee" = "#377EB8", "Bonobo" = "#4DAF4A",
                    "Gorilla" = "#984EA3", "B. Orangutan" = "#FF7F00",
                    "S. Orangutan" = "#FFD92F", "Siamang" = "#A65628")

# Four histograms: raw |score| and length-normalized |score|/length (per aligned
# base), for HPRCv2 (pooled, all "Human") and T2T apes (by target genome).
sc_hprc_raw <- read_tsv(file.path(score_dir, "hprcv2-25k.input-score-frequency.tsv"), show_col_types = FALSE) %>%
  transmute(neg = -score, frequency, dataset = "HPRCv2", metric = RAW_LAB, species = "Human")
sc_hprc_pb  <- read_tsv(file.path(score_dir, "hprcv2-25k.input-scorePerBase-frequency.tsv"), show_col_types = FALSE) %>%
  transmute(neg = -score_per_base, frequency, dataset = "HPRCv2", metric = PB_LAB, species = "Human")
sc_t2t_raw  <- read_tsv(file.path(score_dir, "t2t-ape-pangenome.input-score-frequency.by-target.tsv"), show_col_types = FALSE) %>%
  transmute(neg = -score, frequency, dataset = "T2T-Primates", metric = RAW_LAB, species = tgt_species(target))
sc_t2t_pb   <- read_tsv(file.path(score_dir, "t2t-ape-pangenome.input-scorePerBase-frequency.by-target.tsv"), show_col_types = FALSE) %>%
  transmute(neg = -score_per_base, frequency, dataset = "T2T-Primates", metric = PB_LAB, species = tgt_species(target))

df_scores <- bind_rows(sc_hprc_raw, sc_hprc_pb, sc_t2t_raw, sc_t2t_pb) %>% filter(neg > 0)

# Log-spaced bins per metric (shared across datasets so x is comparable within each
# column); then stack species within each bin.
df_scores_binned <- df_scores %>%
  group_by(metric) %>%
  group_modify(~{
    edges <- 10 ^ seq(log10(min(.x$neg)), log10(max(.x$neg)), length.out = 81)
    idx <- findInterval(.x$neg, edges, all.inside = TRUE)
    mutate(.x, bin_low = edges[idx], bin_high = edges[idx + 1])
  }) %>%
  ungroup() %>%
  group_by(metric, dataset, bin_low, bin_high, species) %>%
  summarise(frequency = sum(frequency), .groups = "drop") %>%
  mutate(species = factor(species, levels = species_levels)) %>%
  arrange(metric, dataset, bin_low, species) %>%
  group_by(metric, dataset, bin_low, bin_high) %>%
  mutate(ymax = cumsum(frequency), ymin = ymax - frequency) %>%
  ungroup() %>%
  mutate(dataset = factor(dataset, levels = c("HPRCv2", "T2T-Primates")),
         metric  = factor(metric,  levels = c(RAW_LAB, PB_LAB)))

p_scores <- ggplot(df_scores_binned,
                   aes(xmin = bin_low, xmax = bin_high, ymin = ymin, ymax = ymax, fill = species)) +
  geom_rect() +
  ggh4x::facet_grid2(dataset ~ metric, scales = "free", independent = "y") +
  scale_x_log10(breaks = scales::breaks_log(n = 6),
                labels = scales::label_number(scale_cut = scales::cut_short_scale(), drop0trailing = TRUE)) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = species_colors, name = "Target genome") +
  labs(x = NULL, y = "Frequency") +
  common_theme +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(fig_dir, "figS3_score_distributions.png"), p_scores, width = 12, height = 7, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS3_score_distributions.pdf"), p_scores, width = 12, height = 7, bg = "white")

message("Figure S3 saved: figS3_score_distributions.pdf")

# =============================================================================
# Figure S7: Absolute file sizes (simulated data) — complements main Fig 2
# =============================================================================
library(scales)

df_sizes <- data %>%
  filter(memory_mode == "high") %>%
  filter((tp_type == "fastga" & mc == 100) |
         (tp_type == "fastga-native" & mc == 100) |
         (tp_type == "standard" & mc == 32)) %>%
  filter(e %in% c(0.001, 0.01, 0.05, 0.10, 0.20)) %>%
  mutate(
    bgzip_bytes = size_cigar_bgzip_bytes,
    tpa_bytes   = size_tpa_bytes,
    method = case_when(tp_type == "fastga" ~ "FL-TP",
                       tp_type == "fastga-native" ~ "FL-TP FASTGA",
                       cm == "edit-distance" ~ "EB-TP",
                       cm == "diagonal-distance" ~ "DB-TP"),
    length_label = factor(
      case_when(l == 100 ~ "100 bp", l == 1000 ~ "1 Kbp",
                l == 10000 ~ "10 Kbp", l == 100000 ~ "100 Kbp"),
      levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp")),
    error_label = factor(paste0(e * 100, "%"), levels = c("0.1%", "1%", "5%", "10%", "20%"))
  )

fltp_s <- df_sizes %>% filter(method == "FL-TP") %>%
  select(l, e, length_label, error_label, bgzip_bytes, fltp_cigzip_tpa = tpa_bytes)
fltp_native_s <- df_sizes %>% filter(method == "FL-TP FASTGA") %>%
  select(l, e, length_label, error_label, fltp_fastga_tpa = tpa_bytes)
ebtp_s <- df_sizes %>% filter(method == "EB-TP") %>%
  select(l, e, length_label, error_label, ebtp_tpa = tpa_bytes)
dbtp_s <- df_sizes %>% filter(method == "DB-TP") %>%
  select(l, e, length_label, error_label, dbtp_tpa = tpa_bytes)

sizes_wide <- fltp_s %>%
  left_join(fltp_native_s, by = c("l", "e", "length_label", "error_label")) %>%
  left_join(ebtp_s, by = c("l", "e", "length_label", "error_label")) %>%
  left_join(dbtp_s, by = c("l", "e", "length_label", "error_label"))

sizes_long <- sizes_wide %>%
  pivot_longer(cols = c(bgzip_bytes, fltp_cigzip_tpa, fltp_fastga_tpa, ebtp_tpa, dbtp_tpa),
               names_to = "format", values_to = "bytes") %>%
  mutate(format_label = factor(
    case_when(format == "bgzip_bytes" ~ "CG BGZIP",
              format == "fltp_cigzip_tpa" ~ "FL-TP TPA",
              format == "fltp_fastga_tpa" ~ "FL-TP 1aln",
              format == "ebtp_tpa" ~ "EB-TP TPA",
              format == "dbtp_tpa" ~ "DB-TP TPA"),
    levels = c("CG BGZIP", "FL-TP TPA", "FL-TP 1aln", "EB-TP TPA", "DB-TP TPA")))

size_colors <- scale_fill_manual(
  values = c("CG BGZIP" = "#1b9e77", "FL-TP TPA" = "#7570b3",
             "FL-TP 1aln" = "#9e9ac8", "EB-TP TPA" = "#d95f02", "DB-TP TPA" = "#377eb8"),
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

# (The former Figure S8, FL-TP compression ratio across trace spacings, was
# dropped: figS1 now shows the same FL-TP spacing sweep as its top method row.)

message("\nAll supplementary figures (S1, S2, S3, S7) generated successfully!")
