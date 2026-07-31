#!/usr/bin/env Rscript
# Plot ES-2bit benchmark results from the all.tsv produced by es2bit_bench.
#
# Usage:
#   Rscript plot-es2bit-bench.R <all.tsv> <output_dir>

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: plot-es2bit-bench.R <all.tsv> <output_dir>")
tsv_path <- args[1]
out_dir  <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read and reshape --------------------------------------------------------
raw <- read_tsv(tsv_path, show_col_types = FALSE)

raw <- raw %>%
  mutate(
    n   = as.integer(str_extract(basename(file), "(?<=set_)\\d+")),
    eps_code = str_extract(basename(file), "(?<=_)\\d+(?=\\.paf)"),
    eps = as.numeric(paste0("0.", str_sub(eps_code, 2)))
  )

methods <- c("cigar", "es2bit", "fltp", "fltp_nd", "ebtp", "dbtp")
labels  <- c("CIGAR TEXT", "ES-2bit", "FL-TP", "FL-TP nd", "EB-TP", "DB-TP")

pivot_metric <- function(df, suffix, col_name) {
  df %>%
    select(n, eps, mean_e, ends_with(suffix)) %>%
    pivot_longer(cols = ends_with(suffix), names_to = "method", values_to = col_name) %>%
    mutate(method = str_remove(method, paste0("_", suffix)))
}

df <- pivot_metric(raw, "disk_bits_per_e", "disk_bits_per_e") %>%
  left_join(pivot_metric(raw, "disk_B_per_e",    "disk_B_per_e"),    by = c("n","eps","mean_e","method")) %>%
  left_join(pivot_metric(raw, "gz_B_per_e",      "gz_B_per_e"),      by = c("n","eps","mean_e","method")) %>%
  left_join(pivot_metric(raw, "enc_ns_per_e",    "enc_ns_per_e"),    by = c("n","eps","mean_e","method")) %>%
  left_join(pivot_metric(raw, "dec_ns_per_e",    "dec_ns_per_e"),    by = c("n","eps","mean_e","method")) %>%
  left_join(pivot_metric(raw, "total_dec_ms",    "total_dec_ms"),    by = c("n","eps","mean_e","method")) %>%
  mutate(
    gz_bits_per_e = gz_B_per_e * 8,
    method = factor(method, levels = methods, labels = labels),
    eps_label = paste0(eps * 100, "%"),
    eps_label = fct_reorder(eps_label, eps),
    n_label = factor(
      case_when(
        n == 100 ~ "100 bp",
        n == 1000 ~ "1 Kbp",
        n == 10000 ~ "10 Kbp",
        n == 100000 ~ "100 Kbp"
      ),
      levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp")
    )
  )

# --- Colors and theme (matching the paper) -----------------------------------
method_colors <- c(
  "CIGAR TEXT" = "#1b9e77",  # teal (same as CG BGZIP in paper)
  "ES-2bit"    = "#984ea3",  # dark magenta
  "FL-TP"    = "#7570b3",  # purple
  "FL-TP nd" = "#9e9ac8",  # light purple
  "EB-TP"    = "#d95f02",  # orange
  "DB-TP"    = "#377eb8"   # blue
)

common_theme <- theme_bw(base_size = 15) +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 13),
    axis.text         = element_text(size = 13),
    axis.title        = element_text(size = 15, face = "bold"),
    strip.text        = element_text(size = 15, face = "bold"),
    panel.grid.minor  = element_blank(),
    legend.key.size   = unit(0.5, "cm"),
    legend.margin     = margin(t = -5, unit = "pt")
  )

save_plot <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("  saved", name, "\n")
}

# Shared bar-plot skeleton. Uses linear y-axis with free_y per facet.
bar_plot <- function(df, y_col, y_lab, y_labels = waiver()) {
  ggplot(df, aes(x = eps_label, y = .data[[y_col]], fill = method)) +
    geom_col(position = position_dodge(0.85), width = 0.8) +
    facet_wrap(~ n_label, nrow = 2, scales = "free_y") +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 8), labels = y_labels) +
    labs(x = "Error rate", y = y_lab, fill = "Method") +
    common_theme +
    guides(fill = guide_legend(nrow = 1))
}

# --- Combined: disk (solid) + gzip (alpha) in one 2x2 grid per unit ----------

# 8 bars per error rate: 4 methods × 2 storage types (disk=solid, gzip=faded).
combo_plot <- function(df, disk_col, gz_col, y_lab) {
  long <- df %>%
    select(n_label, eps, eps_label, method, disk = !!sym(disk_col), gzip = !!sym(gz_col)) %>%
    pivot_longer(c(disk, gzip), names_to = "storage", values_to = "value") %>%
    mutate(storage = factor(storage, levels = c("disk", "gzip")))

  # Map method + storage to a combined key for fill and alpha.
  method_lvls <- levels(df$method)
  combo_lvls <- paste(rep(method_lvls, each = 2), c("disk", "gzip"), sep = ".")
  long$combo <- factor(paste(long$method, long$storage, sep = "."), levels = combo_lvls)

  fill_vals  <- setNames(rep(method_colors[method_lvls], each = 2), combo_lvls)
  alpha_vals <- setNames(rep(c(1.0, 0.4), length(method_lvls)), combo_lvls)

  ggplot(long, aes(x = eps_label, y = value, fill = combo, alpha = combo)) +
    geom_col(position = position_dodge(0.9), width = 0.85) +
    facet_wrap(~ n_label, nrow = 2, scales = "free_y") +
    scale_fill_manual(values = fill_vals,
                      breaks = paste(method_lvls, "disk", sep = "."),
                      labels = method_lvls) +
    scale_alpha_manual(values = alpha_vals, guide = "none") +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 8)) +
    labs(x = "Error rate", y = y_lab, fill = "Method",
         caption = "Solid = uncompressed, faded = gzip") +
    common_theme +
    guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)))
}

# Bits per edit uses a BROKEN y-axis
bits_long <- df %>%
  select(n_label, eps, eps_label, method,
         disk = disk_bits_per_e, gzip = gz_bits_per_e) %>%
  pivot_longer(c(disk, gzip), names_to = "storage", values_to = "value") %>%
  mutate(storage = factor(storage, levels = c("disk", "gzip")))
method_lvls <- levels(df$method)
combo_lvls  <- paste(rep(method_lvls, each = 2), c("disk", "gzip"), sep = ".")
bits_long$combo <- factor(paste(bits_long$method, bits_long$storage, sep = "."), levels = combo_lvls)
fill_vals  <- setNames(rep(method_colors[method_lvls], each = 2), combo_lvls)
alpha_vals <- setNames(rep(c(1.0, 0.4), length(method_lvls)), combo_lvls)

low_max <- 3; high_max <- 170
low_vis_end <- 0.55; gap <- 0.03; high_vis_start <- low_vis_end + gap
low_scale  <- low_vis_end / low_max
high_scale <- (1 - high_vis_start) / (high_max - low_max)
tf <- function(y) ifelse(y <= low_max, y * low_scale, high_vis_start + (y - low_max) * high_scale)
bits_long <- bits_long %>% mutate(value_plot = tf(pmin(value, high_max)))

low_ticks  <- c(0, 1, 2, 3)
high_ticks <- c(25, 50, 75, 100, 125, 150)
y_breaks   <- c(tf(low_ticks), tf(high_ticks))
y_labels   <- c(as.character(low_ticks), as.character(high_ticks))

p_bits <- ggplot(bits_long, aes(x = eps_label, y = value_plot, fill = combo, alpha = combo)) +
  geom_col(position = position_dodge(0.9), width = 0.85) +
  geom_hline(yintercept = tf(2), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  annotate("text", x = Inf, y = tf(2), label = "2 bits", hjust = 1.1, vjust = -0.5,
           size = 3.5, color = "grey40") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = low_vis_end, ymax = high_vis_start, fill = "white") +
  geom_hline(yintercept = c(low_vis_end, high_vis_start), linewidth = 0.4) +
  facet_wrap(~ n_label, nrow = 2) +
  scale_fill_manual(values = fill_vals,
                    breaks = paste(method_lvls, "disk", sep = "."), labels = method_lvls) +
  scale_alpha_manual(values = alpha_vals, guide = "none") +
  scale_y_continuous(breaks = y_breaks, labels = y_labels, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Error rate", y = "Bits per edit operation", fill = "Method",
       caption = "Solid = uncompressed, faded = gzip") +
  common_theme +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)))
save_plot(p_bits, "es2bit_bits_per_edit", w = 10, h = 8)
save_plot(combo_plot(df, "disk_B_per_e", "gz_B_per_e", "Bytes per edit operation"),
          "es2bit_bytes_per_edit", w = 10, h = 8)

# --- Timing plots (exclude CIGAR, which is always 0) -------------------------
df_no_cigar <- df %>% filter(method != "CIGAR TEXT")
fmt_k <- scales::label_number(scale_cut = scales::cut_short_scale())
save_plot(bar_plot(df_no_cigar, "dec_ns_per_e",  "Decode ns per edit operation",      y_labels = fmt_k), "es2bit_decode_ns_per_edit", w = 10, h = 8)
save_plot(bar_plot(df_no_cigar, "enc_ns_per_e",  "Encode ns per edit operation",      y_labels = fmt_k), "es2bit_encode_ns_per_edit", w = 10, h = 8)
save_plot(bar_plot(df_no_cigar, "total_dec_ms",  "Total decode time (ms)",  y_labels = fmt_k), "es2bit_total_decode_ms", w = 10, h = 8)

cat("Done.\n")
