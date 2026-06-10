#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)

# Configuration: set directories via environment variables or use default relative paths
data_dir <- Sys.getenv("TRACEPOINTS_DATA_DIR", unset = "data")
fig_dir <- Sys.getenv("TRACEPOINTS_FIG_DIR", unset = "paper/figures")

# Read data
data <- read_tsv(file.path(data_dir, "simulated-data", "benchmark.results.tsv"), show_col_types = FALSE)

# Filter for δ=32 (mc=32) and ultralow memory mode
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
    error_label = factor(
      paste0(e * 100, "%"),
      levels = c("1%", "5%", "10%", "20%")
    )
  )

# Separate FL-TP data
fltp_data <- df %>%
  filter(method == "FL-TP") %>%
  select(l, e, length_label, error_label, ratio_bgzip, fltp_tpa = ratio_tpa)

# Get EB-TP and DB-TP data
ebtp_data <- df %>%
  filter(method == "EB-TP") %>%
  select(l, e, length_label, error_label, ebtp_tpa = ratio_tpa)

dbtp_data <- df %>%
  filter(method == "DB-TP") %>%
  select(l, e, length_label, error_label, dbtp_tpa = ratio_tpa)

# Join all data
plot_data <- fltp_data %>%
  left_join(ebtp_data, by = c("l", "e", "length_label", "error_label")) %>%
  left_join(dbtp_data, by = c("l", "e", "length_label", "error_label"))

# Pivot to long format with all 4 methods as bars
plot_long <- plot_data %>%
  pivot_longer(
    cols = c(ratio_bgzip, fltp_tpa, ebtp_tpa, dbtp_tpa),
    names_to = "format",
    values_to = "ratio"
  ) %>%
  mutate(
    format_label = factor(
      case_when(
        format == "ratio_bgzip" ~ "CG BGZIP",
        format == "fltp_tpa" ~ "FL-TP TPA",
        format == "ebtp_tpa" ~ "EB-TP TPA",
        format == "dbtp_tpa" ~ "DB-TP TPA"
      ),
      levels = c("CG BGZIP", "FL-TP TPA", "EB-TP TPA", "DB-TP TPA")
    )
  )

# Common theme
common_theme <- theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 15, face = "bold"),
    strip.text = element_text(size = 15, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(t = -5, unit = "pt")
  )

color_scale <- scale_fill_manual(
  values = c(
    "CG BGZIP" = "#1b9e77",
    "FL-TP TPA" = "#7570b3",
    "EB-TP TPA" = "#d95f02",
    "DB-TP TPA" = "#377eb8"
  ),
  name = "Method"
)

# =============================================================================
# Figure 1: Compression ratios (3 error rates: 1%, 10%, 20%)
# =============================================================================

library(grid)
library(gridExtra)

plot_long_3rates <- plot_long %>% filter(error_label != "5%")

plot_long_3r_row1 <- plot_long_3rates %>% filter(length_label %in% c("100 bp", "1 Kbp"))
plot_long_3r_row2 <- plot_long_3rates %>% filter(length_label %in% c("10 Kbp", "100 Kbp"))

p_3r_row1 <- ggplot(plot_long_3r_row1, aes(x = error_label, y = ratio, fill = format_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ length_label, nrow = 1) +
  color_scale +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Compression ratio") +
  common_theme +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Broken y-axis for row 2 (10 Kb, 100 Kb): expand the 0-0.05 range visually
# so small tracepoint ratios are legible, while BGZIP/FL-TP bars (up to ~0.3)
# remain visible above the break.
low_max <- 0.05
high_max <- 0.30
low_visual_end <- 0.18
gap_visual <- 0.01
high_visual_start <- low_visual_end + gap_visual
high_visual_end <- 0.30
low_scale <- low_visual_end / low_max
high_scale <- (high_visual_end - high_visual_start) / (high_max - low_max)
transform_y <- function(y) {
  ifelse(y <= low_max, y * low_scale,
         high_visual_start + (y - low_max) * high_scale)
}

plot_long_3r_row2 <- plot_long_3r_row2 %>%
  mutate(ratio_plot = transform_y(pmin(ratio, high_max)))

break_data <- data.frame(length_label = factor(c("10 Kbp", "100 Kbp"),
  levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp")))

low_ticks <- seq(0, low_max, by = 0.01)
high_ticks <- c(0.1, 0.2, 0.3)
breaks_r2 <- c(transform_y(low_ticks), transform_y(high_ticks))
labels_r2 <- c(sprintf("%.2f", low_ticks), sprintf("%.1f", high_ticks))

p_3r_row2 <- ggplot(plot_long_3r_row2, aes(x = error_label, y = ratio_plot, fill = format_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_rect(data = break_data,
            aes(xmin = -Inf, xmax = Inf,
                ymin = low_visual_end, ymax = high_visual_start),
            inherit.aes = FALSE, fill = "white") +
  geom_hline(data = break_data, aes(yintercept = low_visual_end), linewidth = 0.4) +
  geom_hline(data = break_data, aes(yintercept = high_visual_start), linewidth = 0.4) +
  facet_wrap(~ length_label, nrow = 1) +
  color_scale +
  scale_y_continuous(breaks = breaks_r2, labels = labels_r2,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Error rate", y = "Compression ratio") +
  common_theme +
  guides(fill = guide_legend(nrow = 1))

g1_3r <- ggplotGrob(p_3r_row1)
g2_3r <- ggplotGrob(p_3r_row2)
maxWidth_3r <- unit.pmax(g1_3r$widths, g2_3r$widths)
g1_3r$widths <- maxWidth_3r
g2_3r$widths <- maxWidth_3r
p_comp_3rates <- grid.arrange(g1_3r, g2_3r, ncol = 1, heights = c(0.45, 0.55))

ggsave(file.path(fig_dir, "compression_ratios.png"), p_comp_3rates, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "compression_ratios.pdf"), p_comp_3rates, width = 10, height = 6, bg = "white")

message("Figure saved: compression_ratios.png")

# =============================================================================
# Figure 2: Decoding cost (runtime + memory)
# =============================================================================

library(patchwork)
library(cowplot)

df_runtime <- data %>%
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
    error_label = factor(paste0(e * 100, "%"), levels = c("1%", "10%", "20%"))
  )

method_levels <- c("ORIGINAL", "CG BGZIP", "FL-TP", "EB-TP", "DB-TP")

# Prepare runtime data for ORIGINAL (original alignment computation) - use FL-TP rows as reference
df_align_runtime <- df_runtime %>%
  filter(method == "FL-TP") %>%
  mutate(method = factor("ORIGINAL", levels = method_levels)) %>%
  select(l, e, method, length_label, error_label, runtime = align_runtime_sec)

# Prepare runtime data for CG BGZIP (bgzip decompression) - use FL-TP rows as reference
df_bgzip_runtime <- df_runtime %>%
  filter(method == "FL-TP") %>%
  mutate(method = factor("CG BGZIP", levels = method_levels)) %>%
  select(l, e, method, length_label, error_label, runtime = bgzip_decompress_runtime_sec)

# Prepare runtime data: TPA -> PAF = decompress + decode
df_fltp_runtime <- df_runtime %>%
  filter(method == "FL-TP") %>%
  mutate(runtime = decompress_runtime_sec + decode_runtime_sec) %>%
  select(l, e, method, length_label, error_label, runtime)

df_ebdb_runtime <- df_runtime %>%
  filter(method %in% c("EB-TP", "DB-TP")) %>%
  mutate(runtime = decompress_runtime_sec + decode_heuristic_runtime_sec) %>%
  select(l, e, method, length_label, error_label, runtime)

df_combined_runtime <- bind_rows(df_align_runtime, df_bgzip_runtime, df_fltp_runtime, df_ebdb_runtime) %>%
  mutate(method = factor(method, levels = method_levels))

# Prepare memory data for ORIGINAL
df_align_memory <- df_runtime %>%
  filter(method == "FL-TP") %>%
  mutate(method = factor("ORIGINAL", levels = method_levels)) %>%
  select(l, e, method, length_label, error_label, memory_kb = align_memory_kb) %>%
  mutate(memory_mb = memory_kb / 1024)

# Prepare memory data for CG BGZIP
df_bgzip_memory <- df_runtime %>%
  filter(method == "FL-TP") %>%
  mutate(method = factor("CG BGZIP", levels = method_levels)) %>%
  select(l, e, method, length_label, error_label, memory_kb = bgzip_decompress_memory_kb) %>%
  mutate(memory_mb = memory_kb / 1024)

# Prepare memory data: TPA -> PAF = max(decompress, decode) for peak memory
df_fltp_memory <- df_runtime %>%
  filter(method == "FL-TP") %>%
  mutate(memory_kb = pmax(decompress_memory_kb, decode_memory_kb)) %>%
  select(l, e, method, length_label, error_label, memory_kb) %>%
  mutate(memory_mb = memory_kb / 1024)

df_ebdb_memory <- df_runtime %>%
  filter(method %in% c("EB-TP", "DB-TP")) %>%
  mutate(memory_kb = pmax(decompress_memory_kb, decode_heuristic_memory_kb)) %>%
  select(l, e, method, length_label, error_label, memory_kb) %>%
  mutate(memory_mb = memory_kb / 1024)

df_combined_memory <- bind_rows(df_align_memory, df_bgzip_memory, df_fltp_memory, df_ebdb_memory) %>%
  mutate(method = factor(method, levels = method_levels))

method_colors_5 <- c("ORIGINAL" = "#e41a1c", "CG BGZIP" = "#1b9e77",
                      "FL-TP" = "#7570b3", "EB-TP" = "#d95f02", "DB-TP" = "#377eb8")

library(ggh4x)

# --- Broken axis for 100 Kb runtime panel ---
# Transform data: squish the 50-200 gap into a 5-unit visual gap
break_lower <- 35
break_upper <- 210
gap_size <- 5
offset <- break_upper - break_lower - gap_size  # 145

df_runtime_plot <- df_combined_runtime %>%
  mutate(
    runtime_plot = ifelse(length_label == "100 Kbp" & runtime > break_lower,
                          runtime - offset, runtime)
  )

# Break indicator (white band + ~ marks on y-axis) only in 100 Kb facet
break_data <- data.frame(
  length_label = factor("100 Kbp", levels = c("100 bp", "1 Kbp", "10 Kbp", "100 Kbp"))
)


# Custom breaks/labels for 100 Kb y-axis
breaks_100kb <- c(seq(0, 30, by = 10), 35, break_lower + gap_size + seq(0, 10, by = 5))
labels_100kb <- c(as.character(seq(0, 30, by = 10)), "35", as.character(seq(210, 220, by = 5)))

# Runtime panel (top - no x-axis labels/title)
p_dc_runtime <- ggplot(df_runtime_plot, aes(x = error_label, y = runtime_plot, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  # Break indicators: white band + dashed boundary lines (only in 100 Kb facet)
  geom_rect(data = break_data,
            aes(xmin = -Inf, xmax = Inf, ymin = break_lower, ymax = break_lower + gap_size),
            inherit.aes = FALSE, fill = "white") +
  geom_hline(data = break_data, aes(yintercept = break_lower),
             linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(data = break_data, aes(yintercept = break_lower + gap_size),
             linetype = "solid", color = "black", linewidth = 0.5) +
  facet_wrap(~ length_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = method_colors_5, name = "Method") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  facetted_pos_scales(
    y = list(NULL, NULL, NULL,
             scale_y_continuous(breaks = breaks_100kb, labels = labels_100kb,
                                expand = expansion(mult = c(0, 0.05))))
  ) +
  labs(x = NULL, y = "Time (sec)") +
  common_theme +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Memory panel (bottom - with x-axis, unchanged)
p_dc_memory <- ggplot(df_combined_memory, aes(x = error_label, y = memory_mb, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ length_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = method_colors_5, name = "Method") +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.05)), labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = "Error rate", y = "Peak memory (MB)") +
  common_theme +
  theme(strip.text = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

# Convert to grobs and align panel widths
g1_dc <- ggplotGrob(p_dc_runtime)
g2_dc <- ggplotGrob(p_dc_memory + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)))

# Interrupt the y-axis line at the break zone for the 100 Kb panel
# Panel border is rendered on top of geoms, so we overlay white rects on the gtable
idx_100kb <- which(g1_dc$layout$name == "panel-4-1")
t_100kb <- g1_dc$layout$t[idx_100kb]
l_100kb <- g1_dc$layout$l[idx_100kb]

# Convert break zone to NPC (normalized panel coordinates)
max_y_100kb <- max(df_runtime_plot$runtime_plot[df_runtime_plot$length_label == "100 Kbp"])
y_expand <- 0.05
y_range <- max_y_100kb * (1 + y_expand)
y_npc_lo <- break_lower / y_range
y_npc_hi <- (break_lower + gap_size) / y_range

# White rects covering left and right panel borders at the break zone
left_cover <- rectGrob(
  x = unit(0, "npc"), y = unit(y_npc_lo, "npc"),
  width = unit(1.5, "mm"), height = unit(y_npc_hi - y_npc_lo, "npc"),
  just = c("centre", "bottom"), gp = gpar(fill = "white", col = NA)
)
right_cover <- rectGrob(
  x = unit(1, "npc"), y = unit(y_npc_lo, "npc"),
  width = unit(1.5, "mm"), height = unit(y_npc_hi - y_npc_lo, "npc"),
  just = c("centre", "bottom"), gp = gpar(fill = "white", col = NA)
)
g1_dc <- gtable::gtable_add_grob(g1_dc, left_cover, t = t_100kb, l = l_100kb, name = "break-cover-left")
g1_dc <- gtable::gtable_add_grob(g1_dc, right_cover, t = t_100kb, l = l_100kb, name = "break-cover-right")

maxWidth_dc <- unit.pmax(g1_dc$widths, g2_dc$widths)
g1_dc$widths <- maxWidth_dc
g2_dc$widths <- maxWidth_dc
p_decoding_cost <- grid.arrange(g1_dc, g2_dc, ncol = 1, heights = c(0.45, 0.55))

ggsave(file.path(fig_dir, "decoding_cost_broken.png"), p_decoding_cost, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "decoding_cost_broken.pdf"), p_decoding_cost, width = 10, height = 6, bg = "white")

message("Figure saved: decoding_cost_broken.png/pdf")

# --- Log-scale version in milliseconds (Figure 3 in the manuscript) ---
# Reviewer requested log time axes; ms keeps all values >1 so log bars grow upward.
log_ticks_major <- as.vector(outer(c(1, 3), 10^(0:6)))
log_ticks_minor <- as.vector(outer(c(2, 5), 10^(0:6)))

df_runtime_ms <- df_combined_runtime %>% mutate(runtime_ms = runtime * 1000)

p_dc_runtime_log <- ggplot(df_runtime_ms, aes(x = error_label, y = runtime_ms, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ length_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = method_colors_5, name = "Method") +
  scale_y_log10(expand = expansion(mult = c(0, 0.05)),
                breaks = log_ticks_major, minor_breaks = log_ticks_minor,
                labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = NULL, y = "Time (ms)") +
  common_theme +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_dc_memory_log <- ggplot(df_combined_memory, aes(x = error_label, y = memory_mb, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ length_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = method_colors_5, name = "Method") +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.05)),
                breaks = log_ticks_major, minor_breaks = log_ticks_minor,
                labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = "Error rate", y = "Peak memory (MB)") +
  common_theme +
  theme(strip.text = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

g1_dc_log <- ggplotGrob(p_dc_runtime_log)
g2_dc_log <- ggplotGrob(p_dc_memory_log + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1)))
maxWidth_dc_log <- unit.pmax(g1_dc_log$widths, g2_dc_log$widths)
g1_dc_log$widths <- maxWidth_dc_log
g2_dc_log$widths <- maxWidth_dc_log
p_decoding_cost_log <- grid.arrange(g1_dc_log, g2_dc_log, ncol = 1, heights = c(0.45, 0.55))

ggsave(file.path(fig_dir, "decoding_cost.png"), p_decoding_cost_log, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "decoding_cost.pdf"), p_decoding_cost_log, width = 10, height = 6, bg = "white")

message("Figure saved: decoding_cost.png/pdf")
