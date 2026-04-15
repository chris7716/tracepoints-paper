#!/usr/bin/env Rscript

# Supplementary Figure S4: Primate Pangenome Performance
# 2x2 grid showing TPA size, tracepoints, decode time, and peak memory
# across 16 target genomes ordered by phylogeny

library(tidyverse)
library(ggplot2)
library(patchwork)

# Configuration: set directories via environment variables or use default relative paths
data_dir <- Sys.getenv("TRACEPOINTS_DATA_DIR", unset = "data")
fig_dir <- Sys.getenv("TRACEPOINTS_FIG_DIR", unset = "paper/figures")

# Read primate benchmark data
data <- read_tsv(file.path(data_dir, "real-data", "t2t-ape-pangenome.benchmark.results.tsv"), show_col_types = FALSE)

# Process data
df <- data %>%
  mutate(
    # Extract target name
    target = gsub(".p70.aln", "", paf_file),
    # Assign species
    species = case_when(
      grepl("^chm13", target) ~ "Human",
      grepl("^grch38", target) ~ "Human",
      grepl("^hg002", target) ~ "Human",
      grepl("^mPanTro", target) ~ "Chimpanzee",
      grepl("^mPanPan", target) ~ "Bonobo",
      grepl("^mGorGor", target) ~ "Gorilla",
      grepl("^mPonAbe", target) ~ "B. Orangutan",
      grepl("^mPonPyg", target) ~ "S. Orangutan",
      grepl("^mSymSyn", target) ~ "Siamang",
      TRUE ~ "Other"
    ),
    # Method label (includes threshold parameter)
    method = case_when(
      cm == "edit-distance" ~ paste0("EB-TP (\u03b4=", mc, ")"),
      cm == "diagonal-distance" ~ paste0("DB-TP (b=", mc, ")")
    ),
    # Compute metrics
    tpa_size_mb = size_tpa_bytes / 1e6,
    tracepoints_m = num_tracepoints / 1e6,
    decode_time_sec = decode_runtime_sec,
    peak_memory_gb = decode_memory_kb / 1024 / 1024
  )

# Order targets by phylogeny
target_order <- c(
  "chm13#1", "grch38#1", "hg002#M", "hg002#P",
  "mPanTro3#1", "mPanTro3#2", "mPanPan1#M", "mPanPan1#P",
  "mGorGor1#M", "mGorGor1#P",
  "mPonAbe1#1", "mPonAbe1#2", "mPonPyg2#1", "mPonPyg2#2",
  "mSymSyn1#1", "mSymSyn1#2"
)

# Create shorter labels for x-axis
target_labels <- c(
  "chm13#1" = "CHM13", "grch38#1" = "GRCh38", "hg002#M" = "HG002m", "hg002#P" = "HG002p",
  "mPanTro3#1" = "PTro1", "mPanTro3#2" = "PTro2", "mPanPan1#M" = "PPanM", "mPanPan1#P" = "PPanP",
  "mGorGor1#M" = "GorM", "mGorGor1#P" = "GorP",
  "mPonAbe1#1" = "PAbe1", "mPonAbe1#2" = "PAbe2", "mPonPyg2#1" = "PPyg1", "mPonPyg2#2" = "PPyg2",
  "mSymSyn1#1" = "Sym1", "mSymSyn1#2" = "Sym2"
)

df <- df %>%
  mutate(
    target = factor(target, levels = target_order),
    target_label = factor(target_labels[as.character(target)], levels = target_labels),
    method = factor(method, levels = c(
      "EB-TP (\u03b4=32)", "EB-TP (\u03b4=64)", "EB-TP (\u03b4=128)",
      "DB-TP (b=32)"
    )),
    species = factor(species, levels = c("Human", "Chimpanzee", "Bonobo", "Gorilla", "B. Orangutan", "S. Orangutan", "Siamang"))
  )

# Common theme (matching main figures)
common_theme <- theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 15, face = "bold"),
    strip.text = element_text(size = 15, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(t = -5, unit = "pt")
  )

# Color scale: warm gradient for EB-TP, cool gradient for DB-TP
color_scale <- scale_fill_manual(
  values = c(
    "EB-TP (\u03b4=32)" = "#d95f02", "EB-TP (\u03b4=64)" = "#e6950a", "EB-TP (\u03b4=128)" = "#f0c050",
    "DB-TP (b=32)" = "#1b5e9b"
  ),
  name = "Method"
)

# Species group separators (vertical lines between groups)
# Human(4) | Chimp+Bonobo(4) | Gorilla(2) | Orangutan(4) | Siamang(2)
species_separators <- geom_vline(
  xintercept = c(4.5, 8.5, 10.5, 14.5),
  linetype = "dashed", color = "gray50", linewidth = 0.5
)

# Species labels data frame (positions are centers of each group)
species_labels_df <- data.frame(
  x = c(2.5, 6.5, 9.5, 12.5, 15.5),
  label = c("Human", "Chimp./Bonobo", "Gorilla", "Orangutan", "Siamang")
)

# Panel A: TPA size (MB)
p_tpa_size <- ggplot(df, aes(x = target_label, y = tpa_size_mb, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
  species_separators +
  geom_text(data = species_labels_df, aes(x = x, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3.5, fontface = "plain") +
  color_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "TPA size (MB)") +
  common_theme +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel B: Tracepoints (millions)
p_tracepoints <- ggplot(df, aes(x = target_label, y = tracepoints_m, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
  species_separators +
  geom_text(data = species_labels_df, aes(x = x, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3.5, fontface = "plain") +
  color_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of tracepoints (M)") +
  common_theme +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel C: Decode time (sec) - log scale
p_decode_time <- ggplot(df, aes(x = target_label, y = decode_time_sec, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
  species_separators +
  geom_text(data = species_labels_df, aes(x = x, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3.5, fontface = "plain") +
  color_scale +
  scale_y_log10(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Target genome", y = "Decode time (sec)") +
  common_theme +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# Panel D: Peak memory (GB) - log scale
p_peak_memory <- ggplot(df, aes(x = target_label, y = peak_memory_gb, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
  species_separators +
  geom_text(data = species_labels_df, aes(x = x, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3.5, fontface = "plain") +
  color_scale +
  scale_y_log10(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Target genome", y = "Peak memory (GB)") +
  common_theme +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# Combine into 2x2 grid with shared legend
p_combined <- (p_tpa_size + p_tracepoints) / (p_decode_time + p_peak_memory) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Add panel labels
p_combined <- p_combined +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Save
ggsave(file.path(fig_dir, "figS4_primate_performance.png"),
       p_combined, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS4_primate_performance.pdf"),
       p_combined, width = 12, height = 8, bg = "white", device = cairo_pdf)

# Variant for reviewer: all four panels on log10 scale
p_tpa_size_log <- p_tpa_size + scale_y_log10(expand = expansion(mult = c(0, 0.15)))
p_tracepoints_log <- p_tracepoints + scale_y_log10(expand = expansion(mult = c(0, 0.15)))
p_combined_alllog <- (p_tpa_size_log + p_tracepoints_log) / (p_decode_time + p_peak_memory) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_combined_alllog <- p_combined_alllog +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(fig_dir, "figS4_primate_performance_alllog.png"),
       p_combined_alllog, width = 12, height = 8, dpi = 300, bg = "white")

message("Figure saved: figS4_primate_performance.png/pdf")

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
df %>%
  group_by(method, mc) %>%
  summarise(
    tpa_size_mb = sum(tpa_size_mb),
    tracepoints_m = sum(tracepoints_m),
    decode_time_sec = sum(decode_time_sec),
    avg_memory_gb = mean(peak_memory_gb),
    .groups = "drop"
  ) %>%
  print()
