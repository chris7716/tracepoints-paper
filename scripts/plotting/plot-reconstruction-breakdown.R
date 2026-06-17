#!/usr/bin/env Rscript

# Reviewer request: show how reconstruction (TPA -> PAF) time splits into its two
# stages for each case: decompression (reading/decoding the tracepoint streams) and
# re-alignment (per-segment WFA to rebuild the CIGAR). Stacked bars, default params
# (FL-TP l=100, EB-TP delta=32, DB-TP b=32), faceted method x sequence length.
#
# Input:  data/simulated-data/benchmark.results.tsv
# Output: paper/figures/figS9_reconstruction_breakdown.{png,pdf}
# Deps:   tidyverse, ggh4x

suppressPackageStartupMessages({ library(tidyverse); library(ggh4x) })

data_dir <- Sys.getenv("TRACEPOINTS_DATA_DIR", unset = "data")
fig_dir  <- Sys.getenv("TRACEPOINTS_FIG_DIR", unset = "paper/figures")
data <- read_tsv(file.path(data_dir, "simulated-data", "benchmark.results.tsv"), show_col_types = FALSE)

mk <- function(d) d %>% mutate(
  length_label = factor(case_when(l==100~"100 bp", l==1000~"1 Kbp", l==10000~"10 Kbp", l==100000~"100 Kbp"),
                        levels = c("100 bp","1 Kbp","10 Kbp","100 Kbp")),
  error_label = factor(paste0(e*100,"%"), levels = c("0.1%","1%","5%","10%","20%")))

# Re-alignment time = decode (FL-TP) / banded heuristic decode (EB-TP, DB-TP), matching
# the reconstruction reported in Figure 3.
df <- data %>% filter(memory_mode == "high") %>%
  filter((tp_type == "fastga" & mc == 100) | (tp_type == "standard" & mc == 32)) %>%
  mutate(method = factor(case_when(tp_type == "fastga" ~ "FL-TP TPA",
                                   cm == "edit-distance" ~ "EB-TP TPA",
                                   cm == "diagonal-distance" ~ "DB-TP TPA"),
                         levels = c("FL-TP TPA","EB-TP TPA","DB-TP TPA")),
         Decompression  = decompress_runtime_sec,
         `Re-alignment` = if_else(tp_type == "fastga", decode_runtime_sec, decode_heuristic_runtime_sec)) %>%
  mk()

long <- df %>% select(method, length_label, error_label, Decompression, `Re-alignment`) %>%
  pivot_longer(c(Decompression, `Re-alignment`), names_to = "stage", values_to = "t") %>%
  mutate(stage = factor(stage, levels = c("Re-alignment","Decompression")))

p <- ggplot(long, aes(error_label, t, fill = stage)) +
  geom_col(width = 0.7) +
  ggh4x::facet_grid2(method ~ length_label, scales = "free_y", independent = "y") +
  scale_fill_manual(values = c("Decompression" = "#bdbdbd", "Re-alignment" = "#4575b4"),
                    name = NULL, breaks = c("Decompression","Re-alignment")) +
  labs(x = "Error rate", y = "Reconstruction time (sec)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "figS9_reconstruction_breakdown.png"), p, width = 11, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figS9_reconstruction_breakdown.pdf"), p, width = 11, height = 7.5, bg = "white")
message("Figure saved: ", file.path(fig_dir, "figS9_reconstruction_breakdown.{png,pdf}"))
