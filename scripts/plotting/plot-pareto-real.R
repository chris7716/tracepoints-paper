#!/usr/bin/env Rscript
# Figure 3: Pareto frontier on the real pangenomes.
# Table 2 definitions: time summed over files,
# memory = highest peak (decompress or decode), gap-affine2p reconstruction.
# Rows: reconstruction time, peak memory; columns: human, primates; x = TPA size.
# Dashed: frontier; main labels frontier only
# Usage: Rscript scripts/plotting/plot-pareto-real.R [all]   (from the repo root)
# Output: paper/figures/pareto_real.pdf (l: 6 values, delta/b: 16, 32, 128, 256, 512, 1024; one column)
#         paper/figures/pareto_real_all.pdf (every tested value, Supp Fig S13; only with 'all')
suppressPackageStartupMessages({ library(tidyverse); library(ggrepel); library(ggh4x) })
data_dir <- Sys.getenv("TRACEPOINTS_DATA_DIR", unset = "data")
fig_dir  <- Sys.getenv("TRACEPOINTS_FIG_DIR", unset = "paper/figures")

# ggrepel 0.9.6 parse: "" must stay an obstacle
assignInNamespace("not_empty", function(xs) vapply(seq_along(xs), function(i)
  if (is.expression(xs)) !(length(xs[[i]]) == 1 && is.na(xs[[i]])) else xs[i] != "", logical(1)), "ggrepel")

lev  <- c("EB-TP", "DB-TP", "FL-TP", "FL-TP nd")
cols <- c("EB-TP" = "#d95f02", "DB-TP" = "#377eb8", "FL-TP" = "#7570b3", "FL-TP nd" = "#9e9ac8")   # plot-main-figures.R
shp  <- c("EB-TP" = 16, "DB-TP" = 17, "FL-TP" = 15, "FL-TP nd" = 0)
meth <- c("edit-distance" = "EB-TP", "diagonal-distance" = "DB-TP", "fastga" = "FL-TP", "fastga-no-diff" = "FL-TP nd")

# One pangenome: totals per setting
summ <- function(f, label) {
  message("Reading: ", f)
  d <- read_tsv(f, show_col_types = FALSE) %>% filter(distance == "gap-affine2p", cm %in% names(meth)) %>%
    mutate(target = sub("\\.fltp$", "", paf_file))
  nt <- n_distinct(d$target)
  if (anyNA(d$decode_runtime_sec) || anyNA(d$size_tpa_bytes)) stop(label, ": NA sizes or times")
  s <- d %>% group_by(cm, mc) %>% summarise(n = n_distinct(target), k = n(), x = sum(size_tpa_bytes) / 2^30,
                                            time = sum(decompress_runtime_sec + decode_runtime_sec) / 3600,
                                            mem = max(pmax(decompress_memory_kb, decode_memory_kb)) / 2^20, .groups = "drop")
  bad <- filter(s, n != nt | k != nt)
  if (nrow(bad)) { print(bad); stop(label, ": settings without exactly one row per file (", nt, " files)") }
  message(label, ": ", nt, " files, ", nrow(s), " settings: ",
          paste(s %>% group_by(cm) %>% summarise(v = paste(sort(mc), collapse = ",")) %>% transmute(paste0(meth[cm], " ", v)) %>% pull(), collapse = "; "))
  s %>% mutate(method = factor(meth[cm], levels = lev), ds = label)
}
r <- bind_rows(summ(file.path(data_dir, "real-data/hprcv2-25k.benchmark.results.tsv"), "Human"),
               summ(file.path(data_dir, "real-data/t2t-ape-pangenome.benchmark.results.tsv"), "Primates"))
rows <- c("Time (h)", "Peak memory (GiB)")
# Hand-placed labels (main figure): position, justification
man <- tribble(~ds, ~rowlab, ~method, ~mc, ~lx, ~ly, ~hj, ~vj,
  "Human", rows[1], "DB-TP", 512, 16.5, 36, 0, 0.5,    "Human", rows[1], "DB-TP", 256, 17.5, 26, 0, 0.5,
  "Human", rows[1], "DB-TP", 128, 18.5, 13, 0, 0.5,    "Human", rows[1], "EB-TP", 1024, 13.9, 25, 1, 0.5,
  "Human", rows[1], "EB-TP", 512, 14.0, 14, 1, 0.5,    "Human", rows[1], "EB-TP", 256, 15.6, 7, 1, 0.5,
  "Human", rows[1], "EB-TP", 128, 19.5, 3.3, 1, 0.5,   "Human", rows[2], "DB-TP", 512, 16.5, 110, 0, 0.5,
  "Human", rows[2], "EB-TP", 256, 18.0, 53, 0.5, 0,    "Human", rows[2], "EB-TP", 128, 19.5, 30, 0.5, 1) %>%
  mutate(ds = factor(ds, levels = c("Human", "Primates")), rowlab = factor(rowlab, levels = rows), method = factor(method, levels = lev))
# Hand-placed labels (Supp S13, every setting)
man_all <- tribble(~ds, ~rowlab, ~method, ~mc, ~lx, ~ly, ~hj, ~vj,
  "Human", rows[1], "DB-TP", 512, 16.5, 36, 0, 0.5,    "Human", rows[1], "DB-TP", 256, 17.5, 25, 0, 0.5,
  "Human", rows[1], "DB-TP", 128, 19.5, 14, 0, 0.5,    "Human", rows[1], "DB-TP", 112, 20.0, 11.5, 0, 0.5,
  "Human", rows[1], "DB-TP", 96, 20.5, 9.4, 0, 0.5,    "Human", rows[1], "DB-TP", 80, 21.0, 7.7, 0, 0.5,
  "Human", rows[1], "DB-TP", 64, 21.5, 6.3, 0, 0.5,    "Human", rows[1], "DB-TP", 48, 22.5, 5.15, 0, 0.5,
  "Human", rows[1], "DB-TP", 16, 27.5, 4.1, 0, 0.5,    "Human", rows[1], "EB-TP", 1024, 13.9, 25, 1, 0.5,
  "Human", rows[1], "EB-TP", 512, 14.6, 14, 1, 0.5,    "Human", rows[1], "EB-TP", 256, 16.6, 7, 1, 0.5,
  "Human", rows[1], "EB-TP", 128, 19.5, 3.5, 1, 0.5,   "Human", rows[1], "EB-TP", 112, 20.0, 2.3, 0.5, 0.5,
  "Human", rows[1], "EB-TP", 96, 23.0, 1.9, 0.5, 0.5,   "Human", rows[1], "EB-TP", 80, 26.5, 2.3, 0.5, 0.5,
  "Human", rows[1], "EB-TP", 64, 30.5, 1.9, 0.5, 0.5,   "Human", rows[1], "EB-TP", 48, 35.0, 2.3, 0.5, 0.5,
  "Human", rows[1], "EB-TP", 32, 40.5, 1.9, 0.5, 0.5,   "Human", rows[1], "EB-TP", 16, 55, 2.05, 0.5, 0.5,
  "Primates", rows[1], "DB-TP", 1024, 0.105, 12, 0.5, 0, "Primates", rows[1], "DB-TP", 512, 0.12, 3.4, 1, 0.5,
  "Primates", rows[1], "DB-TP", 256, 0.165, 1.6, 1, 0.5, "Primates", rows[1], "EB-TP", 1024, 0.21, 3.2, 0, 0.5,
  "Primates", rows[1], "DB-TP", 128, 0.205, 1.05, 1, 0.5, "Primates", rows[1], "DB-TP", 112, 0.21, 0.84, 1, 0.5,
  "Primates", rows[1], "DB-TP", 96, 0.22, 0.67, 1, 0.5,  "Primates", rows[1], "DB-TP", 80, 0.235, 0.535, 1, 0.5,
  "Primates", rows[1], "DB-TP", 64, 0.26, 0.43, 1, 0.5,  "Primates", rows[1], "DB-TP", 48, 0.30, 0.345, 1, 0.5,
  "Primates", rows[1], "DB-TP", 32, 0.53, 0.33, 0, 0.5) %>%
  mutate(ds = factor(ds, levels = c("Human", "Primates")), rowlab = factor(rowlab, levels = rows), method = factor(method, levels = lev))
# One figure from per-setting totals
draw <- function(r, out, w, h, hand) {
message("== ", out, ": ", nrow(r), " settings")
d <- bind_rows(r %>% mutate(y = time, rowlab = rows[1]), r %>% mutate(y = mem, rowlab = rows[2])) %>%
  mutate(rowlab = factor(rowlab, levels = rows), ds = factor(ds, levels = c("Human", "Primates")))

# Frontier: no other setting smaller and faster
front <- d %>% group_by(ds, rowlab) %>% arrange(x, y, .by_group = TRUE) %>% filter(y < lag(cummin(y), default = Inf)) %>% ungroup()
message("Frontier settings per panel:"); print(front %>% count(ds, rowlab, method) %>% pivot_wider(names_from = method, values_from = n, values_fill = 0))
print(front %>% select(ds, rowlab, method, mc, x, y), n = Inf)

# Green area under the frontier steps
under <- front %>% group_by(ds, rowlab) %>% arrange(x, .by_group = TRUE) %>%
  mutate(xmin = if_else(row_number() == 1, 0, x), xmax = lead(x, default = Inf)) %>% ungroup()
rng <- d %>% group_by(ds, rowlab) %>% summarise(lx0 = min(log10(x)), lx1 = max(log10(x)), ly0 = min(log10(y)), ly1 = max(log10(y)), .groups = "drop") %>%
  group_by(rowlab) %>% mutate(ly0 = if_else(rowlab == rows[2], min(ly0), ly0), ly1 = if_else(rowlab == rows[2], max(ly1), ly1)) %>% ungroup()
# Obstacles over the "better" text (panel expansions 7%, 12%)
corner <- rng %>% mutate(X0 = lx0 - 0.07 * (lx1 - lx0), X1 = lx1 + 0.07 * (lx1 - lx0), Y0 = ly0 - 0.12 * (ly1 - ly0), Y1 = ly1 + 0.12 * (ly1 - ly0)) %>%
  crossing(f = seq(0.02, 0.18, by = 0.02), g = c(0.03, 0.07)) %>% transmute(ds, rowlab, x = 10^(X0 + f * (X1 - X0)), y = 10^(Y0 + g * (Y1 - Y0)))
shade <- under %>% mutate(ymax = y)

# Labels: frontier settings only
seg <- function(x0, y0, x1, y1, n = 15) { t <- seq(0, 1, length.out = n)[-c(1, n)]
  tibble(x = 10^(log10(x0) + t * (log10(x1) - log10(x0))), y = 10^(log10(y0) + t * (log10(y1) - log10(y0)))) }
along <- d %>% group_by(ds, rowlab, method) %>% arrange(x, .by_group = TRUE) %>%
  reframe(bind_rows(Map(seg, head(x, -1), head(y, -1), tail(x, -1), tail(y, -1))))
steps <- front %>% group_by(ds, rowlab) %>% arrange(x, .by_group = TRUE) %>%
  reframe(bind_rows(Map(seg, head(x, -1), head(y, -1), tail(x, -1), head(y, -1)), Map(seg, tail(x, -1), head(y, -1), tail(x, -1), tail(y, -1))))
# Hand labels: frontier points only, gap at the label end
k <- c("ds", "rowlab", "method", "mc")
mlab <- inner_join(front %>% select(all_of(k), x, y), hand, by = k) %>%
  mutate(xe = 10^(log10(x) + 0.85 * (log10(lx) - log10(x))), ye = 10^(log10(y) + 0.85 * (log10(ly) - log10(y))))
labd <- d %>% left_join(front %>% distinct(ds, rowlab, method, mc) %>% mutate(onf = TRUE), by = c("ds", "rowlab", "method", "mc")) %>%
  mutate(onf = coalesce(onf, FALSE), lab = if_else(onf, as.character(mc), ""), psize = 2.4) %>%
  mutate(lab = if_else(paste(ds, rowlab, method, mc) %in% with(mlab, paste(ds, rowlab, method, mc)), "", lab)) %>%
  bind_rows(along %>% mutate(lab = "", onf = FALSE, psize = 1), steps %>% mutate(lab = "", onf = FALSE, psize = 1, method = factor(NA, levels = lev)),
            corner %>% mutate(lab = "", onf = FALSE, psize = 4, method = factor(NA, levels = lev)),
            mlab %>% transmute(ds, rowlab, x = lx, y = ly, lab = "", onf = FALSE, psize = 3, method = factor(NA, levels = lev))) %>%
  mutate(lab_col = if_else(lab == "", "", paste0("'", lab, "'")))
message("Labels: ", sum(labd$lab != ""), " (", sum(labd$onf & labd$lab != ""), " on the frontier); line obstacles: ", nrow(along) + nrow(steps))

ticks <- unique(signif(as.vector(outer(c(1, 1.5, 2, 3, 5, 7), 10^(-5:6))), 3))
brk <- function(lim) { b <- ticks[ticks >= lim[1] & ticks <= lim[2]]; if (length(b) >= 3) b else pretty(lim, 4) }
# y: 1-2-5 ticks above 12 breaks
brk_y <- function(lim) { b <- brk(lim); if (length(b) > 12) b[signif(b / 10^floor(log10(b)), 2) %in% c(1, 2, 5)] else b }
lab <- function(x) sapply(x, format, scientific = FALSE, drop0trailing = TRUE, big.mark = ",")
# Repelled labels, lines and points as obstacles
rep <- function(col, ...) geom_text_repel(data = labd, aes(label = .data[[col]], point.size = psize), parse = TRUE, size = 3.4, show.legend = FALSE,
                                          min.segment.length = 0.1, segment.size = 0.25, box.padding = 0.5, point.padding = 0.3, force = 4, force_pull = 0.3,
                                          max.iter = 100000, max.time = 60, max.overlaps = Inf, seed = 42, ...)
# "better" text, 6 pt from both axes
better <- grid::textGrob("better", x = unit(6, "pt"), y = unit(6, "pt"), just = c("left", "bottom"), gp = grid::gpar(col = "#3d6612", fontsize = 12, fontface = "bold.italic"))
# Same y range for the memory row
ext <- d %>% filter(rowlab == rows[2]) %>% mutate(y0 = min(y), y1 = max(y)) %>% group_by(rowlab, ds) %>% reframe(x = min(x), y = c(first(y0), first(y1)))
p <- ggplot(d, aes(x, y, color = method, shape = method)) +
  geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax), inherit.aes = FALSE, fill = "#66a61e", alpha = 0.12, show.legend = FALSE) +
  annotation_custom(better, xmin = 0, xmax = Inf, ymin = 0, ymax = Inf) +
  geom_step(data = front, aes(group = interaction(ds, rowlab)), color = "grey40", linetype = "dashed", direction = "hv", linewidth = 0.6, show.legend = FALSE) +
  geom_line(aes(group = method), linewidth = 0.8) + geom_point(size = 2.4) +
  rep("lab_col") +
  geom_segment(data = mlab, aes(x, y, xend = xe, yend = ye, color = method), linewidth = 0.25, inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = mlab, aes(lx, ly, label = mc, hjust = hj, vjust = vj, color = method), size = 3.4, inherit.aes = FALSE, show.legend = FALSE) +
  geom_blank(data = ext, aes(x, y), inherit.aes = FALSE) +
  facet_grid2(rowlab ~ ds, scales = "free", independent = "y", switch = "y") +
  scale_x_log10(breaks = brk, labels = lab, expand = expansion(mult = 0.07)) + scale_y_log10(breaks = brk_y, labels = lab, expand = expansion(mult = 0.12)) +
  scale_color_manual(values = cols, breaks = lev, name = "Method", drop = FALSE) + scale_shape_manual(values = shp, breaks = lev, name = "Method", drop = FALSE) +
  labs(x = "File size (GiB)", y = NULL) +
  # Text sizes as plot-main-figures.R
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom", legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 13),
        legend.key.size = unit(0.5, "cm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"), axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = -1)), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 15, face = "bold"), strip.text = element_text(size = 15, face = "bold"),
        aspect.ratio = 1, strip.placement.x = "inside", strip.placement.y = "outside", panel.spacing.y = unit(2, "pt"), panel.spacing.x = unit(2, "pt"), strip.switch.pad.grid = unit(0, "pt"), strip.background.y = element_blank(), strip.text.y.left = element_text(size = 15, face = "bold", angle = 90, margin = margin(r = 2, l = 0)))
# x = 0 is the panel edge: mute that warning
n_inf <- 0
withCallingHandlers(ggsave(file.path(fig_dir, paste0(out, ".pdf")), p, width = w, height = h, bg = "white"),
  warning = function(w) if (grepl("transformation introduced infinite values", conditionMessage(w))) { n_inf <<- n_inf + 1; invokeRestart("muffleWarning") })
message("Muted ", n_inf, " expected 'infinite values' warnings")
if (system2("bash", c("scripts/plotting/crop-pdf.sh", file.path(fig_dir, paste0(out, ".pdf")))) != 0) stop("crop-pdf.sh failed")
message("Figure saved: ", file.path(fig_dir, paste0(out, ".pdf")))
}
# Main text: subset incl. min and max
six <- r %>% filter(if_else(method %in% c("EB-TP", "DB-TP"), mc %in% c(16, 32, 128, 256, 512, 1024), mc %in% c(100, 200, 500, 1000, 2000, 3000)))
draw(six, "pareto_real", 7.5, 7.5, man)
# Every tested value: only on request
if ("all" %in% commandArgs(TRUE)) draw(r, "pareto_real_all", 7.5, 7.5, man_all) else message("Skipped pareto_real_all (pass 'all' to draw it)")
