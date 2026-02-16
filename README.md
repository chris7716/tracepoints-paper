# Tracepoints paper

Code and data for reproducing the figures and analyses in the tracepoints paper.

## Structure

- `code.md` — step-by-step commands for all experiments (simulated and real data)
- `data/` — benchmark results (TSV)
  - `simulated-data/` — simulated alignment benchmarks
  - `real-data/` — HPRCv2 human and T2T ape pangenome benchmarks
- `scripts/plotting/` — R and Python scripts for generating figures
- `scripts/` — helper scripts for data processing and statistics

## Reproducing figures

```bash
Rscript scripts/plotting/plot-main-figures.R
Rscript scripts/plotting/plot-supp-figures.R
Rscript scripts/plotting/plot-supp-figure-primates.R
Rscript scripts/compute_compression_stats.R
```

Figures are saved to `paper/figures/` by default. Override with `TRACEPOINTS_FIG_DIR` and `TRACEPOINTS_DATA_DIR` environment variables.
