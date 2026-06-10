# Tracepoints paper

Code and data to reproduce the figures and statistics in the
"Adaptive Tracepoints for Pangenome Alignment Compression" publication.

## Structure

- [`code.md`](code.md): step-by-step commands for every experiment (simulated and real data), from sequence generation to TPA encoding.
- [`data/`](data/): benchmark tables (TSV) the figures are built from
  - [`simulated-data/`](data/simulated-data/): simulated alignment benchmarks (Figures 2 and 3)
  - [`real-data/`](data/real-data/): HPRCv2 human and T2T ape pangenome benchmarks (Figure S4, supplementary statistics)
- [`scripts/plotting/`](scripts/plotting/): figure scripts
- [`scripts/`](scripts/): helpers

## Requirements

R with `tidyverse`, `ggplot2`, `ggh4x`, `scales`, `patchwork`, `cowplot`, `gridExtra`:

```r
install.packages(c("tidyverse", "ggh4x", "scales", "patchwork", "cowplot", "gridExtra"))
```

The from-scratch pipeline below additionally needs [WFA2-lib](https://github.com/smarco/WFA2-lib) (`generate_dataset`, `align_benchmark`), [cigzip](https://github.com/AndreaGuarracino/cigzip), and Python 3.

## Reproduce the figures from the included data

Figures 2 (compression ratios) and 3 (reconstruction cost) both come from `data/simulated-data/benchmark.results.tsv`:

```bash
mkdir -p paper/figures
Rscript scripts/plotting/plot-main-figures.R           # -> compression_ratios.pdf (Fig 2), decoding_cost.pdf (Fig 3)
Rscript scripts/plotting/plot-supp-figures.R           # -> Figures S1, S2, S3, S7
Rscript scripts/plotting/plot-supp-figure-primates.R   # -> Figure S4 (T2T apes)
Rscript scripts/compute_compression_stats.R            # prints the statistics cited in the text
```

## Regenerate the simulated benchmark from scratch (Figures 2 and 3)

This rebuilds `data/simulated-data/benchmark.results.tsv` from raw sequences.

```bash
# 1. Generate sequences and optimal alignments for the 16 conditions
#    (4 lengths x 4 divergences, 10,000 pairs each)
for l in 100 1000 10000 100000; do
  for e in 0.01 0.05 0.10 0.20; do
    generate_dataset -n 10000 -l $l -e $e -o set_${l}_${e}.seq
    align_benchmark -i set_${l}_${e}.seq -a edit-wfa --wfa-memory high \
      --output set_${l}_${e}.out -t 1
    python3 scripts/out-to-paf.py set_${l}_${e}.out set_${l}_${e}.paf
  done
done

# 2. Encode/compress/decode each PAF with cigzip and collect the per-condition
#    metrics into benchmark.results.tsv. The full benchmarking loop (all methods
#    and thresholds, with timing and memory) is in the "Simulated data" section
#    of code.md.

# 3. Plot Figures 2 and 3 from the regenerated table
mkdir -p paper/figures
Rscript scripts/plotting/plot-main-figures.R
```
