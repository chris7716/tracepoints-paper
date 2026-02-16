# Tracepoints Paper

This repository contains benchmarking scripts and reports for evaluating the performance of encoding/decoding algorithms in comparison with different kind of alignments.

## Overview

The benchmark pipeline performs the following steps:
1. Generate sequence datasets using generate_dataset tool included in WFA2-lib
2. Convert sequences to FASTA format
3. Perform self-alignment using FastGA
4. Encode alignments to tracepoints using lib_tracepoints library. We use the cigzip tool to interact with lib_tracepoints library.
5. Decode tracepoints back to alignments
6. Verify correctness through diff comparison

## Project Structure

```
├── run_steps.sh              # Main execution script
├── run_steps_report.sh       # Report-generating version
├── cigzip-datasets/          # Output PAF files from cigzip
├── fastga-datasets/          # Output PAF files from FastGA
├── wf-dataset/              # Generated sequences and FASTA files
├── tmp/                     # Intermediate files
└── reports/                 # Benchmark reports
```

## Dependencies

- **WFA2**: `../WFA2-lib/bin/generate_dataset`
- **FastGA**: `../FASTGA/FastGA`  
- **cigzip**: `../cigzip/target/debug/cigzip`

## Usage

### Generate Report

```bash
./run_steps_report.sh -r <records> -L <length> -E <error> -N <id> -T <type>
```

### Parameters

| Parameter | Description | Example Values |
|-----------|-------------|----------------|
| `-r` | Number of sequence pairs to generate | `10000` |
| `-L` | Sequence length | `100`, `1000`, `10000` |
| `-E` | Error rate | `0.001`, `0.01`, `0.05`, `0.1` |
| `-N` | Dataset identifier | `1`, `2`, `3`, etc. |
| `-T` | Encoding type | `standard` |

## Configuration

The scripts use the following default settings:
- **Complexity metric**: `edit-distance`
- **Max complexity**: `32`
- **Threads**: `8`
- **FastGA temp directory**: `/home/hasitha/data/projects/fastga-tmp`

## Generated Files

For each run with parameters `-N X`, the following files are created:

| File Pattern | Description |
|--------------|-------------|
| `wf-dataset/sample.dataset.X.seq` | Raw WFA2 output |
| `wf-dataset/dataset.X.fa` | FASTA sequences |
| `fastga-datasets/dataset.X.paf` | FastGA alignments |
| `tmp/tmp.standard.X.paf` | cigzip intermediate |
| `cigzip-datasets/dataset.standard.X.paf` | cigzip decoded alignments |
| `cigzip-datasets/dataset.standard.X.tp.mc32.paf` | cigzip tracepoints |
| `cigzip-datasets/dataset.standard.X.tp.mc32.decompressed.paf` | Verified output |

## Benchmark Results

The [reports](reports/) directory contains detailed benchmark results for various parameter combinations:

### Dataset Parameters Tested

| Dataset | Records | Length | Error Rate | Report |
|---------|---------|--------|------------|---------|
| N1 | 10,000 | 100 | 0.001 | [run_standard_N1_20251109-030004.md](reports/run_standard_N1_20251109-030004.md) |
| N2 | 10,000 | 100 | 0.01 | [run_standard_N2_20251109-030652.md](reports/run_standard_N2_20251109-030652.md) |
| N3 | 10,000 | 100 | 0.05 | [run_standard_N3_20251109-032155.md](reports/run_standard_N3_20251109-032155.md) |
| N4 | 10,000 | 100 | 0.1 | [run_standard_N4_20251109-033201.md](reports/run_standard_N4_20251109-033201.md) |
| N5 | 10,000 | 1,000 | 0.001 | [run_standard_N5_20251109-033931.md](reports/run_standard_N5_20251109-033931.md) |
| N6 | 10,000 | 1,000 | 0.01 | [run_standard_N6_20251109-034744.md](reports/run_standard_N6_20251109-034744.md) |
| N7 | 10,000 | 1,000 | 0.05 | [run_standard_N7_20251109-035521.md](reports/run_standard_N7_20251109-035521.md) |
| N8 | 10,000 | 1,000 | 0.1 | [run_standard_N8_20251109-040902.md](reports/run_standard_N8_20251109-040902.md) |
| N9 | 10,000 | 10,000 | 0.001 | [run_standard_N9_20251109-041450.md](reports/run_standard_N9_20251109-041450.md) |
| N10 | 10,000 | 10,000 | 0.01 | [run_standard_N10_20251109-045054.md](reports/run_standard_N10_20251109-045054.md) |
| N11 | 10,000 | 10,000 | 0.05 | [run_standard_N11_20251109-045709.md](reports/run_standard_N11_20251109-045709.md) |
| N12 | 10,000 | 10,000 | 0.1 | [reports/run_standard_N12_20251109-050305.md](reports/run_standard_N12_20251109-050305.md) |

## Report Format

Each report includes:
- **Input Parameters**: Dataset configuration
- **Tool Binaries**: Version information for all tools
- **Pipeline Steps**: Detailed command execution with outputs
- **Performance Metrics**: Runtime and memory usage via `/usr/bin/time -v`
- **File Sizes**: Compression analysis of intermediate and final outputs
- **Verification**: Diff results confirming correctness

## Pipeline Steps Detail

1. **WFA2 Dataset Generation**: Creates synthetic sequence pairs with specified error rates
2. **FASTA Conversion**: Transforms WFA2 format to standard FASTA
3. **FastGA Self-Alignment**: Generates PAF alignments from sequences
4. **cigzip Standard Encoding**: Converts PAF to tracepoint format
5. **cigzip Standard Decoding**: Reconstructs PAF from tracepoints
6. **cigzip Tracepoint Encoding**: Creates compressed tracepoints (max-complexity=32)
7. **cigzip Verification**: Decodes and verifies against original alignments

## Example Usage

Generate a benchmark for 1000bp sequences with 1% error rate:

```bash
./run_steps_report.sh -r 10000 -L 1000 -E 0.01 -N 6 -T standard
```

This creates report: `reports/run_standard_N6_<timestamp>.md`

## Performance Analysis

The benchmark evaluates:
- **Compression ratio**: Comparison of file sizes across pipeline stages
- **Runtime performance**: Encoding/decoding speed with multi-threading
- **Memory usage**: Peak memory consumption during processing
- **Correctness**: Verification that decoded alignments match originals