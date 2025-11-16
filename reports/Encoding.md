# Information about decoding

This report explains about the stats on PAF encoding using tracepoints.

### cigzip encode (max-complexity=32, threads=8) → cigzip-datasets/dataset.standard.N.tp.mc32.paf

**Command**:
```bash
cigzip encode -p cigzip-datasets/dataset.standard.N.paf --type standard --complexity-metric edit-distance --max-complexity 32 -t 8 > cigzip-datasets/dataset.standard.N.tp.mc32.paf
```

### Decoding Stats

| Dataset | Record Count | Sequence Length | Error | Tracepoint Type | CPU Time (s) | Average CPU Time (Per Aignment) (ms) | Run Time (s) | Average Runtime (per aignment) (ms) | Peak Memory (KB) | Total Tracepoints | Tracepoints in Bytes |
|---------|--------------|-----------------|-------|-----------------|--------------|--------------------------------------|--------------|-------------------------------------|------------------|--------------|----------|
| N1 | 10000 | 100   | 0.001 | Standard | 0.47 | 0.047 | 0.17 | 0.017 | 6.55 | 57952 | 927232 |
| N2 | 10000 | 100   | 0.01  | Standard | 0.46 | 0.046 | 0.16 | 0.016 | 6.66 | 57876 | 926016 |
| N3 | 10000 | 100   | 0.05  | Standard | 0.31 | 0.031 | 0.10 | 0.010 | 6.56 | 40108 | 641728 |
| N4 | 10000 | 100   | 0.1   | Standard | 0.33 | 0.033 | 0.11 | 0.011 | 6.54 | 40000 | 640000 |
| N5 | 10000 | 1000  | 0.001 | Standard | 0.64 | 0.064 | 0.22 | 0.022 | 8.41 | 80000 | 1280000 |
| N6 | 10000 | 1000  | 0.01  | Standard | 0.74 | 0.074 | 0.22 | 0.022 | 8.46 | 80000 | 1280000 |
| N7 | 10000 | 1000  | 0.05  | Standard | 1.06 | 0.106 | 0.23 | 0.023 | 8.73 | 120000 | 1920000 |
| N8 | 10000 | 1000  | 0.1   | Standard | 0.80 | 0.080 | 0.16 | 0.016 | 6.86 | 80000 | 1280000 |
| N9 | 10000 | 10000 | 0.001 | Standard | 0.78 | 0.078 | 0.25 | 0.025 | 6.55 | 57952 | 927232 |
| N10 | 10000 | 10000 | 0.01  | Standard | 1.62 | 0.162 | 0.27 | 0.027 | 6.98 | 199200 | 3187200 |
| N11 | 10000 | 10000 | 0.05  | Standard | 4.67 | 0.467 | 0.64 | 0.064 | 9.85 | 659748 | 10555968 |
| N12 | 10000 | 10000 | 0.1   | Standard | 7.96 | 0.796 | 1.07 | 0.107 | 11.04 | 1209552 | 19352832 |
