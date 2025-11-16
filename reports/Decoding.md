# Information about decoding

## Decoding Steps

The PAF file which was encoded using cigzip, is decoded using the following step. At the same time we are calculating the decoding stats and verify the diff with the original PAF.

### cigzip decode + diff verification

**Command**:
```bash
\time -v cigzip decode -p cigzip-datasets/dataset.standard.<N>.tp.mc32.paf --type standard --complexity-metric edit-distance --sequence-files wf-dataset/dataset.<N>.fa --max-complexity 32 -t 8 > cigzip-datasets/dataset.standard.<N>.tp.mc32.decompressed.paf && diff <(sort cigzip-datasets/dataset.standard.<N>.paf) <(sort cigzip-datasets/dataset.standard.<N>.tp.mc32.decompressed.paf) | wc -l
```

### Decoding Stats

| Dataset | Record Count | Sequence Length | Error | Tracepoint Type | CPU Time (s) | Average CPU Time (Per Aignment) (ms) | Run Time (s) | Average Runtime (per aignment) (ms) | Peak Memory (KB) | 
|---------|--------------|-----------------|-------|-----------------|--------------|--------------------------------------|--------------|-------------------------------------|------------------|
| N1 | 10000 | 100   | 0.001 | Standard | 367.93 | 36.79 | 46.22 | 4.62 | 33.4 |
| N2 | 10000 | 100   | 0.01  | Standard | 367.92 | 36.79 | 46.23 | 4.62 | 35.5 |
| N3 | 10000 | 100   | 0.05  | Standard | 254.52 | 25.45 | 31.95 | 3.20 | 33.4 |
| N4 | 10000 | 100   | 0.1   | Standard | 254.10 | 25.41 | 31.93 | 3.19 | 33.4 |
| N5 | 10000 | 1000  | 0.001 | Standard | 523.24 | 52.32 | 65.74 | 6.57 | 33.4 |
| N6 | 10000 | 1000  | 0.01  | Standard | 524.30 | 52.43 | 65.85 | 6.59 | 35.4 |
| N7 | 10000 | 1000  | 0.05  | Standard | 524.33 | 52.43 | 65.88 | 6.59 | 33.6 |
| N8 | 10000 | 1000  | 0.1   | Standard | 383.52 | 38.35 | 48.20 | 4.82 | 34.3 |
| N9 | 10000 | 10000 | 0.001 | Standard | 551.58 | 55.16 | 69.27 | 6.93 | 33.6 |
| N10 | 10000 | 10000 | 0.01  | Standard | 552.66 | 55.27 | 69.39 | 6.94 | 34.9 |
| N11 | 10000 | 10000 | 0.05  | Standard | 558.39 | 55.84 | 70.19 | 7.02 | 35.2 |
| N12 | 10000 | 10000 | 0.1   | Standard | 565.66 | 56.57 | 71.14 | 7.11 | 37.3 |

### File Sizes

| Dataset | Record Count | Sequence Length | Error | Tracepoint Type | Seq File (MB) | FASTA file (MB) | FastGA PAF (MB) | Tracepoint Encoded PAF (MB) |
|---------|--------------|-----------------|-------|-----------------|---------------|-----------------|-----------------|-----------------------------|
| N1 | 10000 | 100   | 0.001 | Standard | 1.95 | 2.16 | 2.33 | 3.49 |
| N2 | 10000 | 100   | 0.01  | Standard | 1.95 | 2.16 | 2.33 | 3.49 |
| N3 | 10000 | 100   | 0.05  | Standard | 1.95 | 2.16 | 1.60 | 2.41 |
| N4 | 10000 | 100   | 0.1   | Standard | 1.95 | 2.16 | 1.59 | 2.40 |
| N5 | 10000 | 1000  | 0.001 | Standard | 19.11 | 19.33 | 3.54 | 5.14 |
| N6 | 10000 | 1000  | 0.01  | Standard | 19.11 | 19.33 | 4.41 | 5.14 |
| N7 | 10000 | 1000  | 0.05  | Standard | 19.11 | 19.33 | 7.51 | 5.26 |
| N8 | 10000 | 1000  | 0.1   | Standard | 19.11 | 19.33 | 5.81 | 3.87 |
| N9 | 10000 | 10000 | 0.001 | Standard | 190.78 | 191.00 | 4.88 | 5.45 |
| N10 | 10000 | 10000 | 0.01  | Standard | 190.78 | 191.00 | 13.60 | 5.94 |
| N11 | 10000 | 10000 | 0.05  | Standard | 190.78 | 191.00 | 44.94 | 7.58 |
| N12 | 10000 | 10000 | 0.1   | Standard | 190.78 | 191.00 | 76.79 | 9.69 |
