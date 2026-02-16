# cigzip PAF Benchmark Report

- Timestamp: `2026-01-12T21:18:06-06:00`
- Input PAF: `/export2/data/hasitha/projects/tp-datasets/data/simulated/cigzip-datasets/dataset.standard.1001.paf`
- Input FASTA: `/export2/data/hasitha/projects/tp-datasets/data/simulated/wf-dataset/dataset.10001.fa`

## Input Parameters

| Parameter | Value |
|---|---|
| PAF file (-p) | ./data/simulated/cigzip-datasets/dataset.standard.1001.paf |
| FASTA file (-f) | data/simulated/wf-dataset/dataset.10001.fa |
| Benchmark ID (-N) | 25001 |
| Complexity metric (-C) | diagonal-distance |
| Max complexity (-M) | 32 |
| Tracepoint type (-T) | standard |
| Threads (-t) | 8 |
| Output directory (-o) | data/real |

## Input File Information

| File | Size (bytes) | Size (MB) | Records/Lines |
|---|---|---|---|
| ./data/simulated/cigzip-datasets/dataset.standard.1001.paf | 3620827 | 3.45 | 29032 |
| data/simulated/wf-dataset/dataset.10001.fa | 220708944 | 210.48 | 20000 |

## Tool Information

```text
CIGZIP_BIN : ../cigzip/target/debug/cigzip
TIME_BIN   : /usr/bin/time -v

# cigzip version/help:
cigzip 0.1.0
```

## Benchmark Pipeline


### 1) cigzip encode (baseline, no complexity limit)

**Command**:
```bash
/usr/bin/time -v ../cigzip/target/debug/cigzip encode -p "./data/simulated/cigzip-datasets/dataset.standard.1001.paf" --type "standard" --complexity-metric "diagonal-distance" --max-complexity 999999999 -t "8" > "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.encoded.paf"
```
**Output**:
```text
	Command being timed: "../cigzip/target/debug/cigzip encode -p ./data/simulated/cigzip-datasets/dataset.standard.1001.paf --type standard --complexity-metric diagonal-distance --max-complexity 999999999 -t 8"
	User time (seconds): 0.32
	System time (seconds): 0.22
	Percent of CPU this job got: 300%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.18
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 6528
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 1
	Minor (reclaiming a frame) page faults: 464
	Voluntary context switches: 22253
	Involuntary context switches: 32
	Swaps: 0
	File system inputs: 0
	File system outputs: 7176
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

### 2) cigzip decode (verification)

**Command**:
```bash
/usr/bin/time -v ../cigzip/target/debug/cigzip decode -p "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.encoded.paf" --type "standard" --complexity-metric "diagonal-distance" --sequence-files "data/simulated/wf-dataset/dataset.10001.fa" -t "8" > "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.decoded.paf"
```
**Output**:
```text
	Command being timed: "../cigzip/target/debug/cigzip decode -p data/real/dataset.standard.1001.standard.diagonal-distance.mc32.encoded.paf --type standard --complexity-metric diagonal-distance --sequence-files data/simulated/wf-dataset/dataset.10001.fa -t 8"
	User time (seconds): 427.33
	System time (seconds): 11.20
	Percent of CPU this job got: 793%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:55.25
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 36256
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 52901
	Voluntary context switches: 1517
	Involuntary context switches: 673
	Swaps: 0
	File system inputs: 0
	File system outputs: 8504
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

### 4) cigzip encode with max-complexity=32 (main benchmark)

**Command**:
```bash
/usr/bin/time -v ../cigzip/target/debug/cigzip encode -p "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.decoded.paf" --type "standard" --complexity-metric "diagonal-distance" --minimal --max-complexity "32" -t "8" > "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.minimal.paf"
```
**Output**:
```text
	Command being timed: "../cigzip/target/debug/cigzip encode -p data/real/dataset.standard.1001.standard.diagonal-distance.mc32.decoded.paf --type standard --complexity-metric diagonal-distance --minimal --max-complexity 32 -t 8"
	User time (seconds): 0.34
	System time (seconds): 0.15
	Percent of CPU this job got: 360%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.13
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 6424
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 487
	Voluntary context switches: 23630
	Involuntary context switches: 4
	Swaps: 0
	File system inputs: 0
	File system outputs: 4496
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

### 5) cigzip decode with max-complexity=32

**Command**:
```bash
/usr/bin/time -v ../cigzip/target/debug/cigzip decode -p "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.minimal.paf" --type "standard" --complexity-metric "diagonal-distance" --sequence-files "data/simulated/wf-dataset/dataset.10001.fa" --max-complexity "32" -t "8" > "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.minimal.decoded.paf"
```
**Output**:
```text
[2026-01-13T03:19:02Z WARN  cigzip] Ignoring --max-complexity=32 because --heuristic was not requested
	Command being timed: "../cigzip/target/debug/cigzip decode -p data/real/dataset.standard.1001.standard.diagonal-distance.mc32.minimal.paf --type standard --complexity-metric diagonal-distance --sequence-files data/simulated/wf-dataset/dataset.10001.fa --max-complexity 32 -t 8"
	User time (seconds): 414.75
	System time (seconds): 9.19
	Percent of CPU this job got: 795%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:53.27
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 38620
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 129375
	Voluntary context switches: 1719
	Involuntary context switches: 915
	Swaps: 0
	File system inputs: 0
	File system outputs: 8544
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

### 6) Verification: diff original vs max-complexity decoded

**Command**:
```bash
diff <(sort "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.decoded.paf") <(sort "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.minimal.decoded.paf") | wc -l
```
**Output**:
```text
0
```

### 7) Compression Analysis

**Analysis**:
```text
Original PAF size: 3620827 bytes (3.45 MB)
Encoded size (no limit): 3672612 bytes (3.50 MB)
Encoded size (mc=32): 2298822 bytes (2.19 MB)

Compression ratios:
- No limit: 1.0143 (-1.43% reduction)
- Max complexity 32: 0.6349 (36.51% reduction)

Space savings with max complexity: 1.31 MB
```

## Performance Summary

| Metric | Value |
|---|---|
| Original PAF Records | 29032 |
| Original Size (MB) | 3.45 |
| Encoded Size (no limit, MB) | 3.50 |
| Encoded Size (mc=32, MB) | 2.19 |
| Compression Ratio (no limit) | 1.0143 |
| Compression Ratio (mc=32) | 0.6349 |
| Complexity Metric | diagonal-distance |
| Max Complexity Threshold | 32 |
| Tracepoint Type | standard |
| Threads Used | 8 |

## Output Files

| File | Description | Size (bytes) | Size (MB) |
|---|---|---|---|
| data/real/dataset.standard.1001.standard.diagonal-distance.mc32.encoded.paf | Encoded (no complexity limit) | 3672612 | 3.50 |
| data/real/dataset.standard.1001.standard.diagonal-distance.mc32.decoded.paf | Decoded (verification) | 4343674 | 4.14 |
| data/real/dataset.standard.1001.standard.diagonal-distance.mc32.minimal.paf | Encoded (mc=32) | 2298822 | 2.19 |
| data/real/dataset.standard.1001.standard.diagonal-distance.mc32.minimal.decoded.paf | Decoded (mc=32) | 4343674 | 4.14 |
