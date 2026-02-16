# cigzip PAF Benchmark Report

- Timestamp: `2026-01-12T21:00:28-06:00`
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
	User time (seconds): 0.34
	System time (seconds): 0.22
	Percent of CPU this job got: 301%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.18
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 6444
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 461
	Voluntary context switches: 22459
	Involuntary context switches: 8
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
	User time (seconds): 415.45
	System time (seconds): 11.19
	Percent of CPU this job got: 794%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:53.67
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 36124
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 7
	Minor (reclaiming a frame) page faults: 52100
	Voluntary context switches: 1694
	Involuntary context switches: 734
	Swaps: 0
	File system inputs: 0
	File system outputs: 8504
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

### 3) Verification: diff original vs decoded

**Command**:
```bash
diff <(sort "./data/simulated/cigzip-datasets/dataset.standard.1001.paf") <(sort "data/real/dataset.standard.1001.standard.diagonal-distance.mc32.decoded.paf") | wc -l
```
**Output**:
```text
35370
