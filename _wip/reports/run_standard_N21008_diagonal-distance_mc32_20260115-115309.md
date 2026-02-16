# cigzip / FastGA run

- Timestamp: `2026-01-15T11:53:09-06:00`

## Inputs

| Parameter | Value |
|---|---|
| record_count (-r) | 10000 |
| length (-L) | 10000 |
| error (-E) | 0.1 |
| N (-N) | 21008 |
| type (-T) | standard |
| complexity metric (-C) | diagonal-distance |
| max complexity (-M) | 32 |
| threads (-t) | 8 |

## Tool Binaries

```text
WFA2_BIN   : ../WFA2-lib/bin/generate_dataset1
FASTGA_BIN : ../FASTGA/FastGA
CIGZIP_BIN : ../cigzip/target/debug/cigzip

# WFA2 help/version (if any):
USE: ./generate-datasets [OPTIONS]...
      Options::
        --output|o         PATH        Output path of the generated sequences
        --num-patterns|n   INT         Total number of sequence-pairs generated
        --length|l         INT         Length of the pattern (pattern.length) 
(no output)

# FastGA banner (no -h):

Usage: FastGA [-vkMS] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR)>] [<format(-paf)>]
              [-f<int(10)>] [-c<int(85)> [-s<int(1000)>] [-l<int(100)>] [-i<float(.7)]
              <source1:path>[<precursor>] [<source2:path>[<precursor>]]

(no output)

# cigzip help (if any):
Usage: cigzip <COMMAND>

Commands:
  encode  Encode alignments into tracepoints
  decode  Decode tracepoints back to CIGAR
```

## Pipeline


### 1) WFA2 generate_dataset

**Command**:
```bash
../WFA2-lib/bin/generate_dataset1 -n 10000 -l 10000 -e 0.1 -o data/simulated/wf-dataset/sample.dataset.21008.seq --unbalanced-indels 100,1,6 --diagonal-bandwidth 0.8
```
**Output**:
```text
```

### 2) Convert to FASTA (awk)

**Command**:
```bash
awk '
BEGIN { counter=1 }
/^>/ { 
    print ">pattern_" counter
    print substr($0, 2)
}
/^</ { 
    print ">text_" counter
    print substr($0, 2)
    counter++
}' data/simulated/wf-dataset/sample.dataset.21008.seq > data/simulated/wf-dataset/dataset.21008.fa
```
**Output**:
```text
FASTA written: data/simulated/wf-dataset/dataset.21008.fa
```

### 3) FastGA self-alignment → data/simulated/fastga-datasets/dataset.21008.paf

**Command**:
```bash
../FASTGA/FastGA -P/home/hasitha/data/projects/fastga-tmp -pafx data/simulated/wf-dataset/dataset.21008.fa data/simulated/wf-dataset/dataset.21008.fa > data/simulated/fastga-datasets/dataset.21008.paf
```
**Output**:
```text
```

### 4) cigzip encode (type=standard) → data/simulated/tmp/tmp.standard.21008.paf

**Command**:
```bash
../cigzip/target/debug/cigzip encode -p data/simulated/fastga-datasets/dataset.21008.paf --type standard --complexity-metric edit-distance --max-complexity 99999999 > data/simulated/tmp/tmp.standard.21008.paf
```
**Output**:
```text
```

### 5) cigzip decode (type=standard) → data/simulated/cigzip-datasets/dataset.standard.21008.paf

**Command**:
```bash
../cigzip/target/debug/cigzip decode -p data/simulated/tmp/tmp.standard.21008.paf --type standard --complexity-metric edit-distance --sequence-files data/simulated/wf-dataset/dataset.21008.fa > data/simulated/cigzip-datasets/dataset.standard.21008.paf
```
**Output**:
```text
```

### 6) cigzip encode (max-complexity=32, threads=8) --minimal → cigzip-datasets/dataset.standard.21008.tp.mc32.paf

**Command**:
```bash
/usr/bin/time ../cigzip/target/debug/cigzip encode -p data/simulated/cigzip-datasets/dataset.standard.21008.paf --type standard --complexity-metric diagonal-distance --max-complexity 32 -t 8 --minimal > data/simulated/cigzip-datasets/dataset.standard.21008.tp.mc32.paf
```
**Output**:
```text
2.08user 0.19system 0:00.34elapsed 665%CPU (0avgtext+0avgdata 6988maxresident)k
0inputs+11408outputs (3major+593minor)pagefaults 0swaps
```

### 7) cigzip decode + diff verification

**Command**:
```bash
/usr/bin/time -v ../cigzip/target/debug/cigzip decode -p data/simulated/cigzip-datasets/dataset.standard.21008.tp.mc32.paf --type standard --complexity-metric diagonal-distance --sequence-files data/simulated/wf-dataset/dataset.21008.fa --max-complexity 32 -t 8 > data/simulated/cigzip-datasets/dataset.standard.21008.tp.mc32.decompressed.paf && diff <(sort data/simulated/cigzip-datasets/dataset.standard.21008.paf) <(sort data/simulated/cigzip-datasets/dataset.standard.21008.tp.mc32.decompressed.paf) | wc -l
```
**Output**:
```text
[2026-01-15T17:56:02Z WARN  cigzip] Ignoring --max-complexity=32 because --heuristic was not requested
	Command being timed: "../cigzip/target/debug/cigzip decode -p data/simulated/cigzip-datasets/dataset.standard.21008.tp.mc32.paf --type standard --complexity-metric diagonal-distance --sequence-files data/simulated/wf-dataset/dataset.21008.fa --max-complexity 32 -t 8"
	User time (seconds): 601.18
	System time (seconds): 14.60
	Percent of CPU this job got: 791%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:17.77
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 40720
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 4
	Minor (reclaiming a frame) page faults: 63019
	Voluntary context switches: 2592
	Involuntary context switches: 940
	Swaps: 0
	File system inputs: 0
	File system outputs: 35528
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
0
```

## Artifacts


| File | Size (bytes) |
|---|---|
| data/simulated/wf-dataset/sample.dataset.21008.seq | 203507299 |
| data/simulated/wf-dataset/dataset.21008.fa | 203735087 |
| data/simulated/fastga-datasets/dataset.21008.paf | 12460543 |
| data/simulated/tmp/tmp.standard.21008.paf | 5781435 |
| data/simulated/cigzip-datasets/dataset.standard.21008.paf | 18119153 |
| cigzip-datasets/dataset.standard.21008.tp.mc32.paf | 5838192 |
| cigzip-datasets/dataset.standard.21008.tp.mc32.decompressed.paf | 18119153 |

## Artifacts MB


| data/simulated/wf-dataset/sample.dataset.21008.seq | data/simulated/wf-dataset/dataset.21008.fa | data/simulated/fastga-datasets/dataset.21008.paf | data/simulated/tmp/tmp.standard.21008.paf | data/simulated/cigzip-datasets/dataset.standard.21008.paf | cigzip-datasets/dataset.standard.21008.tp.mc32.paf | cigzip-datasets/dataset.standard.21008.tp.mc32.decompressed.paf |
|---|---|---|---|---|---|---|
| 194.08 | 194.30 | 11.88 | 5.51 | 17.28 | 5.57 | 17.28 |
