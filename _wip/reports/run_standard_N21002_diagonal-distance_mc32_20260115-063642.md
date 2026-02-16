# cigzip / FastGA run

- Timestamp: `2026-01-15T06:36:42-06:00`

## Inputs

| Parameter | Value |
|---|---|
| record_count (-r) | 10000 |
| length (-L) | 10000 |
| error (-E) | 0.1 |
| N (-N) | 21002 |
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
../WFA2-lib/bin/generate_dataset1 -n 10000 -l 10000 -e 0.1 -o data/simulated/wf-dataset/sample.dataset.21002.seq --unbalanced-indels 25,5,6 --diagonal-bandwidth 0.08
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
}' data/simulated/wf-dataset/sample.dataset.21002.seq > data/simulated/wf-dataset/dataset.21002.fa
```
**Output**:
```text
FASTA written: data/simulated/wf-dataset/dataset.21002.fa
```

### 3) FastGA self-alignment → data/simulated/fastga-datasets/dataset.21002.paf

**Command**:
```bash
../FASTGA/FastGA -P/home/hasitha/data/projects/fastga-tmp -pafx data/simulated/wf-dataset/dataset.21002.fa data/simulated/wf-dataset/dataset.21002.fa > data/simulated/fastga-datasets/dataset.21002.paf
```
**Output**:
```text
```

### 4) cigzip encode (type=standard) → data/simulated/tmp/tmp.standard.21002.paf

**Command**:
```bash
../cigzip/target/debug/cigzip encode -p data/simulated/fastga-datasets/dataset.21002.paf --type standard --complexity-metric edit-distance --max-complexity 99999999 > data/simulated/tmp/tmp.standard.21002.paf
```
**Output**:
```text
```

### 5) cigzip decode (type=standard) → data/simulated/cigzip-datasets/dataset.standard.21002.paf

**Command**:
```bash
../cigzip/target/debug/cigzip decode -p data/simulated/tmp/tmp.standard.21002.paf --type standard --complexity-metric edit-distance --sequence-files data/simulated/wf-dataset/dataset.21002.fa > data/simulated/cigzip-datasets/dataset.standard.21002.paf
```
**Output**:
```text
```

### 6) cigzip encode (max-complexity=32, threads=8) --minimal → cigzip-datasets/dataset.standard.21002.tp.mc32.paf

**Command**:
```bash
/usr/bin/time ../cigzip/target/debug/cigzip encode -p data/simulated/cigzip-datasets/dataset.standard.21002.paf --type standard --complexity-metric diagonal-distance --max-complexity 32 -t 8 --minimal > data/simulated/cigzip-datasets/dataset.standard.21002.tp.mc32.paf
```
**Output**:
```text
1.12user 0.27system 0:00.26elapsed 520%CPU (0avgtext+0avgdata 8816maxresident)k
0inputs+8336outputs (1major+524minor)pagefaults 0swaps
```

### 7) cigzip decode + diff verification

**Command**:
```bash
/usr/bin/time -v ../cigzip/target/debug/cigzip decode -p data/simulated/cigzip-datasets/dataset.standard.21002.tp.mc32.paf --type standard --complexity-metric diagonal-distance --sequence-files data/simulated/wf-dataset/dataset.21002.fa --max-complexity 32 -t 8 > data/simulated/cigzip-datasets/dataset.standard.21002.tp.mc32.decompressed.paf && diff <(sort data/simulated/cigzip-datasets/dataset.standard.21002.paf) <(sort data/simulated/cigzip-datasets/dataset.standard.21002.tp.mc32.decompressed.paf) | wc -l
```
**Output**:
```text
[2026-01-15T12:39:29Z WARN  cigzip] Ignoring --max-complexity=32 because --heuristic was not requested
	Command being timed: "../cigzip/target/debug/cigzip decode -p data/simulated/cigzip-datasets/dataset.standard.21002.tp.mc32.paf --type standard --complexity-metric diagonal-distance --sequence-files data/simulated/wf-dataset/dataset.21002.fa --max-complexity 32 -t 8"
	User time (seconds): 570.17
	System time (seconds): 14.26
	Percent of CPU this job got: 795%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:13.47
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 40252
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 2
	Minor (reclaiming a frame) page faults: 114969
	Voluntary context switches: 1976
	Involuntary context switches: 959
	Swaps: 0
	File system inputs: 0
	File system outputs: 20576
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
| data/simulated/wf-dataset/sample.dataset.21002.seq | 200741081 |
| data/simulated/wf-dataset/dataset.21002.fa | 200968869 |
| data/simulated/fastga-datasets/dataset.21002.paf | 7373966 |
| data/simulated/tmp/tmp.standard.21002.paf | 5781671 |
| data/simulated/cigzip-datasets/dataset.standard.21002.paf | 10470083 |
| cigzip-datasets/dataset.standard.21002.tp.mc32.paf | 4267306 |
| cigzip-datasets/dataset.standard.21002.tp.mc32.decompressed.paf | 10470083 |

## Artifacts MB


| data/simulated/wf-dataset/sample.dataset.21002.seq | data/simulated/wf-dataset/dataset.21002.fa | data/simulated/fastga-datasets/dataset.21002.paf | data/simulated/tmp/tmp.standard.21002.paf | data/simulated/cigzip-datasets/dataset.standard.21002.paf | cigzip-datasets/dataset.standard.21002.tp.mc32.paf | cigzip-datasets/dataset.standard.21002.tp.mc32.decompressed.paf |
|---|---|---|---|---|---|---|
| 191.44 | 191.66 | 7.03 | 5.51 | 9.99 | 4.07 | 9.99 |
