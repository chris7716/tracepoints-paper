# cigzip / FastGA run

- Timestamp: `2025-12-28T11:07:18-06:00`

## Inputs

| Parameter | Value |
|---|---|
| record_count (-r) | 10000 |
| length (-L) | 10000 |
| error (-E) | 0.1 |
| N (-N) | 5006 |
| type (-T) | standard |
| complexity metric (-C) | edit-distance |
| max complexity (-M) | 200 |
| threads (-t) | 8 |

## Tool Binaries

```text
WFA2_BIN   : ../WFA2-lib/bin/generate_dataset
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
../WFA2-lib/bin/generate_dataset -n 10000 -l 10000 -e 0.1 -o data/simulated/wf-dataset/sample.dataset.5006.seq
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
}' data/simulated/wf-dataset/sample.dataset.5006.seq > data/simulated/wf-dataset/dataset.5006.fa
```
**Output**:
```text
FASTA written: data/simulated/wf-dataset/dataset.5006.fa
```

### 3) FastGA self-alignment → data/simulated/fastga-datasets/dataset.5006.paf

**Command**:
```bash
../FASTGA/FastGA -P/home/hasitha/data/projects/fastga-tmp -pafx data/simulated/wf-dataset/dataset.5006.fa data/simulated/wf-dataset/dataset.5006.fa > data/simulated/fastga-datasets/dataset.5006.paf
```
**Output**:
```text
```

### 4) cigzip encode (type=standard) → data/simulated/tmp/tmp.standard.5006.paf

**Command**:
```bash
../cigzip/target/debug/cigzip encode -p data/simulated/fastga-datasets/dataset.5006.paf --type standard --complexity-metric edit-distance --max-complexity 99999999 > data/simulated/tmp/tmp.standard.5006.paf
```
**Output**:
```text
```

### 5) cigzip decode (type=standard) → data/simulated/cigzip-datasets/dataset.standard.5006.paf

**Command**:
```bash
../cigzip/target/debug/cigzip decode -p data/simulated/tmp/tmp.standard.5006.paf --type standard --complexity-metric edit-distance --sequence-files data/simulated/wf-dataset/dataset.5006.fa > data/simulated/cigzip-datasets/dataset.standard.5006.paf
```
**Output**:
```text
```

### 6) cigzip encode (max-complexity=200, threads=8) --minimal → cigzip-datasets/dataset.standard.5006.tp.mc200.paf

**Command**:
```bash
/usr/bin/time ../cigzip/target/debug/cigzip encode -p data/simulated/cigzip-datasets/dataset.standard.5006.paf --type standard --complexity-metric edit-distance --max-complexity 200 -t 8 --minimal > data/simulated/cigzip-datasets/dataset.standard.5006.tp.mc200.paf
```
**Output**:
```text
9.92user 0.19system 0:01.34elapsed 755%CPU (0avgtext+0avgdata 10692maxresident)k
0inputs+8888outputs (0major+5447minor)pagefaults 0swaps
```

### 7) cigzip decode + diff verification

**Command**:
```bash
/usr/bin/time -v ../cigzip/target/debug/cigzip decode -p data/simulated/cigzip-datasets/dataset.standard.5006.tp.mc200.paf --type standard --complexity-metric edit-distance --sequence-files data/simulated/wf-dataset/dataset.5006.fa --max-complexity 200 -t 8 > data/simulated/cigzip-datasets/dataset.standard.5006.tp.mc200.decompressed.paf && diff <(sort data/simulated/cigzip-datasets/dataset.standard.5006.paf) <(sort data/simulated/cigzip-datasets/dataset.standard.5006.tp.mc200.decompressed.paf) | wc -l
```
**Output**:
```text
[2025-12-28T17:11:14Z WARN  cigzip] Ignoring --max-complexity=200 because --heuristic was not requested
	Command being timed: "../cigzip/target/debug/cigzip decode -p data/simulated/cigzip-datasets/dataset.standard.5006.tp.mc200.paf --type standard --complexity-metric edit-distance --sequence-files data/simulated/wf-dataset/dataset.5006.fa --max-complexity 200 -t 8"
	User time (seconds): 636.61
	System time (seconds): 14.97
	Percent of CPU this job got: 794%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:22.04
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 40524
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 3
	Minor (reclaiming a frame) page faults: 65795
	Voluntary context switches: 2142
	Involuntary context switches: 1023
	Swaps: 0
	File system inputs: 0
	File system outputs: 161304
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
| data/simulated/wf-dataset/sample.dataset.5006.seq | 200044415 |
| data/simulated/wf-dataset/dataset.5006.fa | 200272203 |
| data/simulated/fastga-datasets/dataset.5006.paf | 80486010 |
| data/simulated/tmp/tmp.standard.5006.paf | 5708672 |
| data/simulated/cigzip-datasets/dataset.standard.5006.paf | 82467839 |
| cigzip-datasets/dataset.standard.5006.tp.mc200.paf | 4547454 |
| cigzip-datasets/dataset.standard.5006.tp.mc200.decompressed.paf | 82467839 |

## Artifacts MB


| data/simulated/wf-dataset/sample.dataset.5006.seq | data/simulated/wf-dataset/dataset.5006.fa | data/simulated/fastga-datasets/dataset.5006.paf | data/simulated/tmp/tmp.standard.5006.paf | data/simulated/cigzip-datasets/dataset.standard.5006.paf | cigzip-datasets/dataset.standard.5006.tp.mc200.paf | cigzip-datasets/dataset.standard.5006.tp.mc200.decompressed.paf |
|---|---|---|---|---|---|---|
| 190.78 | 190.99 | 76.76 | 5.44 | 78.65 | 4.34 | 78.65 |
