# cigzip / FastGA run

- Timestamp: `2026-01-15T11:46:16-06:00`

## Inputs

| Parameter | Value |
|---|---|
| record_count (-r) | 10000 |
| length (-L) | 10000 |
| error (-E) | 0.1 |
| N (-N) | 21007 |
| type (-T) | standard |
| complexity metric (-C) | edit-distance |
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
../WFA2-lib/bin/generate_dataset1 -n 10000 -l 10000 -e 0.1 -o data/simulated/wf-dataset/sample.dataset.21007.seq --unbalanced-indels 100,1,6 --diagonal-bandwidth 0.8
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
}' data/simulated/wf-dataset/sample.dataset.21007.seq > data/simulated/wf-dataset/dataset.21007.fa
```
**Output**:
```text
FASTA written: data/simulated/wf-dataset/dataset.21007.fa
```

### 3) FastGA self-alignment → data/simulated/fastga-datasets/dataset.21007.paf

**Command**:
```bash
../FASTGA/FastGA -P/home/hasitha/data/projects/fastga-tmp -pafx data/simulated/wf-dataset/dataset.21007.fa data/simulated/wf-dataset/dataset.21007.fa > data/simulated/fastga-datasets/dataset.21007.paf
```
**Output**:
```text
```

### 4) cigzip encode (type=standard) → data/simulated/tmp/tmp.standard.21007.paf

**Command**:
```bash
../cigzip/target/debug/cigzip encode -p data/simulated/fastga-datasets/dataset.21007.paf --type standard --complexity-metric edit-distance --max-complexity 99999999 > data/simulated/tmp/tmp.standard.21007.paf
```
**Output**:
```text
```

### 5) cigzip decode (type=standard) → data/simulated/cigzip-datasets/dataset.standard.21007.paf

**Command**:
```bash
../cigzip/target/debug/cigzip decode -p data/simulated/tmp/tmp.standard.21007.paf --type standard --complexity-metric edit-distance --sequence-files data/simulated/wf-dataset/dataset.21007.fa > data/simulated/cigzip-datasets/dataset.standard.21007.paf
```
**Output**:
```text
```

### 6) cigzip encode (max-complexity=32, threads=8) --minimal → cigzip-datasets/dataset.standard.21007.tp.mc32.paf

**Command**:
```bash
/usr/bin/time ../cigzip/target/debug/cigzip encode -p data/simulated/cigzip-datasets/dataset.standard.21007.paf --type standard --complexity-metric edit-distance --max-complexity 32 -t 8 --minimal > data/simulated/cigzip-datasets/dataset.standard.21007.tp.mc32.paf
```
**Output**:
```text
2.06user 0.18system 0:00.34elapsed 660%CPU (0avgtext+0avgdata 9064maxresident)k
0inputs+10440outputs (0major+600minor)pagefaults 0swaps
```

### 7) cigzip decode + diff verification

**Command**:
```bash
/usr/bin/time -v ../cigzip/target/debug/cigzip decode -p data/simulated/cigzip-datasets/dataset.standard.21007.tp.mc32.paf --type standard --complexity-metric edit-distance --sequence-files data/simulated/wf-dataset/dataset.21007.fa --max-complexity 32 -t 8 > data/simulated/cigzip-datasets/dataset.standard.21007.tp.mc32.decompressed.paf && diff <(sort data/simulated/cigzip-datasets/dataset.standard.21007.paf) <(sort data/simulated/cigzip-datasets/dataset.standard.21007.tp.mc32.decompressed.paf) | wc -l
```
**Output**:
```text
[2026-01-15T17:49:13Z WARN  cigzip] Ignoring --max-complexity=32 because --heuristic was not requested
	Command being timed: "../cigzip/target/debug/cigzip decode -p data/simulated/cigzip-datasets/dataset.standard.21007.tp.mc32.paf --type standard --complexity-metric edit-distance --sequence-files data/simulated/wf-dataset/dataset.21007.fa --max-complexity 32 -t 8"
	User time (seconds): 590.30
	System time (seconds): 16.13
	Percent of CPU this job got: 794%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:16.31
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 40944
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 4
	Minor (reclaiming a frame) page faults: 60719
	Voluntary context switches: 3738
	Involuntary context switches: 915
	Swaps: 0
	File system inputs: 0
	File system outputs: 35560
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
| data/simulated/wf-dataset/sample.dataset.21007.seq | 203507741 |
| data/simulated/wf-dataset/dataset.21007.fa | 203735529 |
| data/simulated/fastga-datasets/dataset.21007.paf | 12467907 |
| data/simulated/tmp/tmp.standard.21007.paf | 5781370 |
| data/simulated/cigzip-datasets/dataset.standard.21007.paf | 18138087 |
| cigzip-datasets/dataset.standard.21007.tp.mc32.paf | 5344676 |
| cigzip-datasets/dataset.standard.21007.tp.mc32.decompressed.paf | 18138087 |

## Artifacts MB


| data/simulated/wf-dataset/sample.dataset.21007.seq | data/simulated/wf-dataset/dataset.21007.fa | data/simulated/fastga-datasets/dataset.21007.paf | data/simulated/tmp/tmp.standard.21007.paf | data/simulated/cigzip-datasets/dataset.standard.21007.paf | cigzip-datasets/dataset.standard.21007.tp.mc32.paf | cigzip-datasets/dataset.standard.21007.tp.mc32.decompressed.paf |
|---|---|---|---|---|---|---|
| 194.08 | 194.30 | 11.89 | 5.51 | 17.30 | 5.10 | 17.30 |
