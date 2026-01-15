# cigzip / FastGA run

- Timestamp: `2026-01-15T06:15:55-06:00`

## Inputs

| Parameter | Value |
|---|---|
| record_count (-r) | 10000 |
| length (-L) | 10000 |
| error (-E) | 0.1 |
| N (-N) | 20001 |
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
../WFA2-lib/bin/generate_dataset1 -n 10000 -l 10000 -e 0.1 -o data/simulated/wf-dataset/sample.dataset.20001.seq --balanced-indels 20,5 --diagonal-bandwidth 0.04
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
}' data/simulated/wf-dataset/sample.dataset.20001.seq > data/simulated/wf-dataset/dataset.20001.fa
```
**Output**:
```text
FASTA written: data/simulated/wf-dataset/dataset.20001.fa
```

### 3) FastGA self-alignment → data/simulated/fastga-datasets/dataset.20001.paf

**Command**:
```bash
../FASTGA/FastGA -P/home/hasitha/data/projects/fastga-tmp -pafx data/simulated/wf-dataset/dataset.20001.fa data/simulated/wf-dataset/dataset.20001.fa > data/simulated/fastga-datasets/dataset.20001.paf
```
**Output**:
```text
sh: 1: FAtoGDB: not found

FastGA: Call to FAtoGDB failed
