# Information about the dataset

The following simulated data is generated using the [generate_dataset tool from WFA2-lib](https://github.com/smarco/WFA2-lib/tree/main/tools).

## Dataset Generation

### 1) WFA2 generate_dataset

**Command**:
```bash
generate_dataset -n <Record Count> -l <Sequence Length> -e <Error> -o wf-dataset/sample.dataset.<N>.seq
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
}' wf-dataset/sample.dataset.<N>.seq > wf-dataset/dataset.<N>.fa
```

### 3) FastGA self-alignment → fastga-datasets/dataset.N.paf

**Command**:
```bash
FastGA -pafx wf-dataset/dataset.<N>.fa wf-dataset/dataset.<N>.fa > fastga-datasets/dataset.<N>.paf
```

### 4) cigzip encode (type=standard) → tmp/tmp.standard.N.paf

**Command**:
```bash
cigzip encode -p fastga-datasets/dataset.<N>.paf --type standard --complexity-metric edit-distance --max-complexity 99999999 > tmp/tmp.standard.<N>.paf
```

### 5) cigzip decode (type=standard) → cigzip-datasets/dataset.standard.N.paf

**Command**:
```bash
cigzip decode -p tmp/tmp.standard.<N>.paf --type standard --complexity-metric edit-distance --sequence-files wf-dataset/dataset.<N>.fa > cigzip-datasets/dataset.standard.<N>.paf
```

| Dataset | Record Count | Sequence Length | Error | Tracepoint Type |
|---------|--------------|--------|-------|-----------------|
| N1 | 10000 | 100   | 0.001 | Standard |
| N2 | 10000 | 100   | 0.01  | Standard |
| N3 | 10000 | 100   | 0.05  | Standard |
| N4 | 10000 | 100   | 0.1   | Standard |
| N5 | 10000 | 1000  | 0.001 | Standard |
| N6 | 10000 | 1000  | 0.01  | Standard |
| N7 | 10000 | 1000  | 0.05  | Standard |
| N8 | 10000 | 1000  | 0.1   | Standard |
| N9 | 10000 | 10000 | 0.001 | Standard |
| N10 | 10000 | 10000 | 0.01  | Standard |
| N11 | 10000 | 10000 | 0.05  | Standard |
| N12 | 10000 | 10000 | 0.1   | Standard |