#!/usr/bin/env python3
"""Print Supplementary Table S3 (tracepoints per alignment) as a LaTeX tabular."""
import csv
import os
import sys

LENGTHS = [(100, "100 bp"), (1000, "1 Kbp"), (10000, "10 Kbp"), (100000, "100 Kbp")]
METHODS = ["FL-TP", "EB-TP", "DB-TP"]
ERRORS = [0.001, 0.01, 0.05, 0.1, 0.2]


def log(msg):
    print(f"[make_table_s3] {msg}", file=sys.stderr)


def method_of(r):
    # FL-TP l=100; EB/DB threshold 32
    if r["tp_type"] == "fastga" and r["mc"] == "100":
        return "FL-TP"
    if r["tp_type"] == "standard" and r["mc"] == "32":
        return {"edit-distance": "EB-TP", "diagonal-distance": "DB-TP"}.get(r["cm"])
    return None


def load_means(path):
    means, rows = {}, 0
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows += 1
            m = method_of(r)
            if r["memory_mode"] != "high" or m is None:
                continue
            key = (int(r["l"]), m, float(r["e"]))
            if key in means:
                sys.exit(f"[make_table_s3] ERROR: duplicate row for {key}")
            n_aln = int(r["num_input_alignments"])
            if n_aln <= 0:
                sys.exit(f"[make_table_s3] ERROR: {n_aln} alignments for {key}")
            means[key] = (int(r["num_tracepoints"]) / n_aln, n_aln)
    log(f"read {rows} rows from {path}, kept {len(means)}")
    return means


def main():
    data_dir = os.environ.get("TRACEPOINTS_DATA_DIR", "data")
    means = load_means(os.path.join(data_dir, "simulated-data", "benchmark.results.tsv"))
    missing = [(l, m, e) for l, _ in LENGTHS for m in METHODS for e in ERRORS if (l, m, e) not in means]
    if missing:
        sys.exit(f"[make_table_s3] ERROR: missing cells {missing}")
    log("alignments per condition: " + ", ".join(sorted({str(n) for _, n in means.values()})))
    fmt = lambda x: f"{x:,.2f}".replace(",", "{,}")
    out = ["\\begin{tabular}{llrrrrr}", "\\hline",
           "\\textbf{Length} & \\textbf{Method} & " + " & ".join(f"\\textbf{{{e * 100:g}\\%}}" for e in ERRORS) + " \\\\"]
    for l, label in LENGTHS:
        out.append("\\hline")
        for j, m in enumerate(METHODS):
            first = f"\\multirow{{3}}{{*}}{{{label}}}" if j == 0 else ""
            out.append(f"{first:24s} & {m} & " + " & ".join(fmt(means[(l, m, e)][0]) for e in ERRORS) + " \\\\")
    out += ["\\hline", "\\end{tabular}"]
    print("\n".join(out))
    log(f"printed {len(LENGTHS) * len(METHODS) * len(ERRORS)} cells")


if __name__ == "__main__":
    main()
