#!/usr/bin/env bash
# validate_fastga_tracepoints.sh
#
# For every simulated dataset, runs FASTGA's PAFtoALN on the original alignment
# PAF and compares its tracepoints against:
#   1. cigzip --type fastga           (FL-TP with diffs)    – encoded on-the-fly
#   2. cigzip --type fastga-no-diff   (FL-TP without diffs) – uses pre-existing
#                                                             files from encode/
#
# Usage:
#   ./validate_fastga_tracepoints.sh
#
# All paths can be overridden via environment variables (see CONFIGURATION).

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

dir_base="${DIR_BASE:-/home/hasitha/data/projects/tracepoints-paper/output/}"
cigzip="${CIGZIP:-/home/hasitha/data/projects/cigzip/target/release/cigzip}"
fastga_dir="${FASTGA_DIR:-/home/hasitha/data/projects/FASTGA}"
trace_spacing=100

paftoaln="$fastga_dir/PAFtoALN"
alnvtopaf="$fastga_dir/ALNtoPAF"
oneview="$fastga_dir/ONEview"

# =============================================================================
# HELPERS
# =============================================================================

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

die() { echo "ERROR: $*" >&2; exit 1; }

n_pass=0
n_fail=0
n_skip=0
failed_prefixes=()

# =============================================================================
# TOOL CHECKS
# =============================================================================

for tool_var in cigzip paftoaln alnvtopaf oneview; do
    path="${!tool_var}"
    [[ -x "$path" ]] || die "$tool_var not found or not executable: $path"
done

log "Configuration:"
log "  DIR_BASE    = $dir_base"
log "  CIGZIP      = $cigzip"
log "  FASTGA_DIR  = $fastga_dir"
log "  trace_spacing = $trace_spacing"
log ""

# =============================================================================
# COMPARISON (Python, inlined)
# =============================================================================

# Accepts three positional args written to a temp file by the caller:
#   $1  fastga tp PAF (cigzip --type fastga)
#   $2  fastga-no-diff tp PAF (cigzip --type fastga-no-diff)
#   $3  ONEview text output
#   $4  trace spacing (int)
#   $5  prefix label (for messages)
COMPARE_PY=$(cat <<'PYEOF'
import sys, re, collections

fga_paf    = sys.argv[1]
fga_nd_paf = sys.argv[2]
oneview    = sys.argv[3]
spacing    = int(sys.argv[4])
label      = sys.argv[5]

def paf_tp_field(line):
    for f in line.split('\t'):
        if f.startswith('tp:Z:'):
            return f[5:].rstrip()
    return None

def parse_fastga_tp(tp_str):
    pairs = []
    for seg in tp_str.split(';'):
        seg = seg.strip()
        if not seg:
            continue
        d, l = seg.split(',')
        pairs.append((int(d), int(l)))
    return pairs

def parse_fastga_no_diff_tp(tp_str):
    return [int(x) for x in tp_str.split(';') if x.strip()]

def key_from_paf(line):
    f = line.split('\t')
    return (f[0], f[5], f[4], int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[6]))

# -- load cigzip fastga tracepoints --
fga_records    = {}
fga_nd_records = {}

with open(fga_paf) as fh:
    for line in fh:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        k  = key_from_paf(line)
        tp = paf_tp_field(line)
        if tp:
            fga_records[k] = parse_fastga_tp(tp)

with open(fga_nd_paf) as fh:
    for line in fh:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        k  = key_from_paf(line)
        tp = paf_tp_field(line)
        if tp:
            fga_nd_records[k] = parse_fastga_no_diff_tp(tp)

# -- load FASTGA tracepoints from ONEview --
scaff_groups  = [[]]
fastga_records = {}
cur_a = None

with open(oneview) as fh:
    for raw in fh:
        line = raw.rstrip('\n')
        if not line:
            continue
        tag = line.split()[0]

        if tag == 'g':
            if scaff_groups[-1]:
                scaff_groups.append([])

        if tag == 'S':
            parts = line.split(None, 2)
            scaff_groups[-1].append(parts[2] if len(parts) > 2 else '')

        elif tag == 'A':
            parts = line.split()
            ai, as_, ae, bi, bs, be = (int(parts[1]), int(parts[2]), int(parts[3]),
                                        int(parts[4]), int(parts[5]), int(parts[6]))
            a_group = scaff_groups[0]
            b_group = scaff_groups[1] if len(scaff_groups) > 1 else scaff_groups[0]
            aname = a_group[ai] if ai < len(a_group) else f'scaff{ai}'
            bname = b_group[bi] if bi < len(b_group) else f'scaff{bi}'
            cur_a = (aname, bname, as_, ae, bs, be)
            fastga_records[cur_a] = {'T': [], 'X': []}

        elif tag == 'R' and cur_a:
            fastga_records[cur_a]['comp'] = True

        elif tag == 'T' and cur_a:
            vals = list(map(int, line.split()[2:]))
            fastga_records[cur_a]['T'] = vals

        elif tag == 'X' and cur_a:
            vals = list(map(int, line.split()[2:]))
            fastga_records[cur_a]['X'] = vals

# -- matching helper --
def find_fastga(qname, tname, strand, qs, qe, ts, te, tlen):
    k = (qname, tname, qs, qe, ts, te)
    if k in fastga_records:
        return fastga_records[k]
    if strand == '-':
        k_rev = (qname, tname, qs, qe, tlen - te, tlen - ts)
        if k_rev in fastga_records:
            return fastga_records[k_rev]
    k2 = (tname, qname, ts, te, qs, qe)
    if k2 in fastga_records:
        return fastga_records[k2]
    if strand == '-':
        k2_rev = (tname, qname, tlen - te, tlen - ts, qs, qe)
        if k2_rev in fastga_records:
            return fastga_records[k2_rev]
    return None

# -- comparisons --
n_records       = 0
n_match_nd      = 0   # fastga T lens == fastga-no-diff lens
n_match_fga_t   = 0   # fastga T == FASTGA T
n_match_fga_x   = 0   # fastga diffs == FASTGA X
mismatches      = []

for k, fga_tp in fga_records.items():
    qname, tname, strand, qs, qe, ts, te, tlen = k
    n_records += 1

    fga_T   = [l for (d, l) in fga_tp]
    fga_D   = [d for (d, l) in fga_tp]
    kl      = f"{qname}→{tname} [{qs}:{qe}, {ts}:{te}]"

    # 1. fastga T lengths must match fastga-no-diff values
    fga_nd = fga_nd_records.get(k)
    if fga_nd is not None:
        if fga_T == fga_nd:
            n_match_nd += 1
        else:
            mismatches.append(f"  fastga-len vs fastga-no-diff MISMATCH: {kl}")
            mismatches.append(f"    fastga T      : {fga_T[:10]}")
            mismatches.append(f"    fastga-no-diff: {fga_nd[:10]}")
    else:
        mismatches.append(f"  fastga-no-diff record not found: {kl}")

    # 2 & 3. fastga T/X must match FASTGA T/X
    fastga = find_fastga(qname, tname, strand, qs, qe, ts, te, tlen)
    if fastga is not None:
        ft = fastga.get('T', [])
        fx = fastga.get('X', [])
        if fga_T == ft:
            n_match_fga_t += 1
        else:
            mismatches.append(f"  fastga T vs FASTGA T MISMATCH: {kl}")
            mismatches.append(f"    cigzip T : {fga_T[:10]}")
            mismatches.append(f"    FASTGA T : {ft[:10]}")
        if fga_D == fx:
            n_match_fga_x += 1
        else:
            mismatches.append(f"  fastga D vs FASTGA X MISMATCH: {kl}")
            mismatches.append(f"    cigzip D : {fga_D[:10]}")
            mismatches.append(f"    FASTGA X : {fx[:10]}")
    else:
        mismatches.append(f"  FASTGA alignment not found: {kl}")
        candidates = [(k2, v) for k2, v in fastga_records.items()
                      if (k2[0]==qname and k2[1]==tname) or (k2[0]==tname and k2[1]==qname)]
        for ck, cv in candidates[:3]:
            mismatches.append(f"    FASTGA has: {ck[0]}→{ck[1]} [{ck[2]}:{ck[3]}, {ck[4]}:{ck[5]}]  T={cv.get('T',[])[:5]}")

# -- report --
def result(ok, desc, n, total):
    status = "PASS" if ok else "FAIL"
    print(f"  {status}  {desc}: {n}/{total}")
    return ok

print(f"[{label}]  records: {n_records}  FASTGA alignments: {len(fastga_records)}")
all_ok = True
all_ok &= result(n_match_nd    == n_records, "cigzip fastga lengths == fastga-no-diff lengths (boundary consistency)", n_match_nd,    n_records)
all_ok &= result(n_match_fga_t == n_records, "cigzip fastga lengths == FASTGA T row (boundary vs reference)         ", n_match_fga_t, n_records)
all_ok &= result(n_match_fga_x == n_records, "cigzip fastga diffs   == FASTGA X row (diff counts vs reference)      ", n_match_fga_x, n_records)

if mismatches:
    print("  Mismatch details (first 30 lines):")
    for m in mismatches[:30]:
        print(m)
    if len(mismatches) > 30:
        print(f"  ... ({len(mismatches)-30} more lines omitted)")

sys.exit(0 if all_ok else 1)
PYEOF
)

# =============================================================================
# MAIN LOOP
# =============================================================================

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

for l in 100 1000 10000 100000; do
    for e in 0.001 0.01 0.05 0.10 0.20; do
        l_nodot=$(echo "$l" | sed 's/\.//g')
        e_nodot=$(echo "$e" | sed 's/\.//g')
        prefix="set_${l_nodot}_${e_nodot}"

        input_paf="$dir_base/simulated-data/alns/$prefix.paf"
        seq_fa="$dir_base/simulated-data/seqs/$prefix.fa"
        fga_nd_paf="$dir_base/simulated-data/encode/$prefix.fastga-no-diff.none.100.high.tp.paf"

        log "Processing: $prefix"

        # Check required input files exist
        missing=0
        for f in "$input_paf" "$seq_fa" "$fga_nd_paf"; do
            if [[ ! -f "$f" ]]; then
                log "  SKIP – missing file: $f"
                missing=1
            fi
        done
        if [[ "$missing" -eq 1 ]]; then
            n_skip=$((n_skip + 1))
            continue
        fi

        job="$TMP_DIR/$prefix"
        mkdir -p "$job"

        # -- Encode cigzip fastga (FL-TP with diffs) --
        fga_tp_paf="$job/fastga.tp.paf"
        "$cigzip" encode \
            --paf "$input_paf" \
            --type fastga \
            --max-complexity "$trace_spacing" \
            > "$fga_tp_paf" 2>/dev/null

        # -- Run PAFtoALN --
        cp "$input_paf" "$job/input.paf"
        "$paftoaln" "$job/input.paf" "$seq_fa" 2>/dev/null
        aln_file="$job/input.1aln"
        if [[ ! -f "$aln_file" ]]; then
            log "  SKIP – PAFtoALN did not produce $aln_file"
            n_skip=$((n_skip + 1))
            rm -rf "$job"
            continue
        fi

        # -- Extract FASTGA tracepoints via ONEview --
        oneview_out="$job/oneview.txt"
        "$oneview" -t aln "$aln_file" 2>/dev/null > "$oneview_out"

        # -- Compare --
        if python3 - "$fga_tp_paf" "$fga_nd_paf" "$oneview_out" \
                "$trace_spacing" "$prefix" <<< "$COMPARE_PY" ; then
            n_pass=$((n_pass + 1))
        else
            n_fail=$((n_fail + 1))
            failed_prefixes+=("$prefix")
        fi

        rm -rf "$job"
        log ""
    done
done

# =============================================================================
# SUMMARY
# =============================================================================

log "════════════════════════════════════════════════════════════"
log "Validation complete"
log "  Passed : $n_pass"
log "  Failed : $n_fail"
log "  Skipped: $n_skip"

if [[ "${#failed_prefixes[@]}" -gt 0 ]]; then
    log "  Failed prefixes:"
    for p in "${failed_prefixes[@]}"; do
        log "    $p"
    done
    log "════════════════════════════════════════════════════════════"
    exit 1
fi

log "════════════════════════════════════════════════════════════"
