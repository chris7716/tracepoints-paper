#!/usr/bin/env bash
# validate_fastga_tracepoints.sh
#
# For every simulated dataset, compare FASTGA PAFtoALN tracepoints against
# cigzip-encoded tracepoints from the original WFA2 PAF.
#
# Modes (default: both):
#   --fastga-no-diff   compare FASTGA T vs cigzip fastga-no-diff lengths
#                      (uses pre-existing encode/ files)
#   --fastga           compare FASTGA T/X vs cigzip fastga lengths/diffs
#                      (re-encodes with --type fastga on the fly)
#
# Usage:
#   ./validate_fastga_tracepoints.sh [--fastga-no-diff] [--fastga]
#
# All paths can be overridden via environment variables (see CONFIGURATION).

set -euo pipefail

# =============================================================================
# ARGS
# =============================================================================

check_no_diff=0
check_fastga=0

for arg in "$@"; do
    case "$arg" in
        --fastga-no-diff) check_no_diff=1 ;;
        --fastga)         check_fastga=1  ;;
        *) echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done

# default: run both
if [[ $check_no_diff -eq 0 && $check_fastga -eq 0 ]]; then
    check_no_diff=1
    check_fastga=1
fi

# =============================================================================
# CONFIGURATION
# =============================================================================

dir_base="${DIR_BASE:-/home/hasitha/data/projects/tracepoints-paper/output/}"
cigzip="${CIGZIP:-/home/hasitha/data/projects/cigzip/target/release/cigzip}"
fastga_dir="${FASTGA_DIR:-/home/hasitha/data/projects/FASTGA}"
trace_spacing=100

# suffix used by the benchmark script for fastga-no-diff encode outputs
fgand_suffix="fastga-no-diff.none.100.high"

paftoaln="$fastga_dir/PAFtoALN"
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

for tool_var in paftoaln oneview; do
    path="${!tool_var}"
    [[ -x "$path" ]] || die "$tool_var not found or not executable: $path"
done

if [[ $check_fastga -eq 1 ]]; then
    [[ -x "$cigzip" ]] || die "cigzip not found or not executable: $cigzip"
fi

log "Configuration:"
log "  DIR_BASE           = $dir_base"
log "  CIGZIP             = $cigzip"
log "  FASTGA_DIR         = $fastga_dir"
log "  trace_spacing      = $trace_spacing"
log "  check fastga-no-diff = $check_no_diff  (uses pre-encoded files)"
log "  check fastga         = $check_fastga   (re-encodes on the fly)"
log ""

# =============================================================================
# COMPARISON (Python, inlined)
#
#   $1  cigzip --type fastga PAF   (or empty string to skip fastga checks)
#   $2  cigzip --type fastga-no-diff PAF  (or empty string to skip)
#   $3  ONEview text of the .1aln
#   $4  prefix label
# =============================================================================

COMPARE_PY=$(cat <<'PYEOF'
import sys

fga_paf    = sys.argv[1]   # empty string  → skip fastga checks
fga_nd_paf = sys.argv[2]   # empty string  → skip fastga-no-diff checks
oneview    = sys.argv[3]
label      = sys.argv[4]

# ── PAF helpers ──────────────────────────────────────────────────────────────

def paf_tp(line):
    for f in line.split('\t'):
        if f.startswith('tp:Z:'):
            return f[5:].rstrip()
    return None

def parse_fastga(tp_str):
    pairs = []
    for seg in tp_str.split(';'):
        seg = seg.strip()
        if seg:
            d, l = seg.split(',')
            pairs.append((int(d), int(l)))
    return pairs

def parse_no_diff(tp_str):
    return [int(x) for x in tp_str.split(';') if x.strip()]

def paf_key(line):
    f = line.split('\t')
    # (qname, tname, strand, qstart, qend, tstart, tend)
    return (f[0], f[5], f[4], int(f[2]), int(f[3]), int(f[7]), int(f[8]))

def load_paf(path, parser):
    if not path:
        return {}
    records = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            tp = paf_tp(line)
            if tp:
                records[paf_key(line)] = parser(tp)
    return records

fga_records    = load_paf(fga_paf,    parse_fastga)
fga_nd_records = load_paf(fga_nd_paf, parse_no_diff)

# ── ONEview parser ────────────────────────────────────────────────────────────
# A record: A aIdx aStart aLen bIdx bStart bLen
#   aStart = alignment start on a
#   aLen   = total length of a  (= alignment end; alignment always runs to end)
#   bStart = alignment start on b
#   bLen   = total length of b  (= alignment end)
# S record: S <name_len> <name>  (listed in index order)
# T record: T <count> <v1> ...   (segment lengths on b)
# X record: X <count> <v1> ...   (differences per segment)

seq_names   = []
fastga_recs = {}   # (aname, bname, aStart, aLen, bStart, bLen) → {T, X}
cur_key     = None

with open(oneview) as fh:
    for raw in fh:
        line = raw.rstrip('\n')
        if not line:
            continue
        parts = line.split()
        tag   = parts[0]

        if tag == 'S':
            seq_names.append(parts[2] if len(parts) > 2 else '')

        elif tag == 'A':
            ai, a_start, a_len, bi, b_start, b_len = (
                int(parts[1]), int(parts[2]), int(parts[3]),
                int(parts[4]), int(parts[5]), int(parts[6]))
            aname   = seq_names[ai] if ai < len(seq_names) else f'seq{ai}'
            bname   = seq_names[bi] if bi < len(seq_names) else f'seq{bi}'
            cur_key = (aname, bname, a_start, a_len, b_start, b_len)
            fastga_recs[cur_key] = {'T': [], 'X': []}

        elif tag == 'T' and cur_key:
            fastga_recs[cur_key]['T'] = list(map(int, parts[2:]))

        elif tag == 'X' and cur_key:
            fastga_recs[cur_key]['X'] = list(map(int, parts[2:]))

# ── lookup by FASTGA coordinates ─────────────────────────────────────────────
# cigzip also trims to the suffix-suffix overlap, so its PAF uses the same
# coordinates as PAFtoALN.  We use PAFtoALN's coords as the authoritative key.

def lookup_cigzip(records, aname, bname, a_start, a_len, b_start, b_len):
    k = (aname, bname, '+', a_start, a_len, b_start, b_len)
    if k in records:
        return records[k]
    k2 = (bname, aname, '+', b_start, b_len, a_start, a_len)
    if k2 in records:
        return records[k2]
    return None

# ── compare ──────────────────────────────────────────────────────────────────

n_total  = 0
n_nd_ok  = 0
n_t_ok   = 0
n_x_ok   = 0
n_cigzip_missing = 0
mismatches = []

do_fga    = bool(fga_records)
do_fga_nd = bool(fga_nd_records)

for fk, fv in fastga_recs.items():
    aname, bname, a_start, a_len, b_start, b_len = fk
    ft = fv.get('T', [])
    fx = fv.get('X', [])

    # skip zero-length alignments (PAFtoALN edge case for contained reads)
    if a_start == a_len or b_start == b_len:
        continue

    n_total += 1
    kl = f"{aname}→{bname} [{a_start}:{a_len}, {b_start}:{b_len}]"

    # 1. FASTGA T == cigzip fastga-no-diff lengths
    if do_fga_nd:
        fga_nd_tp = lookup_cigzip(fga_nd_records, aname, bname,
                                  a_start, a_len, b_start, b_len)
        if fga_nd_tp is None:
            mismatches.append(f"  cigzip fastga-no-diff record not found: {kl}")
        elif ft == fga_nd_tp:
            n_nd_ok += 1
        else:
            mismatches.append(f"  FASTGA T vs cigzip fastga-no-diff MISMATCH: {kl}")
            mismatches.append(f"    FASTGA T        : {ft[:8]}")
            mismatches.append(f"    cigzip no-diff  : {fga_nd_tp[:8]}")

    # 2. FASTGA T/X == cigzip fastga lengths/diffs
    if do_fga:
        fga_tp = lookup_cigzip(fga_records, aname, bname,
                               a_start, a_len, b_start, b_len)
        if fga_tp is None:
            n_cigzip_missing += 1
            mismatches.append(f"  cigzip fastga record not found: {kl}")
        else:
            fga_T = [l for (d, l) in fga_tp]
            fga_D = [d for (d, l) in fga_tp]

            if ft == fga_T:
                n_t_ok += 1
            else:
                mismatches.append(f"  FASTGA T vs cigzip fastga T MISMATCH: {kl}")
                mismatches.append(f"    FASTGA T  : {ft[:8]}")
                mismatches.append(f"    cigzip T  : {fga_T[:8]}")

            if fx == fga_D:
                n_x_ok += 1
            else:
                mismatches.append(f"  FASTGA X vs cigzip fastga D MISMATCH: {kl}")
                mismatches.append(f"    FASTGA X  : {fx[:8]}")
                mismatches.append(f"    cigzip D  : {fga_D[:8]}")

# ── report ────────────────────────────────────────────────────────────────────

def check(desc, n, total):
    ok = (n == total)
    print(f"  {'PASS' if ok else 'FAIL'}  {desc}: {n}/{total}")
    return ok

print(f"[{label}]  FASTGA records: {len(fastga_recs)}  compared: {n_total}"
      f"  cigzip-fastga missing: {n_cigzip_missing}")

all_ok = True
if do_fga_nd:
    all_ok &= check("FASTGA T == cigzip fastga-no-diff lengths", n_nd_ok, n_total)
if do_fga:
    all_ok &= check("FASTGA T == cigzip fastga lengths         ", n_t_ok,  n_total)
    all_ok &= check("FASTGA X == cigzip fastga diffs           ", n_x_ok,  n_total)

if mismatches:
    print(f"  Mismatch details (first 30 lines):")
    for m in mismatches[:30]:
        print(m)
    if len(mismatches) > 30:
        print(f"  ... ({len(mismatches) - 30} more lines omitted)")

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

        log "Processing: $prefix"

        for f in "$input_paf" "$seq_fa"; do
            if [[ ! -f "$f" ]]; then
                log "  SKIP – missing file: $f"
                n_skip=$((n_skip + 1))
                continue 2
            fi
        done

        job="$TMP_DIR/$prefix"
        mkdir -p "$job"

        # ── fastga-no-diff: use pre-existing encode file ──────────────────
        fga_nd_tp_paf=""
        if [[ $check_no_diff -eq 1 ]]; then
            pre="$dir_base/simulated-data/encode/$prefix.$fgand_suffix.tp.paf"
            if [[ -f "$pre" ]]; then
                fga_nd_tp_paf="$pre"
            else
                log "  WARN – fastga-no-diff encode file not found, skipping no-diff check: $pre"
            fi
        fi

        # ── fastga: encode the original WFA2 PAF on the fly ──────────────
        fga_tp_paf=""
        if [[ $check_fastga -eq 1 ]]; then
            fga_tp_paf="$job/cigzip_fastga.tp.paf"
            "$cigzip" encode \
                --paf "$input_paf" \
                --type fastga \
                --max-complexity "$trace_spacing" \
                > "$fga_tp_paf" 2>/dev/null
        fi

        # ── PAFtoALN ─────────────────────────────────────────────────────
        cp "$input_paf" "$job/input.paf"
        "$paftoaln" "$job/input.paf" "$seq_fa" 2>/dev/null || true
        aln_file="$job/input.1aln"
        if [[ ! -f "$aln_file" ]]; then
            log "  SKIP – PAFtoALN did not produce $aln_file"
            n_skip=$((n_skip + 1))
            rm -rf "$job"
            continue
        fi

        # ── ONEview ──────────────────────────────────────────────────────
        oneview_out="$job/oneview.txt"
        "$oneview" -t aln "$aln_file" 2>/dev/null > "$oneview_out"

        # ── compare ──────────────────────────────────────────────────────
        if python3 - "$fga_tp_paf" "$fga_nd_tp_paf" "$oneview_out" \
                "$prefix" <<< "$COMPARE_PY"; then
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
