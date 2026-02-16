#!/usr/bin/env python3
"""CIGAR dot-plot visualization with tracepoint sampling grids.

Reads a PAF file with cg:Z: CIGAR tags, calls ``cigzip encode`` to compute
tracepoint boundaries using EditDistance, DiagonalDistance, and FastGA modes,
and emits dot-plot images showing the alignment path overlaid with the
tracepoint grid.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class CigarOp:
    length: int
    op: str


@dataclass
class PafAlignment:
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    cigar_ops: List[CigarOp]


@dataclass
class DotPlotSegment:
    t_start: int
    q_start: int
    t_end: int
    q_end: int
    op_type: str  # '=', 'M', 'X', 'I', 'D'


@dataclass
class TracepointBoundary:
    query_pos: int
    target_pos: int
    a_len: int  # query length of this segment
    b_len: int  # target length of this segment


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_cigar(cigar_str: str) -> List[CigarOp]:
    """Parse a CIGAR string into a list of CigarOp."""
    ops: List[CigarOp] = []
    num: List[str] = []
    for ch in cigar_str:
        if ch.isdigit():
            num.append(ch)
        else:
            if num:
                ops.append(CigarOp(int("".join(num)), ch))
            num.clear()
    return ops


def parse_paf_line(line: str) -> Optional[PafAlignment]:
    """Parse a single PAF line. Returns None when there is no cg:Z: tag."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None

    cigar_str: Optional[str] = None
    for f in fields[12:]:
        if f.startswith("cg:Z:"):
            cigar_str = f[5:]
            break
    if cigar_str is None:
        return None

    return PafAlignment(
        query_name=fields[0],
        query_length=int(fields[1]),
        query_start=int(fields[2]),
        query_end=int(fields[3]),
        strand=fields[4],
        target_name=fields[5],
        target_length=int(fields[6]),
        target_start=int(fields[7]),
        target_end=int(fields[8]),
        cigar_ops=parse_cigar(cigar_str),
    )


# ---------------------------------------------------------------------------
# Dot-plot segments from CIGAR
# ---------------------------------------------------------------------------


def cigar_to_dotplot_segments(
    ops: List[CigarOp],
    target_start: int,
    query_start: int,
) -> List[DotPlotSegment]:
    """Convert CIGAR ops to line segments for the dot-plot path."""
    segments: List[DotPlotSegment] = []
    t = target_start
    q = query_start
    for cigar_op in ops:
        length = cigar_op.length
        op = cigar_op.op
        if op in ("=", "M", "X"):
            segments.append(DotPlotSegment(t, q, t + length, q + length, op))
            t += length
            q += length
        elif op == "I":
            segments.append(DotPlotSegment(t, q, t, q + length, "I"))
            q += length
        elif op == "D":
            segments.append(DotPlotSegment(t, q, t + length, q, "D"))
            t += length
    return segments


# ---------------------------------------------------------------------------
# Tracepoint boundaries via cigzip encode
# ---------------------------------------------------------------------------


def _parse_tp_field(tp_str: str) -> List[Tuple[int, int]]:
    """Parse a tp:Z: value ('a,b;a,b;...') into (a_len, b_len) pairs."""
    pairs: List[Tuple[int, int]] = []
    for token in tp_str.split(";"):
        a_str, b_str = token.split(",")
        pairs.append((int(a_str), int(b_str)))
    return pairs


def _parse_fastga_tp(
    tp_str: str, query_consumed: int, fastga_spacing: int
) -> List[Tuple[int, int]]:
    """Parse FastGA tp:Z: (diff,b_len pairs) into (a_len, b_len) pairs.

    FastGA encodes ``diff,b_len`` per segment.  ``a_len`` is the trace
    spacing for every segment except the last, which absorbs the remainder.
    """
    raw = [(int(a), int(b)) for a, b in (t.split(",") for t in tp_str.split(";"))]
    n = len(raw)
    pairs: List[Tuple[int, int]] = []
    for i, (_diff, b_len) in enumerate(raw):
        if i < n - 1:
            a_len = fastga_spacing
        else:
            a_len = query_consumed - (n - 1) * fastga_spacing
        pairs.append((a_len, b_len))
    return pairs


def _fastga_boundaries_from_cigar(
    cigar_ops: List[CigarOp],
    query_start: int,
    target_start: int,
    spacing: int,
) -> List[TracepointBoundary]:
    """Compute FastGA tracepoint boundaries by walking the CIGAR.

    Boundaries are placed at fixed *spacing* intervals on the query.
    This avoids relying on cigzip's coordinate-dependent clipping.
    """
    total_q = sum(op.length for op in cigar_ops if op.op in ("=", "M", "X", "I"))
    # Target positions at each query offset via CIGAR walk
    # Build list of (delta_q, delta_t) micro-steps
    steps: List[tuple] = []
    for op in cigar_ops:
        if op.op in ("=", "M", "X"):
            steps.append((op.length, op.length))
        elif op.op == "I":
            steps.append((op.length, 0))
        elif op.op == "D":
            steps.append((0, op.length))

    # Walk and emit boundaries at spacing multiples
    boundaries: List[TracepointBoundary] = []
    q_acc = 0  # query consumed so far
    t_acc = 0  # target consumed so far
    prev_q = 0
    prev_t = 0
    next_boundary = spacing

    for dq, dt in steps:
        q_rem = dq
        t_rem = dt
        while q_rem > 0 and next_boundary <= total_q:
            need = next_boundary - q_acc
            if need <= q_rem:
                # proportion of target consumed
                if dq > 0:
                    t_frac = int(dt * need / dq) if q_rem == dq else need
                    # for equal-consuming ops (=, X, M) dt==dq so t_frac==need
                    # for insertions dt==0 so t_frac==0
                    if dq == dt:
                        t_frac = need
                    elif dt == 0:
                        t_frac = 0
                else:
                    t_frac = 0
                q_acc += need
                t_acc += t_frac
                q_rem -= need
                t_rem -= t_frac
                a_len = q_acc - prev_q
                b_len = t_acc - prev_t
                boundaries.append(TracepointBoundary(
                    query_start + q_acc, target_start + t_acc, a_len, b_len,
                ))
                prev_q = q_acc
                prev_t = t_acc
                next_boundary += spacing
            else:
                break
        q_acc += q_rem
        t_acc += t_rem

    # Final segment (remainder)
    if q_acc > prev_q or t_acc > prev_t:
        boundaries.append(TracepointBoundary(
            query_start + q_acc, target_start + t_acc,
            q_acc - prev_q, t_acc - prev_t,
        ))
    return boundaries


def _pairs_to_boundaries(
    pairs: List[Tuple[int, int]],
    query_start: int,
    target_start: int,
) -> List[TracepointBoundary]:
    """Convert (a_len, b_len) pairs to cumulative boundary positions."""
    boundaries: List[TracepointBoundary] = []
    q_pos = query_start
    t_pos = target_start
    for a_len, b_len in pairs:
        q_pos += a_len
        t_pos += b_len
        boundaries.append(TracepointBoundary(q_pos, t_pos, a_len, b_len))
    return boundaries


@dataclass
class CigzipResult:
    """Result for one alignment from ``cigzip encode``."""
    pairs: List[Tuple[int, int]]
    query_start: int
    target_start: int


def run_cigzip_encode(
    paf_path: str,
    max_complexity: int,
    cigzip_bin: str,
    *,
    tracepoint_type: str = "standard",
    metric: Optional[str] = None,
    fastga_spacing: Optional[int] = None,
) -> List[CigzipResult]:
    """Run ``cigzip encode`` and return tracepoint results per alignment.

    For *tracepoint_type* ``"standard"``, *metric* must be
    ``"edit-distance"`` or ``"diagonal-distance"``.  For ``"fastga"``,
    *fastga_spacing* is required.

    Returns a list of :class:`CigzipResult` (one per input alignment).
    """
    if tracepoint_type == "standard":
        cmd = [
            cigzip_bin,
            "encode",
            "--paf",
            paf_path,
            "--max-complexity",
            str(max_complexity),
            "--complexity-metric",
            metric,
            "--minimal",
        ]
    else:
        cmd = [
            cigzip_bin,
            "encode",
            "--paf",
            paf_path,
            "--type",
            "fastga",
            "--max-complexity",
            str(max_complexity),
            "--minimal",
        ]

    label = metric if tracepoint_type == "standard" else "fastga"
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"cigzip encode failed ({label}):", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(1)

    results: List[CigzipResult] = []
    for line in result.stdout.splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue

        tp_str: Optional[str] = None
        for field in fields[12:]:
            if field.startswith("tp:Z:"):
                tp_str = field[5:]
                break
        if tp_str is None:
            continue

        q_start = int(fields[2])
        q_end = int(fields[3])
        t_start = int(fields[7])

        if tracepoint_type == "fastga":
            query_consumed = q_end - q_start
            pairs = _parse_fastga_tp(tp_str, query_consumed, fastga_spacing)
        else:
            pairs = _parse_tp_field(tp_str)

        results.append(CigzipResult(pairs=pairs, query_start=q_start, target_start=t_start))
    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

OP_COLORS = {
    "=": "#2ca02c",  # green – matches
    "M": "#2ca02c",
    "X": "#d62728",  # red – mismatches
    "I": "#9467bd",  # purple – insertions
    "D": "#ff7f0e",  # orange – deletions
}

OP_LABELS = {"=": "matches", "X": "mismatches", "I": "insertions", "D": "deletions"}


def _render_cigar_annotation(ax, fig, cigar_ops, font_scale=1.0):
    """Render colored CIGAR string centered at the top of the figure."""
    from matplotlib.offsetbox import HPacker, TextArea, AnchoredOffsetbox

    fs = 30 * font_scale
    fp = dict(fontsize=fs, fontweight="bold", fontfamily="monospace")

    boxes = [TextArea("CIGAR: ", textprops={**fp, "color": "black"})]
    for op in cigar_ops:
        boxes.append(
            TextArea(
                f"{op.length}{op.op}",
                textprops={**fp, "color": OP_COLORS.get(op.op, "gray")},
            )
        )

    hpack = HPacker(children=boxes, align="baseline", pad=0, sep=0)
    ab = AnchoredOffsetbox(
        loc="upper center",
        child=hpack,
        pad=0,
        frameon=False,
        bbox_to_anchor=(0.5, 0.99),
        bbox_transform=fig.transFigure,
    )
    ax.add_artist(ab)


def _render_dotplot_on_ax(
    ax,
    alignment: PafAlignment,
    segments: List[DotPlotSegment],
    boundaries: List[TracepointBoundary],
    metric_label: str,
    show_labels: bool = False,
    font_scale: float = 1.0,
) -> None:
    """Draw the alignment path and tracepoint grid onto *ax*."""
    aln_span = max(
        alignment.target_end - alignment.target_start,
        alignment.query_end - alignment.query_start,
        1,
    )
    path_lw = max(2.5, min(8.0, 3000 / aln_span)) * font_scale

    # -- Alignment path --
    lines_by_color: dict[str, list] = {}
    for seg in segments:
        color = OP_COLORS.get(seg.op_type, "gray")
        lines_by_color.setdefault(color, []).append(
            [(seg.t_start, seg.q_start), (seg.t_end, seg.q_end)]
        )
    for color, lines in lines_by_color.items():
        lc = LineCollection(lines, colors=color, linewidths=path_lw)
        ax.add_collection(lc)

    # -- Tracepoint grid --
    t_bounds = sorted({alignment.target_start} | {b.target_pos for b in boundaries})
    q_bounds = sorted({alignment.query_start} | {b.query_pos for b in boundaries})

    grid_lw = max(1.2, min(4.0, 1600 / aln_span)) * font_scale
    t_start, t_end = alignment.target_start, alignment.target_end
    q_start, q_end = alignment.query_start, alignment.query_end
    # Dashed interior grid lines (stop at the alignment path)
    pad = aln_span * 0.015
    for b in boundaries:
        if b.target_pos != t_start and b.target_pos != t_end:
            ax.plot([b.target_pos, b.target_pos], [q_start - pad, b.query_pos],
                    color="black", linewidth=grid_lw, linestyle="--", alpha=0.45)
        if b.query_pos != q_start and b.query_pos != q_end:
            ax.plot([t_start - pad, b.target_pos], [b.query_pos, b.query_pos],
                    color="black", linewidth=grid_lw, linestyle="--", alpha=0.45)

    # -- Dots at tracepoint boundary positions (on the alignment path) --
    dot_size = max(6, min(16, int(6000 / aln_span))) * font_scale
    dot_t = [alignment.target_start] + [b.target_pos for b in boundaries]
    dot_q = [alignment.query_start] + [b.query_pos for b in boundaries]
    ax.plot(dot_t, dot_q, 'o', color='black', markersize=dot_size, alpha=0.6, zorder=5)

    # -- Cell labels (near boundary point, offset into open quadrant) --
    if show_labels:
        label_fs = max(20, min(36, int(9000 / aln_span))) * font_scale
        label_offset = aln_span * 0.02
        edge_margin = aln_span * 0.30
        for b in boundaries:
            # Flip placement when near the right or bottom edge
            near_right = (t_end - b.target_pos) < edge_margin
            near_bottom = (b.query_pos - q_start) < edge_margin
            if near_right:
                cx = b.target_pos - label_offset * 4
                ha = "right"
            elif near_bottom:
                cx = b.target_pos + label_offset * 3
                ha = "left"
            else:
                cx = b.target_pos + label_offset
                ha = "left"
            if near_bottom:
                cy = b.query_pos + label_offset
                va = "bottom"
            else:
                cy = b.query_pos - label_offset
                va = "top"
            ax.text(
                cx,
                cy,
                f"({b.a_len},{b.b_len})",
                fontsize=label_fs,
                fontweight="bold",
                ha=ha,
                va=va,
                color="black",
                alpha=0.85,
            )

    # -- Axes & title --
    ax.set_xlim(alignment.target_start - pad, alignment.target_end + pad)
    ax.set_ylim(alignment.query_start - pad, alignment.query_end + pad)
    ax.set_xlabel("Target position (bp)", fontsize=30 * font_scale)
    ax.set_ylabel("Query position (bp)", fontsize=30 * font_scale)
    ax.tick_params(labelsize=22 * font_scale, length=8 * font_scale, width=1.5 * font_scale)
    # Clean tick positions at round multiples
    import numpy as np
    from matplotlib.ticker import MultipleLocator
    span = max(alignment.target_end - alignment.target_start,
               alignment.query_end - alignment.query_start)
    # Pick the largest "nice" step giving 3-8 ticks
    for step in [500, 200, 100, 50, 20, 10]:
        if span / step >= 3:
            break
    ax.xaxis.set_major_locator(MultipleLocator(step))
    ax.yaxis.set_major_locator(MultipleLocator(step))
    ax.set_title(metric_label, fontsize=30 * font_scale, pad=35 * font_scale)
    ax.text(
        0.5, 1.01,
        f"num. tracepoints={len(boundaries)}",
        transform=ax.transAxes,
        fontsize=24 * font_scale,
        ha="center", va="bottom",
    )
    ax.set_aspect("equal")
    ax.autoscale_view()


def plot_dotplot(
    alignment: PafAlignment,
    segments: List[DotPlotSegment],
    boundaries: List[TracepointBoundary],
    metric_label: str,
    outpath: Path,
    show_labels: bool = False,
    dpi: int = 150,
) -> None:
    """Render one dot-plot PNG with the alignment path and tracepoint grid."""
    fig, ax = plt.subplots(figsize=(10, 10))
    _render_dotplot_on_ax(
        ax, alignment, segments, boundaries,
        metric_label, show_labels,
    )
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def plot_combined(
    alignment: PafAlignment,
    segments: List[DotPlotSegment],
    panels: List[Tuple[str, List[TracepointBoundary]]],
    outpath: Path,
    show_labels: bool = False,
    dpi: int = 150,
) -> None:
    """Render a 1-row x 3-column combined dot-plot PNG.

    *panels* is a list of ``(metric_label, boundaries)`` tuples, one per
    column.  Order: FL-TP, EB-TP, DB-TP.
    """
    ncols = len(panels)
    fig, axes = plt.subplots(1, ncols, figsize=(11 * ncols, 14))
    if ncols == 1:
        axes = [axes]
    for i, (ax, (metric_label, boundaries)) in enumerate(zip(axes, panels)):
        _render_dotplot_on_ax(
            ax, alignment, segments, boundaries,
            metric_label, show_labels,
            font_scale=2.0,
        )
        # Only show y-axis label on the leftmost panel
        if i > 0:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)
    # Add colored CIGAR annotation at the top
    _render_cigar_annotation(
        axes[ncols // 2], fig, alignment.cigar_ops, font_scale=2.0
    )
    fig.tight_layout(w_pad=0.0)
    fig.subplots_adjust(top=0.78, wspace=0.05)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def _find_cigzip() -> str:
    """Locate the cigzip binary: check PATH, then ./target/release/cigzip."""
    path = shutil.which("cigzip")
    if path:
        return path
    local = Path("target/release/cigzip")
    if local.is_file():
        return str(local)
    print(
        "Error: cigzip binary not found in PATH or ./target/release/cigzip.\n"
        "Build it with `cargo build --release` or pass --cigzip /path/to/cigzip.",
        file=sys.stderr,
    )
    sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CIGAR dot-plot with tracepoint grids",
    )
    parser.add_argument("paf", help="PAF file with cg:Z: CIGAR tags")
    parser.add_argument(
        "--max-complexity",
        "-c",
        type=int,
        default=10,
        help="Maximum complexity for standard tracepoint segmentation (default: 10)",
    )
    parser.add_argument(
        "--fastga-spacing",
        "-f",
        type=int,
        default=100,
        help="Trace spacing for FastGA mode (default: 100)",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        default="dotplots",
        help="Output directory (default: dotplots)",
    )
    parser.add_argument(
        "--labels",
        action="store_true",
        help="Show (a_len, b_len) text at each tracepoint cell center",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Image DPI (default: 150)",
    )
    parser.add_argument(
        "--max-alignments",
        type=int,
        default=None,
        help="Limit number of alignments processed",
    )
    parser.add_argument(
        "--cigzip",
        default=None,
        help="Path to cigzip binary (default: search PATH, then ./target/release/cigzip)",
    )
    args = parser.parse_args()

    cigzip_bin = args.cigzip if args.cigzip else _find_cigzip()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # -- Run cigzip encode for each metric --
    # (code_name, tracepoint_type, cli_metric_or_None, paper_abbrev, param_symbol)
    metrics = [
        ("EditDistance",     "standard", "edit-distance",     "EB-TP", r"\delta"),
        ("DiagonalDistance", "standard", "diagonal-distance", "DB-TP", "b"),
        ("FastGA",          "fastga",   None,                "FL-TP", "l"),
    ]
    tracepoints_by_metric: dict[str, List[CigzipResult]] = {}
    for display_name, tp_type, cli_metric, _abbrev, _sym in metrics:
        label = cli_metric if cli_metric else "fastga"
        print(f"Running cigzip encode ({label})...", file=sys.stderr)
        complexity = args.fastga_spacing if tp_type == "fastga" else args.max_complexity
        tracepoints_by_metric[display_name] = run_cigzip_encode(
            args.paf,
            complexity,
            cigzip_bin,
            tracepoint_type=tp_type,
            metric=cli_metric,
            fastga_spacing=args.fastga_spacing if tp_type == "fastga" else None,
        )

    # -- Read input PAF for CIGAR ops and coordinates --
    alignments: List[PafAlignment] = []
    with open(args.paf) as fh:
        for line in fh:
            aln = parse_paf_line(line)
            if aln is not None:
                alignments.append(aln)

    if args.max_alignments is not None:
        alignments = alignments[: args.max_alignments]

    # -- Generate dot plots --
    # Order for the combined panel: FastGA, EditDistance, DiagonalDistance
    combined_order = ["FastGA", "EditDistance", "DiagonalDistance"]
    n_images = 0

    for idx, aln in enumerate(alignments):
        segments = cigar_to_dotplot_segments(
            aln.cigar_ops, aln.target_start, aln.query_start
        )

        # Collect per-metric boundaries for the combined plot.
        panels: List[Tuple[str, str, List[TracepointBoundary]]] = []

        for display_name, tp_type, _cli_metric, paper_abbrev, param_symbol in metrics:
            tp_list = tracepoints_by_metric[display_name]
            if idx >= len(tp_list):
                print(
                    f"Warning: no cigzip output for alignment {idx} "
                    f"({display_name}), skipping.",
                    file=sys.stderr,
                )
                continue
            cr = tp_list[idx]
            # FastGA may clip the alignment prefix; recompute boundaries
            # from the CIGAR directly to avoid coordinate-dependent clipping.
            if tp_type == "fastga" and (
                cr.query_start != aln.query_start
                or cr.target_start != aln.target_start
            ):
                boundaries = _fastga_boundaries_from_cigar(
                    aln.cigar_ops, aln.query_start, aln.target_start,
                    args.fastga_spacing,
                )
            else:
                boundaries = _pairs_to_boundaries(
                    cr.pairs, cr.query_start, cr.target_start
                )

            complexity = (
                args.fastga_spacing if tp_type == "fastga" else args.max_complexity
            )

            # -- Individual plot --
            metric_label = f"{paper_abbrev} (${param_symbol}$={complexity})"
            safe_q = re.sub(r"[^\w\-.]", "_", aln.query_name)[:40]
            safe_t = re.sub(r"[^\w\-.]", "_", aln.target_name)[:40]
            fname = f"{idx:04d}_{safe_q}_{safe_t}_{paper_abbrev}.pdf"
            outpath = outdir / fname
            plot_dotplot(
                aln,
                segments,
                boundaries,
                metric_label,
                outpath,
                show_labels=args.labels,
                dpi=args.dpi,
            )
            print(f"[{idx}] {outpath}")
            n_images += 1

            panels.append((display_name, metric_label, boundaries))

        # -- Combined 1x3 plot (FastGA, EditDistance, DiagonalDistance) --
        if panels:
            panel_map = {name: (label, b) for name, label, b in panels}
            ordered = [panel_map[k] for k in combined_order if k in panel_map]
            if ordered:
                safe_q = re.sub(r"[^\w\-.]", "_", aln.query_name)[:40]
                safe_t = re.sub(r"[^\w\-.]", "_", aln.target_name)[:40]
                combo_path = outdir / "cigar_dotplot_combined.pdf"
                plot_combined(
                    aln,
                    segments,
                    ordered,
                    combo_path,
                    show_labels=args.labels,
                    dpi=args.dpi,
                )
                print(f"[{idx}] {combo_path}")
                n_images += 1

    print(
        f"\nDone. {len(alignments)} alignment(s) processed, "
        f"{n_images} image(s) written to {outdir}/",
    )


if __name__ == "__main__":
    main()
