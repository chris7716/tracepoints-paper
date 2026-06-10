#!/usr/bin/env bash
set -euo pipefail

# Directories
REPO_DIR="${TRACEPOINTS_REPO_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"

# ── Rebuild LaTeX PDFs ───────────────────────────────────────────────────────
cd "${REPO_DIR}/paper"

echo "--- Building paper-main.pdf ---"
pdflatex -interaction=nonstopmode paper-main.tex
bibtex paper-main
pdflatex -interaction=nonstopmode paper-main.tex
pdflatex -interaction=nonstopmode paper-main.tex

echo "--- Building paper-supplementary.pdf ---"
pdflatex -interaction=nonstopmode paper-supplementary.tex
bibtex paper-supplementary
pdflatex -interaction=nonstopmode paper-supplementary.tex
pdflatex -interaction=nonstopmode paper-supplementary.tex

# ── Verification ─────────────────────────────────────────────────────────────
echo ""
echo "Generated PDFs:"
ls -la "${REPO_DIR}/paper/paper-main.pdf" "${REPO_DIR}/paper/paper-supplementary.pdf"

echo ""
echo "Done."
