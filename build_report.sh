#!/bin/bash
# build_report.sh — Compile docs/project_report.tex to PDF
#
# Uses scratch dir for pdflatex outputs (avoids Dropbox sandbox restrictions),
# then copies the PDF back. Run this script from a regular Terminal.
#
# project_report.tex lives in docs/ and references figures as ../figs/...,
# so the scratch build mirrors that same docs/ + figs/ sibling layout
# rather than flattening everything into one directory.
#
# Usage:
#   bash build_report.sh
#

set -e
WS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$HOME/.gemini/antigravity-ide/scratch/texbuild"
if command -v pdflatex >/dev/null 2>&1; then
    PDFLATEX="$(command -v pdflatex)"
elif [ -f "/Library/TeX/texbin/pdflatex" ]; then
    PDFLATEX="/Library/TeX/texbin/pdflatex"
elif [ -f "/usr/local/bin/pdflatex" ]; then
    PDFLATEX="/usr/local/bin/pdflatex"
elif [ -f "/usr/bin/pdflatex" ]; then
    PDFLATEX="/usr/bin/pdflatex"
else
    echo "Error: pdflatex command was not found. Please install LaTeX." >&2
    exit 1
fi
TEX="project_report"

mkdir -p "$BUILD/docs"
cp "$WS/docs/${TEX}.tex" "$BUILD/docs/"
cp -r "$WS/figs" "$BUILD/" 2>/dev/null || true

cd "$BUILD/docs"
echo "=== Pass 1 ==="
"$PDFLATEX" -interaction=nonstopmode "${TEX}.tex"
echo "=== Pass 2 ==="
"$PDFLATEX" -interaction=nonstopmode "${TEX}.tex"

cp "${TEX}.pdf" "$WS/docs/"
echo ""
echo "✓ PDF compiled and copied to workspace:"
echo "  $WS/docs/${TEX}.pdf"
#open "$WS/docs/${TEX}.pdf"
