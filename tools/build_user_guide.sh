#!/usr/bin/env bash
# Build docs/USER_GUIDE.pdf from the single source USER_GUIDE.md (pandoc + LuaLaTeX).
# Run from anywhere:  bash tools/build_user_guide.sh
# Edit USER_GUIDE.md, never the PDF.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

COMMIT="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
DATE="$(date +%Y-%m-%d)"

pandoc USER_GUIDE.md \
  --from=gfm+tex_math_dollars \
  --pdf-engine=lualatex \
  --toc --toc-depth=2 \
  --number-sections \
  -V documentclass=article -V papersize=a4 -V fontsize=10pt \
  -V geometry:margin=2.2cm \
  -V mainfont="STIX Two Text" -V mathfont="STIX Two Math" \
  -V monofont="Menlo" -V monofontoptions="Scale=0.85" \
  -V colorlinks=true -V linkcolor=blue -V urlcolor=blue \
  -M title="User Guide: Compass-Needle Lattice Simulator" \
  -M author="J. P. Sinnecker" \
  -M date="Generated ${DATE} from USER_GUIDE.md (commit ${COMMIT})" \
  -H tools/user_guide_header.tex \
  -o docs/USER_GUIDE.pdf

echo "wrote docs/USER_GUIDE.pdf (commit ${COMMIT})"
