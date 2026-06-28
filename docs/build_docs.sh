#!/bin/bash
# Build the LaTeX docs to PDF. MacTeX lives at /Library/TeX/texbin (not always on PATH).
set -e
export PATH="/Library/TeX/texbin:$PATH"
cd "$(dirname "$0")"

# regenerate figures first (idempotent)
if command -v python3 >/dev/null; then
  echo "Generating figures ..."
  [ -f make_figs.py ] && python3 make_figs.py            # illustrative (kernel, schematic)
  if [ -f ../work/sn7plus/OUTG11 ]; then
    python3 make_worked_figs.py                          # worked example (real run)
  else
    echo "  (skipping worked figs: run build/bin/{rcn,rcn2,rcg} in work/sn7plus first)"
  fi
fi

for tex in *.tex; do
  echo "Building $tex ..."
  latexmk -pdf -interaction=nonstopmode -halt-on-error "$tex" >/dev/null
  echo "  -> ${tex%.tex}.pdf"
done
latexmk -c >/dev/null 2>&1 || true   # clean aux files, keep the PDFs
echo "Done."
