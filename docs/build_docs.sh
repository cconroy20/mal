#!/bin/bash
# Build the LaTeX docs to PDF. MacTeX lives at /Library/TeX/texbin (not always on PATH).
set -e
export PATH="/Library/TeX/texbin:$PATH"
cd "$(dirname "$0")"

# regenerate figures first (idempotent)
if command -v python3 >/dev/null && [ -f make_figs.py ]; then
  echo "Generating figures ..."
  python3 make_figs.py
fi

for tex in *.tex; do
  echo "Building $tex ..."
  latexmk -pdf -interaction=nonstopmode -halt-on-error "$tex" >/dev/null
  echo "  -> ${tex%.tex}.pdf"
done
latexmk -c >/dev/null 2>&1 || true   # clean aux files, keep the PDFs
echo "Done."
