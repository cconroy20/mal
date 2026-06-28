#!/bin/bash
# Build the LaTeX docs to PDF. MacTeX lives at /Library/TeX/texbin (not always on PATH).
set -e
export PATH="/Library/TeX/texbin:$PATH"
cd "$(dirname "$0")"
for tex in *.tex; do
  echo "Building $tex ..."
  latexmk -pdf -interaction=nonstopmode -halt-on-error "$tex" >/dev/null
  echo "  -> ${tex%.tex}.pdf"
done
latexmk -c >/dev/null 2>&1 || true   # clean aux files, keep the PDFs
echo "Done."
