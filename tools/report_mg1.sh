#!/bin/bash
# Regenerate the Mg I diagnostic report from the current run outputs.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$ROOT/tools/make_report.py" "Mg I" \
  --outg11 "$ROOT/work/mg1/OUTG11" \
  --nist   "$ROOT/data/nist/MgI_levels.tsv" \
  --out    "$ROOT/docs/reports/MgI_report.pdf"
