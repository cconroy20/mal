#!/bin/bash
# Regenerate the Mg I diagnostic report from the current run outputs.
# The gf page shows the combined energy+gf fit (gf_fit.py, lambda=3) when
# work/mg1/OUTG11.gffit exists; otherwise it falls back to the energy-only RCE
# fitted gf (OUTG11.fitted). Produce OUTG11.gffit with run_mg1.sh GFFIT=1 (or
# tools/gf_fit.py directly).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if [ -f "$ROOT/work/mg1/OUTG11.gffit" ]; then
  GF_OUT="$ROOT/work/mg1/OUTG11.gffit"; GF_LABEL="fitted (energy+gf)"
else
  GF_OUT="$ROOT/work/mg1/OUTG11.fitted"; GF_LABEL="fitted (RCE)"
fi
python3 "$ROOT/tools/make_report.py" "Mg I" \
  --outg11          "$ROOT/work/mg1/OUTG11.abinitio" \
  --nist            "$ROOT/data/nist/MgI_levels.tsv" \
  --levels1         "$ROOT/work/mg1/LEVELS1" \
  --ie              "$ROOT/data/nist/ionization.tsv" \
  --gf-fitted       "$GF_OUT" \
  --gf-fitted-label "$GF_LABEL" \
  --nist-lines      "$ROOT/data/nist/MgI_lines.tsv" \
  --out             "$ROOT/docs/reports/MgI_report.pdf"
