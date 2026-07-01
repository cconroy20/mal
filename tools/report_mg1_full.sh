#!/bin/bash
# Regenerate the Mg I diagnostic report for the FULL-BASIS (122-config) model with
# the IDENTITY-MATCHED, ENERGY-ONLY (gf-off) fit -- the current best Mg I model.
# v3 (2026-06-30): the energy cap was raised from ~58000 to 60000 cm^-1, admitting
# the well-placed n<=8 high-Rydberg levels that the basis CAN reach (the cap had
# been needlessly discarding ~38 good levels); the fit now constrains 87 levels
# (vs ~50), median |E-O| ~2 cm^-1, report-convention level-RMS 31.6 cm^-1. Page 1
# now carries a head-to-head vs Bob Kurucz's RCE fit (--compare-cap). The cap stops
# at 60000, NOT the IE: the topmost 3s.9d/3s.10d Rydberg members are beyond the
# 122-config basis's reach (E_calc off by ~2000-12000) and would swamp the RMS.
# This supersedes tools/report_mg1.sh (the older 9-config work/mg1 run).
#
# The fitted spectrum is gf_fit's OUTG11.energyonly (Python optimizer, no RCE
# LEVELS1), so the level/residual pages are fed a fitted-levels JSON built from it
# by tools/fitted_levels_json.py (identity-matched E_fit per NIST level).
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RUN="$ROOT/work/mg1_full"
NIST="$ROOT/data/nist/MgI_levels.tsv"
# Cap at 60000 cm^-1 (not the IE 61671): the top Rydberg members 3s.9d/10d/11d are
# observed below the IE but BEYOND the 122-config basis's reliable reach (E_calc
# off by ~2000-12000), so including them would swamp the report with basis-limit
# artifacts. 60000 admits every level the basis CAN place (n<=8 d, n<=8 p, etc.).
IE_CM=60000

# (re)build the fitted-levels JSON from the energy-only fit's OUTG11
python3 "$ROOT/tools/fitted_levels_json.py" \
  --outg11 "$RUN/OUTG11.energyonly" --nist "$NIST" --max-energy "$IE_CM" \
  --out "$RUN/fitted_levels.json"

# Compute the headline level-fit RMS FROM the fitted JSON (offset-removed), so it
# tracks the current model instead of a hardcoded value that goes stale.
FIT_RMS=$(python3 -c "import json,numpy as np;d=json.load(open('$RUN/fitted_levels.json'));r=np.array([x['E_fit']-x['E_obs'] for x in d if x['E_fit'] is not None]);r=r-np.median(r);print(f'{np.sqrt(np.mean(r**2)):.1f}')")

# Bob Kurucz's own fitted gf + level residuals, added as reference series.
KURUCZ_ARGS=()
[ -f "$ROOT/kurucz_ref/1200/gf1200.pos" ] && KURUCZ_ARGS=(
  --kurucz-lines "$ROOT/kurucz_ref/1200/gf1200.pos" --kurucz-elem 12.00)
KLEV=()
for f in "$ROOT"/kurucz_ref/1200/c1200ez.log "$ROOT"/kurucz_ref/1200/c1200oz.log; do
  [ -f "$f" ] && KLEV+=("$f")
done
[ ${#KLEV[@]} -gt 0 ] && KURUCZ_ARGS+=(--kurucz-levels "${KLEV[@]}")

python3 "$ROOT/tools/make_report.py" "Mg I" \
  --outg11             "$RUN/OUTG11.abinitio" \
  --nist               "$NIST" \
  --fitted-levels-json "$RUN/fitted_levels.json" \
  --fit-rms            "$FIT_RMS" \
  --compare-cap        "$IE_CM" \
  --ie                 "$ROOT/data/nist/ionization.tsv" \
  --gf-fitted          "$RUN/OUTG11.energyonly" \
  --gf-fitted-label    "fitted (energy-only, full basis)" \
  --fitted-label       "fitted (energy-only, full basis)" \
  --nist-lines         "$ROOT/data/nist/MgI_lines.tsv" \
  "${KURUCZ_ARGS[@]}" \
  --out                "$ROOT/docs/reports/MgI_report.pdf"
