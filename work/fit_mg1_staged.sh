#!/bin/bash
# STAGED (two-phase) Mg I energy fit -- the "fit the strong params first, then
# release the weak ones" protocol.
#
# Phase 1 (STAGE=strong): free the EAV only of configs with a DIRECTLY-OBSERVED
#   level <= cap (centroid pinned by data) + 3p2. Freezes the ~35 targetless
#   ("weak") EAVs -- 3p3d, 3p4s, the 3p.nd/3p.ns doubly-excited manifolds -- at
#   scaled-HF. Converges the well-determined structure in a well-conditioned basin.
#   Writes ING11.energyonly_phase1.
# Phase 2: SEED from ING11.energyonly_phase1, free the FULL set, polish + covariance.
#   Writes ING11.energyonly / OUTG11.energyonly (the canonical fitted model).
#
# The hypothesis under test: does converging the strong set first stop the weak
# EAVs (3p3d etc.) from wandering, vs freeing everything at once (baseline: RMS
# 87.5, median 3.3)? Compare the phase-2 result AND the 3p3d residuals to that.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/.."
RUN="$ROOT/work/mg1_full"

COMMON="FIT_GF_SCALES=1 FREE_G=3s.3p,3s.4p,3s.5p,3s.6p FREE_CI=3s2-3p2 \
FREE_CI_SIGMA=1.0 ZETA_SIGMA=2.0 ZETA_CONFIGS=3s.3p,3p2 MAXE=60000 NJAC=6 TOL=1e-6"

echo "=== PHASE 1: strong params only (weak EAVs frozen at scaled-HF) ==="
env $COMMON STAGE=strong SEED=ING11.scaled \
  python3 "$HERE/fit_mg1_full_energyonly.py"
grep -E "free params|FIT " "$RUN/fit_energyonly_phase1_progress.log" | tail -2

echo "=== PHASE 2: release the full set, seeded from phase 1 ==="
# SEED from phase 1 (good basin); ridge CENTRE stays scaled-HF (physical prior).
env $COMMON SEED=ING11.energyonly_phase1 ABINIT=ING11.scaled \
  python3 "$HERE/fit_mg1_full_energyonly.py"
grep -E "free params|FIT " "$RUN/fit_energyonly_progress.log" | tail -2

echo "DONE. Compare phase-2 FIT line above to the single-phase baseline (RMS 87.5 / median 3.3)."
