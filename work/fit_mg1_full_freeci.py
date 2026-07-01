#!/usr/bin/env python3
"""Full-basis (122-config) Mg I fit, Bob-style background + SELECTIVELY FREED CI.

Path 2 of the settled diagnosis (notes/mg1_gf_analysis.md): the persistent ~0.5
dex on the singlet 3s.nd 1D -> 3s3p 1P and 3s2 -> 3s4p 1P lines comes from the
1D/1P EIGENVECTOR MIXING (3s.nd / 3s.np Rydberg with 3p^2 / 3d^2), which is set
by the FROZEN CI (R^k) structure. Bob nails these with curated per-integral CI;
our uniform-scaled CI can't reproduce the mixing, and freezing CI means the fit
can't fix it. Our 9-config win came from FREEING the key low CI -- which directly
reshapes this mixing. Here we do the same on the full basis: keep Bob's
scaled-HF background + the EAV+P ruleset free set, and ADDITIONALLY free the few
DOMINANT low-l CI integrals coupling the singlet-system configs. The HF ridge
prior (PRIOR_SIGMA['CI']=0.15, tightest) keeps every other CI pinned near its
scaled-HF value, so only the named pairs move -- avoiding the divergence of
freeing all correlated CI at once.

Seeds from ING11.scaled (Rule-2 scaled-HF background; fixed the energies).
Run detached:  cd ~/kurucz/mal && nohup python3 work/fit_mg1_full_freeci.py &
Watch:         tail -f work/mg1_full/fit_freeci_progress.log
"""
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import numpy as np
import ing11_params as IP
import make_report as R
from gf_fit import Forward, fit_lm
from gf_table import matched_table

RUN = os.path.join(ROOT, "work", "mg1_full")
NIST = os.path.join(ROOT, "data", "nist", "MgI_levels.tsv")
NL = os.path.join(ROOT, "data", "nist", "MgI_lines.tsv")
# Ridge width on the FREED CI integrals. The default CI prior (0.15) pins them
# near scaled-HF -- too tight to let the fit make the ~60% move our 9-config fit
# found. Loosen via FREECI_SIGMA env var (e.g. 2.0 ~ effectively free) to test
# whether the gf-degrading 1D/1P mixing is correctable when the prior gets out of
# the way (or whether the limiter is the weak data signal, not the prior).
FREECI_SIGMA = float(os.environ.get("FREECI_SIGMA", "0.15"))
TAG = os.environ.get("FREECI_TAG", "freeci")
LOG = os.path.join(RUN, f"fit_{TAG}_progress.log")
nistl = R.load_nist_lines(NL)

# Bob-style scaled-HF background (Rule 2) -- the run that fixed the energies but
# left gf at 0.212. We add free CI on TOP of it.
SEED = os.path.join(RUN, "ING11.scaled")
ABINIT = os.path.join(RUN, "ING11.scaled")   # ridge prior centre = scaled HF

# The DOMINANT low-l CI integrals that set the singlet 1D/1P/1S mixing. Each is a
# frozenset of two config-name strings (matched against the CI param key, which
# carries the readable pair, e.g. "3s2 - 3p2"). Chosen from the magnitude ranking
# of CI integrals coupling the singlet-system configs (3s2 / 3p2 / 3d2 / 3s.nd):
#   3s2-3p2 (22.05, the smoking gun, our 9-config fit cut it 59%)
#   3s2-3d2 (7.73), 3p2-3d2 (13.11)            -- 3s2 / 3d2 channel
#   3p2-{3s3d,3s4d,3s5d,3s6d}                  -- 3p2 perturbs the 1D Rydberg series
#   3d2-{3s3d,3s4d,3s5d,3s6d}                  -- 3d2 perturbs the 1D Rydberg series
ND = ["3s3d", "3s4d", "3s5d", "3s6d"]
FREE_CI_PAIRS = {
    frozenset(("3s2", "3p2")),
    frozenset(("3s2", "3d2")),
    frozenset(("3p2", "3d2")),
} | {frozenset(("3p2", nd)) for nd in ND} \
  | {frozenset(("3d2", nd)) for nd in ND}


def log(msg):
    with open(LOG, "a") as f:
        f.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")


def gfrms(outg11):
    rows, _ = matched_table(outg11, nistl)
    d = [r["d"] for r in rows
         if r["nist"] >= -1 and r["acc"].strip()[:1] in ("A", "B")]
    return (np.sqrt(np.mean(np.square(d))), len(d)) if d else (float("nan"), 0)


def main():
    open(LOG, "w").close()
    obs = {R._cfgkey(n["config"]) for n in R.load_nist(NIST)}
    t0 = time.time()
    fwd = Forward(RUN, SEED, NIST, NL, free_kinds={"EAV", "P"},
                  abinitio_ing11=ABINIT, obs_configs=obs,
                  free_ci_pairs=FREE_CI_PAIRS, free_ci_sigma=FREECI_SIGMA)
    nci = sum(1 for p in fwd.params if p["kind"] == "CI")
    log(f"free CI pairs requested: {sorted(tuple(sorted(p)) for p in FREE_CI_PAIRS)}")
    log(f"  -> {nci} CI params freed; FREECI_SIGMA={FREECI_SIGMA}")
    g0, n0 = gfrms(fwd.outg11)
    e0 = float(np.sqrt(np.mean(fwd.energy_resid() ** 2)))
    log(f"SEED(scaled+freeCI) {len(fwd.params)} params  levelRMS={e0:.1f}  "
        f"gfRMS(A/B)={g0:.3f} N={n0}")

    best, res, scale, cov = fit_lm(fwd, lam=3, ridge=1.0, maxiter=3000,
                                   presolve=False, progress=log)
    e1 = float(np.sqrt(np.mean(fwd.energy_resid() ** 2)))
    g1, n1 = gfrms(fwd.outg11)
    log(f"FIT levelRMS={e1:.1f}  gfRMS(A/B)={g1:.3f} N={n1}  "
        f"nfev={fwd.neval}  cov={'ok' if cov is not None else 'singular'}  "
        f"({time.time() - t0:.0f}s)")

    IP.write(fwd.raw, fwd.params, best, os.path.join(RUN, f"ING11.{TAG}"))
    fwd.run(best)
    os.replace(fwd.outg11, os.path.join(RUN, f"OUTG11.{TAG}"))
    log(f"wrote ING11.{TAG} + OUTG11.{TAG}  -- DONE")


if __name__ == "__main__":
    main()
