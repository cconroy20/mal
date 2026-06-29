#!/usr/bin/env python3
"""Full-basis (122-config) Mg I fit, Bob-style: free EAV + Slater of OBSERVED
configs (170 params), freeze all CI + unobserved-config params at ab-initio,
regularize toward HF with the ridge prior. Writes per-eval progress to
work/mg1_full/fit_progress.log so the (slow ~1s/eval) run is watchable.

Run detached, e.g.:  cd ~/kurucz/mal && nohup python3 work/fit_mg1_full.py &
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
LOG = os.path.join(RUN, "fit_progress.log")
nistl = R.load_nist_lines(NL)


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
    abinit = os.path.join(RUN, "ING11.abinitio")
    t0 = time.time()
    fwd = Forward(RUN, abinit, NIST, NL, free_kinds={"EAV", "P"},
                  abinitio_ing11=abinit, obs_configs=obs)
    g0, n0 = gfrms(fwd.outg11)
    e0 = float(np.sqrt(np.mean(fwd.energy_resid() ** 2)))
    log(f"SEED  {len(fwd.params)} free params  levelRMS={e0:.1f}  "
        f"gfRMS(A/B)={g0:.3f} N={n0}")

    best, res, scale, cov = fit_lm(fwd, lam=3, ridge=1.0, maxiter=2000,
                                   presolve=False, progress=log)
    e1 = float(np.sqrt(np.mean(fwd.energy_resid() ** 2)))
    g1, n1 = gfrms(fwd.outg11)
    log(f"FIT   levelRMS={e1:.1f}  gfRMS(A/B)={g1:.3f} N={n1}  "
        f"nfev={fwd.neval}  cov={'ok' if cov is not None else 'singular'}  "
        f"({time.time() - t0:.0f}s)")

    IP.write(fwd.raw, fwd.params, best, os.path.join(RUN, "ING11.gffit"))
    fwd.run(best)
    os.replace(fwd.outg11, os.path.join(RUN, "OUTG11.gffit"))
    log("wrote ING11.gffit + OUTG11.gffit  -- DONE")


if __name__ == "__main__":
    main()
