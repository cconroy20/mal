#!/usr/bin/env python3
"""Sweep the gf-term weight lambda in the combined energy+gf fit (gf_fit.py) and
report the level / gf RMS at each, so the trade-off curve is visible and the best
lambda can be picked. Saves the fitted ING11 per lambda.

lambda controls how hard the fit chases gf vs energies:
  lambda -> 0   : energy-only (recovers the RCE-style fit; gf degrades)
  lambda large  : gf-dominated (gf improves, levels can drift)
The useful regime is where gf is well below the ab-initio floor AND the level
RMS stays small.

Usage:
  tools/gf_fit_sweep.py --run-dir work/mg1 --seed work/mg1/ING11.fit \
    --nist data/nist/MgI_levels.tsv --nist-lines data/nist/MgI_lines.tsv \
    --lambdas 0.1,0.3,1,3,10 --maxiter 8000 --out-prefix work/mg1/ING11.gffit
"""
import argparse
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ing11_params as IP
from gf_fit import Forward, make_objective


def fit_one(fwd, lam, ridge, maxiter):
    obj, scale = make_objective(fwd, lam, ridge)
    x0 = np.ones_like(fwd.seed)
    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"maxiter": maxiter, "maxfev": maxiter,
                            "xatol": 1e-4, "fatol": 1e-3, "adaptive": True})
    best = res.x * scale
    # evaluate metrics at the optimum
    fwd.run(best)
    eres = fwd.energy_resid()
    gres, gsig = fwd.gf_resid()
    erms = float(np.sqrt(np.mean(eres ** 2))) if len(eres) else float("nan")
    grms = float(np.sqrt(np.mean(gres ** 2))) if len(gres) else float("nan")
    wgrms = (float(np.sqrt(np.sum((gres / gsig) ** 2) / np.sum(1 / gsig ** 2)))
             if len(gres) else float("nan"))
    return best, erms, grms, wgrms, res.fun


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--nist", required=True)
    ap.add_argument("--nist-lines", required=True)
    ap.add_argument("--lambdas", default="0.1,0.3,1,3,10")
    ap.add_argument("--ridge", type=float, default=0.0)
    ap.add_argument("--maxiter", type=int, default=8000)
    ap.add_argument("--out-prefix", default=None,
                    help="write each fit to <prefix>.lam<val>")
    a = ap.parse_args()

    lams = [float(x) for x in a.lambdas.split(",")]
    fwd = Forward(a.run_dir, a.seed, a.nist, a.nist_lines)
    print(f"{len(fwd.params)} params; seed {a.seed}\n")
    print(f"{'lambda':>8} {'levelRMS':>10} {'gfRMS':>8} {'acc-gfRMS':>10} "
          f"{'chi2':>10} {'nfev':>7}")
    print("-" * 60)
    rows = []
    for lam in lams:
        fwd.neval = 0
        best, erms, grms, wgrms, chi2 = fit_one(fwd, lam, a.ridge, a.maxiter)
        print(f"{lam:8.3g} {erms:10.1f} {grms:8.3f} {wgrms:10.3f} "
              f"{chi2:10.2f} {fwd.neval:7d}")
        rows.append((lam, erms, grms, wgrms))
        if a.out_prefix:
            out = f"{a.out_prefix}.lam{lam:g}"
            IP.write(fwd.raw, fwd.params, best, out)
    # also report the seed for reference
    print("-" * 60)
    return rows


if __name__ == "__main__":
    main()
