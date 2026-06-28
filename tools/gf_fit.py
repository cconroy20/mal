#!/usr/bin/env python3
"""Tier-2 prototype: fit the Cowan radial parameters to BOTH observed energy
levels AND well-measured NIST gf values, using the existing RCG as a black-box
forward model wrapped in a Python optimizer.

Motivation: the standard RCE fit (energies only) was found to DEGRADE gf vs the
ab-initio RCG values (notes/mg1_gf_analysis.md), because fitting energies
reshapes eigenvectors in ways that hurt the dipole matrix elements. Here we put
gf into the objective so the fit can't trade away gf accuracy for a marginally
better energy.

Objective (minimized):
    chi2 = sum_i wE_i (E_calc_i - E_obs_i)^2            [levels, cm^-1]
         + lambda * sum_j wgf_j (loggf_calc_j - loggf_NIST_j)^2   [gf, dex]
  wE_i  : 1/sigma_E^2, sigma_E a flat energy tolerance (cm^-1)
  wgf_j : 1/sigma_gf^2 from NIST accuracy class; only log gf >= GF_MIN included
  lambda: relative weight of the gf term (tunable)

Forward model: write a trial parameter vector into ING11 via ing11_params, run
RCG (~8 ms for Mg I), parse levels + gf from OUTG11, match to NIST by
eigenvector-composition identity (make_report).

Usage:
    tools/gf_fit.py --run-dir work/mg1 --seed work/mg1/ING11.abinitio \
        --nist data/nist/MgI_levels.tsv --nist-lines data/nist/MgI_lines.tsv \
        --lambda 1.0 [--ridge 0.0] [--maxiter 4000]
"""
import argparse
import os
import subprocess
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ing11_params as IP
import make_report as R
from parse_cowan import parse_outg11, parse_compositions

# NIST accuracy class -> approx 1-sigma uncertainty in log gf (dex), ASD legend.
ACC_DEX = {"AAA": 0.013, "AA": 0.013, "A+": 0.013, "A": 0.022, "B+": 0.043,
           "B": 0.087, "C+": 0.13, "C": 0.22, "D+": 0.30, "D": 0.43, "E": 0.70}
GF_MIN = -1.0          # only fit gf for lines at least this strong
SIGMA_E = 50.0         # flat energy uncertainty (cm^-1) for the level term


def acc_sigma(a):
    return ACC_DEX.get(a.strip(), 0.5)


class Forward:
    """Black-box forward model: params -> RCG -> (levels, gf), with NIST targets
    precomputed so each evaluation just runs RCG and matches."""

    def __init__(self, run_dir, seed_ing11, nist_path, nist_lines_path):
        self.run_dir = os.path.abspath(run_dir)
        self.ing11 = os.path.join(self.run_dir, "ING11")
        self.outg11 = os.path.join(self.run_dir, "OUTG11")
        self.rcg = os.path.join(os.path.dirname(self.run_dir), "..",
                                "build", "bin", "rcg")
        self.rcg = os.path.abspath(self.rcg)
        self.raw, self.params = IP.parse(seed_ing11)
        self.seed = np.array([p["value"] for p in self.params])
        self.nist = R.load_nist(nist_path)
        self.nist_lines = R.load_nist_lines(nist_lines_path)
        # NIST level lookup by robust identity (cfgkey, termkey, Jkey)
        self.nist_lev = {}
        for nlv in self.nist:
            k = R._level_idkey(nlv["config"], nlv["term"], nlv["J"])
            self.nist_lev.setdefault(k, nlv["E_obs"])
        # ensure cowan.cfg points RCG at run_dir
        with open(os.path.join(self.run_dir, "cowan.cfg"), "w") as f:
            f.write(self.run_dir + "/\n")
        self.neval = 0

    def run(self, values):
        IP.write(self.raw, self.params, values, self.ing11)
        subprocess.run([self.rcg], cwd=self.run_dir, stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.neval += 1
        return None

    def energy_resid(self):
        """Level residuals (E_calc - E_obs), identifying computed levels by
        eigenvector composition (the OUTG11 header config labels are unreliable),
        and removing a single global offset (= the ground-EAV reference, which
        RCG-from-fitted-params carries but observed energies don't)."""
        comp = parse_compositions(self.outg11)
        res = []
        for (par, J), levs in comp.items():
            for L in levs:
                k = R._level_idkey(L["config"], L["term"], J)
                eobs = self.nist_lev.get(k)
                if eobs is not None:
                    res.append(L["E_calc"] - eobs)
        res = np.array(res)
        if len(res):
            res = res - np.median(res)      # absorb the global energy offset
        return res

    def gf_resid(self):
        """(residuals, sigmas) for matched strong lines from the CURRENT OUTG11."""
        pairs = R.match_gf_by_identity(self.outg11, self.nist_lines)
        res, sig = [], []
        # match_gf_by_identity returns (loggf_comp, loggf_nist, lam); we need acc
        # too, so re-run the richer matcher from gf_table.
        from gf_table import matched_table
        rows, _ = matched_table(self.outg11, self.nist_lines)
        for r in rows:
            if r["nist"] < GF_MIN:
                continue
            res.append(r["d"])
            sig.append(acc_sigma(r["acc"]))
        return np.array(res), np.array(sig)


def make_objective(fwd, lam, ridge):
    seed = fwd.seed
    scale = np.where(np.abs(seed) > 1e-6, np.abs(seed), 1.0)

    def obj(x):
        values = x * scale          # x is in seed-relative units
        fwd.run(values)
        eres = fwd.energy_resid()
        chi2 = np.sum((eres / SIGMA_E) ** 2) if len(eres) else 0.0
        gres, gsig = fwd.gf_resid()
        if len(gres):
            chi2 += lam * np.sum((gres / gsig) ** 2)
        if ridge:
            chi2 += ridge * np.sum((x - 1.0) ** 2)   # pull toward seed
        return chi2

    return obj, scale


def report(fwd, values, label):
    fwd.run(values)
    eres = fwd.energy_resid()
    gres, gsig = fwd.gf_resid()
    erms = np.sqrt(np.mean(eres ** 2)) if len(eres) else float("nan")
    grms = np.sqrt(np.mean(gres ** 2)) if len(gres) else float("nan")
    wgrms = (np.sqrt(np.sum((gres / gsig) ** 2) / np.sum(1 / gsig ** 2))
             if len(gres) else float("nan"))
    print(f"{label:>16}: level RMS = {erms:7.1f} cm^-1   "
          f"strong-gf RMS = {grms:.3f}   acc-wtd gf RMS = {wgrms:.3f}   "
          f"(N_E={len(eres)}, N_gf={len(gres)})")
    return erms, grms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--seed", required=True, help="ING11 to seed from.")
    ap.add_argument("--nist", required=True)
    ap.add_argument("--nist-lines", required=True)
    ap.add_argument("--lambda", dest="lam", type=float, default=1.0)
    ap.add_argument("--ridge", type=float, default=0.0)
    ap.add_argument("--maxiter", type=int, default=4000)
    ap.add_argument("--out", default=None, help="write fitted ING11 here.")
    a = ap.parse_args()

    fwd = Forward(a.run_dir, a.seed, a.nist, a.nist_lines)
    print(f"{len(fwd.params)} adjustable params; seeded from {a.seed}")
    report(fwd, fwd.seed, "seed")

    obj, scale = make_objective(fwd, a.lam, a.ridge)
    x0 = np.ones_like(fwd.seed)
    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"maxiter": a.maxiter, "maxfev": a.maxiter,
                            "xatol": 1e-4, "fatol": 1e-3, "adaptive": True})
    best = res.x * scale
    print(f"\noptimizer: {res.message}  (nfev={fwd.neval}, chi2={res.fun:.3f})")
    report(fwd, best, "gf-fit")

    if a.out:
        IP.write(fwd.raw, fwd.params, best, a.out)
        print(f"wrote fitted parameters -> {a.out}")


if __name__ == "__main__":
    main()
