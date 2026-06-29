#!/usr/bin/env python3
"""Apply Bob Kurucz's HF SCREENING SCALE FACTORS to an ING11's ab-initio radial
integrals (ruleset Rule 2: freeze the structural background at SCALED HF, not raw
HF). Bob's per-family scales, read from his c1200*.log FIXEDHF columns:
    EAV  1.0   (centroids -- not screened)
    ZETA 1.0   (spin-orbit -- relativistic, not screened)
    F^k  0.8   (direct Slater)
    G^k  0.6-0.8 (exchange Slater; ~0.7 typical)
    R^k  0.7-0.8 (configuration interaction; ~0.75 typical)
Scaling CI+Slater down corrects HF's systematic overestimate (electron
correlation screens the electron-electron integrals). A SINGLE global scale fails
(F and G want different corrections); per-family is the minimal faithful version.

Caveat: ING11 lumps ZETA with F/G as kind 'P' and doesn't name them, so without
the OUTGINE name map we scale all 'P' (Slater AND zeta) by --slater. Since zeta
is a minority and is freed for observed configs anyway, this is an acceptable
first approximation; a name-aware version can exempt zeta later.

Usage:
  tools/scale_hf.py --in work/mg1_full/ING11.abinitio \
      --out work/mg1_full/ING11.scaled --ci 0.75 --slater 0.8
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ing11_params as IP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ci", type=float, default=0.75, help="R^k (CI) scale")
    ap.add_argument("--slater", type=float, default=0.8,
                    help="F^k/G^k (+zeta) scale for kind 'P'")
    a = ap.parse_args()
    raw, params = IP.parse(a.inp)
    vals, nci, np_ = [], 0, 0
    for p in params:
        if p["kind"] == "CI":
            vals.append(p["value"] * a.ci); nci += 1
        elif p["kind"] == "P":
            vals.append(p["value"] * a.slater); np_ += 1
        else:                                   # EAV: unscaled
            vals.append(p["value"])
    IP.write(raw, params, vals, a.out)
    print(f"wrote {a.out}  (scaled {nci} CI x{a.ci}, {np_} Slater/zeta "
          f"x{a.slater}; EAV unscaled)")


if __name__ == "__main__":
    main()
