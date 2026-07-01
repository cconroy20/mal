#!/usr/bin/env python3
"""Apply Bob Kurucz's HF SCREENING SCALE FACTORS to an ING11's ab-initio radial
integrals (ruleset Rule 2: freeze/centre the structural background at SCALED HF,
not raw HF). PER-FAMILY scales, from his c1200*.log FIXEDHF columns:
    EAV  1.0   (centroids -- not screened)
    ZETA 1.0   (spin-orbit -- relativistic, NOT screened down)
    F^k  0.8   (direct Slater)
    G^k  0.66  (exchange Slater; Bob's 3s.np G^1 ~0.66-0.785, ~0.7 typical)
    R^k  0.8   (configuration interaction)
Scaling the Slater/CI integrals down corrects HF's systematic overestimate
(correlation screens the e-e integrals); spin-orbit must NOT be scaled down.

NAME-AWARE: ING11 lumps ZETA with F/G as kind 'P' and doesn't name the slots, so
we read each slot's PHYSICAL name from the OUTGINE block headers (via param_labels)
and apply the per-family scale by name. This fixes the prior-CENTRE bias: the old
version scaled zeta down with F/G (--slater 0.8), putting the zeta prior centre at
~0.55x of truth, so the ridge then FROZE zeta there (the energy gain from fixing
fine structure was < the ridge penalty to move it 2.7 sigma). Un-scaling zeta puts
its centre at HF (~right) so the fit no longer fights the prior. Slots
param_labels can't name (open-shell overflow / gaps) fall back to --slater.

Usage:
  tools/scale_hf.py --in work/mg1_full/ING11.abinitio \
      --outgine work/mg1_full/OUTGINE.abinitio --out work/mg1_full/ING11.scaled \
      [--ci 0.8 --fdir 0.8 --gexch 0.66 --zeta 1.0]
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ing11_params as IP
import param_labels as PL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--outgine", default=None,
                    help="OUTGINE for physical slot names (enables per-family "
                         "scaling). Without it, all 'P' get --slater (old behavior).")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ci", type=float, default=0.8, help="R^k (CI) scale")
    ap.add_argument("--fdir", type=float, default=0.8, help="F^k direct Slater scale")
    ap.add_argument("--gexch", type=float, default=0.66, help="G^k exchange scale")
    ap.add_argument("--zeta", type=float, default=1.0, help="spin-orbit scale")
    ap.add_argument("--slater", type=float, default=0.8,
                    help="fallback scale for un-nameable 'P' slots.")
    a = ap.parse_args()
    raw, params = IP.parse(a.inp)

    # physical-name map: param key -> family scale, from OUTGINE (if given)
    key_scale = {}
    if a.outgine:
        for p in PL.physical_params(a.inp, a.outgine, strict=False):
            kind = p.get("kind")
            if kind == "ZETA":
                key_scale[p["key"]] = a.zeta
            elif kind == "F":
                key_scale[p["key"]] = a.fdir
            elif kind == "G":
                key_scale[p["key"]] = a.gexch

    vals = []
    tally = {"ZETA": 0, "F": 0, "G": 0, "CI": 0, "P-fallback": 0, "EAV": 0}
    for p in params:
        if p["kind"] == "CI":
            vals.append(p["value"] * a.ci); tally["CI"] += 1
        elif p["kind"] == "P":
            if p["key"] in key_scale:
                s = key_scale[p["key"]]
                vals.append(p["value"] * s)
                for fam, sc in (("ZETA", a.zeta), ("F", a.fdir), ("G", a.gexch)):
                    if s == sc:
                        tally[fam] += 1
                        break
            else:
                vals.append(p["value"] * a.slater); tally["P-fallback"] += 1
        else:                                       # EAV unscaled
            vals.append(p["value"]); tally["EAV"] += 1
    IP.write(raw, params, vals, a.out)
    print(f"wrote {a.out}  per-family: ZETA x{a.zeta} F x{a.fdir} G x{a.gexch} "
          f"CI x{a.ci}; fallback x{a.slater}")
    print(f"  counts: {tally}")


if __name__ == "__main__":
    main()
