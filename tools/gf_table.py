#!/usr/bin/env python3
"""Per-line gf comparison table: fitted (and ab initio) Cowan log gf vs NIST.

Lines are matched to NIST by the EIGENVECTOR-COMPOSITION identity of both end
levels (config + term + J), the same robust identity the RCE level fit is built
on -- NOT by a term-pair bucket + nearest wavelength, which collides Rydberg
series members (e.g. 3p^2 1S vs 3s5s 1S) and manufactures bogus residuals.

A computed line is matched to a NIST line only when BOTH the lower and upper
levels agree in normalized (config, term, J). Computed lines whose endpoints
have no NIST counterpart (e.g. 3p^2 perturber lines absent from our NIST set)
are reported as UNMATCHED rather than mis-paired.

Usage:
  tools/gf_table.py --abinitio work/mg1/OUTG11.abinitio \
                    --fitted   work/mg1/OUTG11.fitted \
                    --nist-lines data/nist/MgI_lines.tsv
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import make_report as R
from parse_cowan import identify_lines


def matched_table(outg11_path, nist_lines):
    """Identify computed lines by eigenvector composition, then match to NIST by
    (lower-level key, upper-level key) -- the SAME identity match_gf_by_identity
    uses, but keeping per-line labels/accuracy for the table. Returns
    (rows, n_unmatched)."""
    buckets = {}
    for d in nist_lines:
        k = R._line_idkey(d["conf_i"], d["term_i"], d["J_i"],
                          d["conf_k"], d["term_k"], d["J_k"])
        buckets.setdefault(k, []).append(d)

    lines = identify_lines(outg11_path)
    rows, unmatched = [], 0
    used = set()
    for d in lines:
        if not (d["config_low"] and d["config_up"]):
            unmatched += 1
            continue
        k = R._line_idkey(d["config_low"], d["term_id_low"], d["J_low"],
                          d["config_up"], d["term_id_up"], d["J_up"])
        cands = [c for c in buckets.get(k, []) if id(c) not in used]
        if not cands:
            unmatched += 1
            continue
        # within an exact identity bucket, pick the closest wavelength (handles
        # the rare J-resolved fine-structure multiplet sharing one identity key)
        n = min(cands, key=lambda c: abs(c["lambda_A"] - d["lambda_A"]))
        used.add(id(n))
        rows.append({
            "lam": d["lambda_A"], "nist_lam": n["lambda_A"],
            "comp": d["loggf"], "nist": n["loggf"],
            "d": d["loggf"] - n["loggf"], "acc": n.get("acc", ""),
            "low": f"{R._cfgkey(d['config_low'])} {R._termkey(d['term_id_low'])}",
            "up": f"{R._cfgkey(d['config_up'])} {R._termkey(d['term_id_up'])}",
        })
    return rows, unmatched


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abinitio")
    ap.add_argument("--fitted")
    ap.add_argument("--nist-lines", required=True)
    a = ap.parse_args()

    nist_lines = R.load_nist_lines(a.nist_lines)
    series = {}
    meta = {}
    if a.abinitio:
        rows, un = matched_table(a.abinitio, nist_lines)
        series["abinitio"] = {r["lam"]: r for r in rows}
        meta["abinitio"] = un
    if a.fitted:
        rows, un = matched_table(a.fitted, nist_lines)
        series["fitted"] = {r["lam"]: r for r in rows}
        meta["fitted"] = un

    fit = series.get("fitted", {})
    ab = series.get("abinitio", {})
    keys = sorted(fit or ab)
    print(f"{'lam_A':>10} {'lower':>10} {'upper':>10} {'NIST':>7} "
          f"{'fit':>7} {'dfit':>7} {'abin':>7} {'dab':>7} {'acc':>5}")
    print("-" * 80)
    for lam in keys:
        rf = fit.get(lam); ra = ab.get(lam)
        base = rf or ra
        nist = base["nist"]
        fcol = f"{rf['comp']:+7.3f}" if rf else "   -   "
        dfit = f"{rf['d']:+7.3f}" if rf else "   -   "
        acol = f"{ra['comp']:+7.3f}" if ra else "   -   "
        dab = f"{ra['d']:+7.3f}" if ra else "   -   "
        print(f"{lam:10.2f} {base['low']:>10} {base['up']:>10} "
              f"{nist:+7.3f} {fcol} {dfit} {acol} {dab} {base['acc']:>5}")

    # NIST accuracy class -> approx fractional uncertainty in gf (ASD legend);
    # used to (a) restrict the RMS to reliably-measured lines and (b) weight by
    # measurement quality, so weak/poorly-known lines don't dominate the metric.
    ACC_DEX = {"AAA": 0.013, "AA": 0.013, "A+": 0.013, "A": 0.022, "B+": 0.043,
               "B": 0.087, "C+": 0.13, "C": 0.22, "D+": 0.30, "D": 0.43,
               "E": 0.70}

    def acc_dex(a):
        return ACC_DEX.get(a.strip(), 0.5)

    for name, s in (("ab initio", ab), ("fitted", fit)):
        if not s:
            continue
        rows = list(s.values())
        dd = np.array([r["d"] for r in rows])
        key = "abinitio" if name == "ab initio" else "fitted"
        # strong + well-measured subset: log gf >= -1 and NIST class A/B
        good = [r for r in rows
                if r["nist"] >= -1.0 and r["acc"].strip()[:1] in ("A", "B")]
        gd = np.array([r["d"] for r in good]) if good else np.array([])
        # inverse-variance weighting by NIST accuracy
        w = np.array([1.0 / acc_dex(r["acc"]) ** 2 for r in rows])
        wrms = np.sqrt(np.sum(w * dd ** 2) / np.sum(w))
        print(f"\n{name:>10}: N={len(dd)}  RMS={np.sqrt(np.mean(dd**2)):.3f}  "
              f"mean={dd.mean():+.3f}  max|d|={np.abs(dd).max():.3f}  "
              f"(unmatched dropped: {meta.get(key, 0)})")
        print(f"{'':>10}  acc-weighted RMS={wrms:.3f}   "
              f"strong+A/B (loggf>=-1): N={len(gd)} "
              f"RMS={np.sqrt(np.mean(gd**2)):.3f}" if len(gd) else
              f"{'':>10}  acc-weighted RMS={wrms:.3f}")


if __name__ == "__main__":
    main()
