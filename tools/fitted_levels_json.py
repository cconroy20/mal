#!/usr/bin/env python3
"""Emit a fitted-levels JSON for make_report --fitted-levels-json from a gf_fit
(Python optimizer) result. gf_fit has no RCE LEVELS1 file, so we reconstruct the
{config, term, J, parity, E_obs, E_fit} rows the report's level/residual pages
want directly from the FITTED OUTG11, using the same eigenvector-identity match
the fit itself uses: each NIST level is paired with the computed level of the
same (cfgkey, termkey, Jstr) identity; E_fit = that computed level's energy.

Usage:
  python3 tools/fitted_levels_json.py --outg11 work/mg1_full/OUTG11.energyonly \
      --nist data/nist/MgI_levels.tsv --max-energy 61671 \
      --out work/mg1_full/fitted_levels.json
"""
import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import make_report as R
from parse_cowan import parse_compositions, _parity_from_config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outg11", required=True, help="fitted-parameter OUTG11.")
    ap.add_argument("--nist", required=True)
    ap.add_argument("--max-energy", type=float, default=None,
                    help="drop NIST levels above this energy (cm^-1).")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    nist = R.load_nist(a.nist)
    # NIST level energy by robust identity (cfgkey, termkey, Jstr)
    nist_lev = {}
    nist_lbl = {}
    for n in nist:
        k = R._level_idkey(n["config"], n["term"], n["J"])
        if k[0] is None:
            continue
        nist_lev.setdefault(k, n["E_obs"])
        nist_lbl.setdefault(k, n)            # keep an authoritative label

    # computed level by identity (dominant eigenvector basis state + block J)
    comp = parse_compositions(a.outg11)
    by_id = {}
    for (par, J), levs in comp.items():
        for L in levs:
            k = R._level_idkey(L["config"], L["term"], J)
            if None in k:
                continue
            # collision: keep the one closest to E_obs (if known)
            prev = by_id.get(k)
            if prev is None:
                by_id[k] = L
            else:
                eobs = nist_lev.get(k)
                if eobs is not None and \
                        abs(L["E_calc"] - eobs) < abs(prev["E_calc"] - eobs):
                    by_id[k] = L

    rows = []
    for k, eobs in nist_lev.items():
        if a.max_energy is not None and eobs > a.max_energy:
            continue
        L = by_id.get(k)
        if L is None:
            continue
        n = nist_lbl[k]
        rows.append({"config": n["config"], "term": n["term"], "J": str(n["J"]),
                     "parity": _parity_from_config(n["config"]),
                     "E_obs": float(eobs), "E_fit": float(L["E_calc"])})

    with open(a.out, "w") as f:
        json.dump(rows, f)
    print(f"wrote {a.out}  ({len(rows)} fitted levels)")


if __name__ == "__main__":
    main()
