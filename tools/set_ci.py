#!/usr/bin/env python3
"""Set a specific CI (configuration-interaction R^k) integral in an ING11 to a
target value in cm^-1 -- e.g. to freeze it at Bob Kurucz's curated value instead
of our scaled-HF. CI physical cm^-1 = stored * 1000 (verified vs Bob's HF: our
3s2-3p2 HF stored 29.4047 = 29405 cm^-1 ~ Bob HF 32465; NOT *SCALE_CI=1e4).

Usage:
  tools/set_ci.py --in work/mg1_full/ING11.scaled --out work/mg1_full/ING11.bobci \
      --pair 3s2-3p2 --value 25972
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ing11_params as IP

CI_TO_CM1 = 1000.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pair", required=True, help="config pair, e.g. '3s2-3p2'.")
    ap.add_argument("--value", type=float, required=True, help="target value, cm^-1.")
    a = ap.parse_args()
    want = tuple(sorted(x.strip() for x in a.pair.split("-")))
    raw, params = IP.parse(a.inp)
    vals = [p["value"] for p in params]
    hits = []
    for i, p in enumerate(params):
        if p["kind"] != "CI":
            continue
        pr = tuple(sorted(x.strip() for x in p["key"].split("|")[0].split("-")))
        if pr == want:
            old = vals[i] * CI_TO_CM1
            vals[i] = a.value / CI_TO_CM1
            hits.append((p["key"], old, a.value))
    if not hits:
        raise SystemExit(f"no CI integral found for pair {a.pair!r}")
    IP.write(raw, params, vals, a.out)
    for key, old, new in hits:
        print(f"  {key}: {old:.1f} -> {new:.1f} cm^-1")
    print(f"wrote {a.out}  ({len(hits)} CI slot(s) set)")


if __name__ == "__main__":
    main()
