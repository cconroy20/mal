#!/usr/bin/env python3
"""Extract a structured free/fixed RADIAL-PARAMETER profile from Bob Kurucz's RCE
fit logs (c<xxyy>{e,o}*.log) across many species, so his fitting choices can be
compared and turned into general rules.

Each parameter row in a c*.log looks like:
   <idx> <grp>  <NAME>  <value>  <sigma/step>  <hf-value> [FIXEDHF] ... <config>
A parameter is FREE when it is NOT marked FIXEDHF/FIXED (its sigma/step is the
fit uncertainty). NAME is EAV / ZETA n / Fk(ij) / Gk(ij) / Rk(...) / ALPHA / BETA
/ T... ; (ij) are 1-based orbital indices within that block's orbital list.

Usage:
  tools/survey_fits.py 1100 1200 2000 ...     # prints a per-species summary
  tools/survey_fits.py --json 1100 1200 ...   # machine-readable for analysis
"""
import json
import os
import re
import sys
from collections import Counter

REF = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "kurucz_ref")

_ROW = re.compile(
    r"^\s*\d+\s+[0-9A-Za-z]+\s+"
    r"(EAV|ZETA\s*\d+|F\d\(\d+\)|G\d\(\d+\)|R\d\([^)]*\)|ALPHA|BETA|T\d*)\s+"
    r"(-?\d+\.\d+)\s+(\d+\.\d+)\b")
# orbital occupations from the block header: "... 1  d 7  s 0  p 0 ..."
_HDR = re.compile(r"NPAR\s*(\d+)\s+NJ\s*(\d+)\s+FIRSTJ\s+([\d.]+)\s+(.*)")


def _kind(name):
    """Coarse parameter family for tallying."""
    if name.startswith("EAV"):
        return "EAV"
    if name.startswith("ZETA"):
        return "ZETA"
    if name.startswith("ALPHA") or name.startswith("BETA"):
        return "ALPHA/BETA"
    if name.startswith("T"):
        return "Trees-T"
    if re.match(r"F\d", name):
        return "F^k(direct)"
    if re.match(r"G\d", name):
        return "G^k(exchange)"
    if re.match(r"R\d", name):
        return "R^k(CI)"
    return "other"


def parse_log(path):
    """Return (header_dict, list_of_param_dicts) for one c*.log."""
    header = {}
    params = []
    for ln in open(path, errors="replace"):
        if not header:
            mh = _HDR.search(ln)
            if mh:
                header = {"npar": int(mh.group(1)), "nj": int(mh.group(2)),
                          "firstj": mh.group(3),
                          "ground_occ": " ".join(mh.group(4).split())}
        m = _ROW.match(ln)
        if not m:
            continue
        name = re.sub(r"\s+", "", m.group(1))
        fixed = ("FIXEDHF" in ln) or re.search(r"\bFIXED\b", ln) is not None
        params.append({"name": name, "kind": _kind(name),
                       "value": float(m.group(2)), "sigma": float(m.group(3)),
                       "fixed": fixed})
    return header, params


def species_profile(xxyy):
    """Aggregate even+odd logs for one species into a profile."""
    d = os.path.join(REF, xxyy)
    logs = sorted(f for f in os.listdir(d)
                  if re.match(rf"c{xxyy}[eo].*\.log$", f)) if os.path.isdir(d) \
        else []
    # prefer the plain e/o; else the first e and first o variant
    chosen = {}
    for f in logs:
        par = "e" if re.match(rf"c{xxyy}e", f) else "o"
        chosen.setdefault(par, f)  # first wins (sorted: plain before y/z)
    prof = {"xxyy": xxyy, "logs": list(chosen.values()),
            "free": Counter(), "fixed": Counter(), "header": {}}
    allfree = []
    for f in chosen.values():
        hdr, params = parse_log(os.path.join(d, f))
        if hdr and not prof["header"]:
            prof["header"] = hdr
        for p in params:
            (prof["free"] if not p["fixed"] else prof["fixed"])[p["kind"]] += 1
            if not p["fixed"]:
                allfree.append(p["name"])
    prof["n_free"] = sum(prof["free"].values())
    prof["n_fixed"] = sum(prof["fixed"].values())
    prof["free_names"] = dict(Counter(allfree).most_common())
    return prof


SPECIES_NAME = {
    "1100": "Na I", "1200": "Mg I", "1300": "Al I", "1400": "Si I",
    "2000": "Ca I", "2001": "Ca II", "2201": "Ti II", "2400": "Cr I",
    "2600": "Fe I", "2601": "Fe II",
}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    as_json = "--json" in sys.argv
    profs = [species_profile(x) for x in args]
    if as_json:
        print(json.dumps([{k: (dict(v) if isinstance(v, Counter) else v)
                           for k, v in p.items()} for p in profs], indent=1))
        return
    for p in profs:
        nm = SPECIES_NAME.get(p["xxyy"], p["xxyy"])
        h = p["header"]
        print(f"\n=== {nm} ({p['xxyy']})  ground={h.get('ground_occ','?')}  "
              f"NPAR~{h.get('npar','?')}  free={p['n_free']} fixed={p['n_fixed']} "
              f"({100*p['n_free']/max(1,p['n_free']+p['n_fixed']):.0f}% free) ===")
        print("  FREE by family:", dict(p["free"].most_common()))


if __name__ == "__main__":
    main()
