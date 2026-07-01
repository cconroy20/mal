#!/usr/bin/env python3
"""Parse Bob Kurucz's FITTED radial parameters from an RCE least-squares log
(c<xxyy>{e,o}z.log). Returns the free AND fixed parameters with their fitted
value, HF reference, screening scale, kind, and config label.

Purpose: the "full RCE reproduction" test -- transcribe Bob's fitted deck into our
ING11, run OUR RCG, and see whether we reproduce his level residuals (isolates
his-parameters vs our-forward-model as the source of the 3s.nd 1D gap).

Bob's parameter-line format (fixed-ish columns), e.g.:
    1    1  EAV        26194.7    3.0   20894.6             1.000  0    2600.0  0   3s3p
    3    1  G1(12)     23423.5   10.0   29818.4  FIXEDHF    0.600  0       0.0  0
   39    I  ZETA 4         0.0    0.0       0.0  FIXEDHF    1.000  0       0.0  0
  fields: idx, code, KIND[(group)], VALUE, UNCERT, HF_REF, [FIXEDHF], SCALE, ...,
          trailing CONFIG label (may be blank for CI/zeta with no single config).
FREE params have no 'FIXEDHF' token AND a nonzero UNCERT; fixed have 'FIXEDHF'.

This does NOT parse the eigenvalue/level table (those lines start with a parity
digit 0/1 and carry an 'R' marker + doubled term labels) -- the KIND regex only
matches genuine parameter kinds (EAV|ZETA|F k|G k|R k).
"""
import re
import sys

# A genuine parameter KIND: EAV, ZETA [n], F/G/R with an integer k and optional
# (group). Anchored so numeric level-table junk can't match.
_KIND = re.compile(r"^(EAV|ZETA(?:\s+\d+)?|[FGR]\d+(?:\(\d+\))?)$")
# One parameter line. Groups: idx, code, kind, value, uncert, hf_ref, rest(has
# FIXEDHF? + scale + trailing config).
_LINE = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+"                       # idx, code
    r"([A-Z]+\d*(?:\s+\d+)?(?:\(\d+\))?)\s+"      # KIND (EAV / ZETA n / G1(13))
    r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)" # value, uncert, hf_ref
    r"(.*)$")                                     # rest

# The 2nd column ("code") is Bob's CONFIG INDEX: every parameter of one config
# (its EAV, ZETA, Gk, Fk) shares it. The EAV line for that code carries the config
# LABEL (e.g. '3s3p'); the Gk/Fk/ZETA lines for the same code do NOT. So we group
# by code and inherit the EAV's label -- no need to decode the occupation-vector
# header. NB: code is per-PARITY (each parity block restarts at 1), so a parse must
# not merge codes across the even/odd logs.


def parse(path):
    """Return a list of dicts {code, kind, group, value, uncert, hf_ref, scale,
    fixed, config, raw} for every parameter line in an RCE log. `config` is filled
    from the line's own trailing label when present (EAV lines) and otherwise
    INHERITED from the EAV of the same config `code` (Gk/Fk/ZETA lines)."""
    out = []
    for ln in open(path, errors="replace"):
        m = _LINE.match(ln)
        if not m:
            continue
        kind_raw = m.group(3).strip()
        # normalize 'ZETA 4' -> kind ZETA; 'G1(13)' -> kind G1, group 13
        knorm = re.sub(r"\s+\d+$", "", kind_raw)          # drop ZETA's n
        if not _KIND.match(kind_raw):
            continue
        gm = re.search(r"\((\d+)\)", kind_raw)
        group = int(gm.group(1)) if gm else None
        # split the letter kind from the Slater rank k: 'G1(13)' -> kind 'G', k 1;
        # 'F2(28)' -> kind 'F', k 2; 'EAV'/'ZETA' -> kind as-is, k None.
        km = re.match(r"([A-Z]+)(\d+)?", knorm)
        base_kind = km.group(1)
        krank = int(km.group(2)) if km.group(2) else None
        value, uncert, hf_ref = float(m.group(4)), float(m.group(5)), float(m.group(6))
        rest = m.group(7)
        fixed = "FIXEDHF" in rest
        sm = re.search(r"(-?\d+\.\d+)", rest.replace("FIXEDHF", ""))
        scale = float(sm.group(1)) if sm else None
        # trailing config label: last run of config-ish tokens (e.g. '3s3p', '3s3d
        # -3p2' for a CI). Take everything after the last numeric flag column.
        cm = re.search(r"([0-9][spdfghik][0-9a-z.\- ]*?)\s*$", rest)
        config = cm.group(1).strip() if cm else ""
        out.append(dict(code=m.group(2), kind=base_kind, k=krank, group=group,
                        value=value, uncert=uncert, hf_ref=hf_ref, scale=scale,
                        fixed=fixed, config=config, raw=ln.rstrip()))
    # INHERIT the config label from each code's EAV line to its Gk/Fk/ZETA lines.
    # (A CI line's config is a "cfgA - cfgB" pair carried on the line itself, so it
    # keeps its own label and is not overwritten.)
    label_of_code = {p["code"]: p["config"] for p in out
                     if p["kind"] == "EAV" and p["config"]}
    for p in out:
        if not p["config"] and p["kind"] in ("G", "F", "ZETA"):
            p["config"] = label_of_code.get(p["code"], "")
    return out


def _selftest():
    """Validate against hand-verified values from c1200oz.log."""
    ps = parse("kurucz_ref/1200/c1200oz.log")
    by = {}
    for p in ps:
        by.setdefault((p["kind"], p["config"]), p)
    checks = [
        (("EAV", "3s3p"), 26194.7, False),
        (("G1", ""), 23423.5, False),   # G1(12), no config label
    ]
    # find the 3s3p EAV and the first free G1
    eav = [p for p in ps if p["kind"] == "EAV" and p["config"] == "3s3p"]
    g1 = [p for p in ps if p["kind"] == "G" and p["k"] == 1 and not p["fixed"]]
    ok = True
    if eav and abs(eav[0]["value"] - 26194.7) < 0.1:
        print(f"PASS  3s3p EAV = {eav[0]['value']} (expected 26194.7)")
    else:
        print(f"FAIL  3s3p EAV: {[e['value'] for e in eav]}"); ok = False
    if g1 and abs(g1[0]["value"] - 23423.5) < 0.1:
        print(f"PASS  first free G1 = {g1[0]['value']} (expected 23423.5, G1(12))")
    else:
        print(f"FAIL  first free G1: {[g['value'] for g in g1[:3]]}"); ok = False
    free = [p for p in ps if not p["fixed"]]
    print(f"total params: {len(ps)}  free: {len(free)}  "
          f"(EAV free: {sum(1 for p in free if p['kind']=='EAV')})")
    return ok


if __name__ == "__main__":
    import os
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.exit(0 if _selftest() else 1)
