#!/usr/bin/env python3
"""Canonical diagnostics for the Mg-I-style fit: ONE validated path for the level
and parameter comparisons that were being hand-rolled in throwaway snippets (and
getting wrong -- wrong eigenvector slot, wrong unit, silent empty match, median-
vs-ground offset confusion). Everything here reuses the SAME identity matching the
fit uses (gf_fit._levels_by_identity) and the SAME physical param labels
(param_labels), and ASSERTS on the failure modes instead of returning garbage.

Commands:
  levels   --outg11 F --nist N [--max-energy E] [--anchor median|ground]
           Per-level residuals (E_calc - E_obs), offset-removed, by config/term/J,
           with the worst offenders and per-term-family means. The headline RMS
           uses the SAME convention as the fit.
  level    --outg11 F --nist N --config C --term T [--J J]
           One level: E_calc, E_obs, residual -- and the eigenvector purity, so
           you can SEE whether a label is trustworthy (the 1D/3D J=2 case).
  splitting --outg11 F --config C --terms 'T1,T2'
           The computed splitting between two terms of a config (e.g. 3p2 1S vs 3P)
           vs the observed -- the right diagnostic for a Slater integral.
  param    --ing11 I --outgine O --config C --param P [--vs VALUE]
           Physical value (cm^-1) of a named integral (F2/G1/EAV/ZETA), correctly
           slotted, optionally compared to a reference value (e.g. Bob's).

Import the functions for use in analysis instead of re-deriving:
  level_residuals(outg11, nist, max_energy=None, anchor='median') -> list of dicts
  eigenvector_purity(outg11, config, term, J) -> (purity, [(label, frac), ...])
"""
import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import numpy as np
import make_report as R
from parse_cowan import parse_compositions, _BASIS_ROW, _CHUNK_HDR, _NUM
import param_labels as PL


# ----------------------------- level residuals ------------------------------

def _levels_by_identity(outg11):
    """{(cfgkey, termkey, Jstr) -> {E_calc, config, term, J}} for the OUTG11,
    using each level's dominant-eigenvector identity -- the SAME match gf_fit uses.
    On an identity collision, keep the lowest-energy level (deterministic)."""
    comp = parse_compositions(outg11)
    by = {}
    for (par, J), levs in comp.items():
        for L in levs:
            k = R._level_idkey(L["config"], L["term"], J)
            if None in k:
                continue
            if k not in by or L["E_calc"] < by[k]["E_calc"]:
                by[k] = {"E_calc": L["E_calc"], "config": L["config"],
                         "term": L["term"], "J": J}
    return by


def level_residuals(outg11, nist_path, max_energy=None, anchor="median"):
    """Per-level (E_calc - E_obs) over levels matched by identity to NIST. Offset
    removed by `anchor`: 'median' (minimizes RMS; the fit's convention) or
    'ground' (the lowest-energy matched level reads 0). Returns a list of dicts
    {cfg, term, J, E_obs, E_calc, resid}, residual already offset-removed.
    RAISES if zero levels match (the silent-empty-match failure)."""
    by = _levels_by_identity(outg11)
    nist = R.load_nist(nist_path)
    nlev = {}
    for n in nist:
        k = R._level_idkey(n["config"], n["term"], n["J"])
        if None not in k:
            nlev.setdefault(k, n["E_obs"])
    rows = []
    for k, eobs in nlev.items():
        if max_energy is not None and eobs > max_energy:
            continue
        L = by.get(k)
        if L is None:
            continue
        rows.append({"cfg": k[0], "term": k[1], "J": k[2],
                     "E_obs": eobs, "E_calc": L["E_calc"],
                     "resid_raw": L["E_calc"] - eobs})
    if not rows:
        raise ValueError("inspect.level_residuals: ZERO levels matched NIST by "
                         "identity -- check the OUTG11/NIST pairing (this is the "
                         "silent-empty-match failure, now loud).")
    draw = np.array([r["resid_raw"] for r in rows])
    if anchor == "median":
        off = float(np.median(draw))
    elif anchor == "ground":
        off = min(rows, key=lambda r: r["E_obs"])["resid_raw"]
    else:
        raise ValueError(f"unknown anchor {anchor!r} (use median|ground)")
    for r in rows:
        r["resid"] = r["resid_raw"] - off
    rows.sort(key=lambda r: r["E_obs"])
    return rows, off


# ----------------------------- eigenvector purity ---------------------------

def eigenvector_purity(outg11, config, term, Jstr):
    """(purity, top_components) for the level whose dominant identity matches
    (config, term, J). purity = |dominant coeff|^2; top_components = up to 4
    (label, frac) pairs. Lets you SEE whether a 1D/3D-type label is trustworthy
    before trusting its residual. Returns (None, []) if not found."""
    want = (R._cfgkey(config), R._termkey(term), R._Jkey(str(Jstr)))
    lines = open(outg11, errors="replace").readlines()
    n = len(lines)
    i = 0
    while i < n:
        m = re.search(r"EIGENVALUES\s*\(J=\s*([-\d.]+)\)", lines[i])
        if not m:
            i += 1
            continue
        Jval = "%g" % float(m.group(1))
        if R._Jkey(Jval) != want[2]:
            i += 1
            continue
        # eigenvalues
        evs = []
        j = i + 1
        while j < n and not lines[j].strip():
            j += 1
        while j < n:
            s = lines[j]
            if any(t in s for t in ("CONFIG. NO.", "EIGENVECTORS", "G-VALUES")):
                break
            nums = _NUM.findall(s)
            if nums:
                evs.extend(float(x) for x in nums)
            elif not s.strip() and evs:
                break
            j += 1
        while j < n and not ("EIGENVECTORS" in lines[j] and "LS" in lines[j]):
            if re.search(r"EIGENVALUES\s*\(J=", lines[j]):
                break
            j += 1
        if j >= n or "EIGENVECTORS" not in lines[j]:
            i = j
            continue
        coeffs, label = {}, {}
        col_base = cw = 0
        k = j + 1
        while k < n:
            s = lines[k].rstrip("\n")
            if "PURITY" in s or re.search(r"EIGENVALUES\s*\(J=", s) \
                    or ("EIGENVECTORS" in s and "JJ" in s):
                break
            mr = _BASIS_ROW.match(s)
            if mr:
                cfg, tm, bk = mr.group(1), mr.group(3), int(mr.group(4))
                cs = [float(x) for x in _NUM.findall(mr.group(5))]
                row = coeffs.setdefault(bk, [])
                if len(row) < col_base:
                    row.extend([0.0] * (col_base - len(row)))
                row.extend(cs)
                label[bk] = (cfg, tm)
                cw = max(cw, len(cs))
            elif _CHUNK_HDR.match(s):
                col_base += cw
                cw = 0
            k += 1
        # for each eigenvector column, dominant label; find the one matching want
        for col in range(len(evs)):
            comps = []
            for bk, cs in coeffs.items():
                if col < len(cs):
                    comps.append((cs[col] ** 2, label[bk]))
            if not comps:
                continue
            comps.sort(reverse=True)
            dom = comps[0][1]
            if (R._cfgkey(dom[0]), R._termkey(dom[1])) == (want[0], want[1]):
                top = [(f"{lab[0]} {lab[1]}", round(frac, 3))
                       for frac, lab in comps[:4]]
                return round(comps[0][0], 3), top
        i = j
    return None, []


# ----------------------------- term splitting -------------------------------

def term_splitting(outg11, config, term_a, term_b):
    """Computed energy splitting E(term_b) - E(term_a) within a config (lowest-J
    representative of each), e.g. 3p2 1S - 3P. RAISES if a term is absent."""
    comp = parse_compositions(outg11)
    found = {}
    for (par, J), levs in comp.items():
        for L in levs:
            if R._cfgkey(L.get("config", "")) == R._cfgkey(config):
                t = R._termkey(L.get("term", ""))
                found.setdefault(t, []).append(L["E_calc"])
    ta, tb = R._termkey(term_a), R._termkey(term_b)
    if ta not in found or tb not in found:
        raise ValueError(f"term_splitting: {config} missing "
                         f"{ta if ta not in found else tb} "
                         f"(have {sorted(found)})")
    ea, eb = min(found[ta]), min(found[tb])
    return eb - ea, ea, eb


# --------------------------------- CLI --------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("levels")
    p.add_argument("--outg11", required=True)
    p.add_argument("--nist", required=True)
    p.add_argument("--max-energy", type=float, default=None)
    p.add_argument("--anchor", choices=["median", "ground"], default="median")

    p = sub.add_parser("level")
    p.add_argument("--outg11", required=True)
    p.add_argument("--nist", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--term", required=True)
    p.add_argument("--J", default=None)

    p = sub.add_parser("splitting")
    p.add_argument("--outg11", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--terms", required=True, help="'T1,T2' e.g. '3P,1S'")
    p.add_argument("--observed", default=None,
                   help="optional 'E_a,E_b' observed energies for comparison.")

    p = sub.add_parser("param")
    p.add_argument("--ing11", required=True)
    p.add_argument("--outgine", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--param", required=True, help="F2 / G1 / EAV / ZETA ...")
    p.add_argument("--vs", type=float, default=None, help="reference value cm^-1.")

    a = ap.parse_args()

    if a.cmd == "levels":
        rows, off = level_residuals(a.outg11, a.nist, a.max_energy, a.anchor)
        res = np.array([r["resid"] for r in rows])
        print(f"{len(rows)} levels matched | anchor={a.anchor} (offset {off:+.1f}) "
              f"| RMS={np.sqrt((res**2).mean()):.1f}  median|d|={np.median(np.abs(res)):.1f}  "
              f"max|d|={np.abs(res).max():.1f}")
        worst = sorted(rows, key=lambda r: -abs(r["resid"]))[:12]
        print("worst:")
        for r in worst:
            print(f"  {r['resid']:+8.1f}  {r['cfg']:8} {r['term']:4} J={r['J']}")
        # per-term-family means
        from collections import defaultdict
        bt = defaultdict(list)
        for r in rows:
            bt[r["term"]].append(r["resid"])
        print("by term (|mean|>15):")
        for t, vs in sorted(bt.items(), key=lambda x: -abs(np.mean(x[1]))):
            if abs(np.mean(vs)) > 15:
                print(f"  {t:4} n={len(vs):2} mean={np.mean(vs):+7.1f}")

    elif a.cmd == "level":
        rows, off = level_residuals(a.outg11, a.nist)
        hits = [r for r in rows if R._cfgkey(r["cfg"]) == R._cfgkey(a.config)
                and R._termkey(r["term"]) == R._termkey(a.term)
                and (a.J is None or r["J"] == R._Jkey(a.J))]
        for r in hits:
            pur, top = eigenvector_purity(a.outg11, a.config, a.term, r["J"])
            print(f"{r['cfg']} {r['term']} J={r['J']}: E_calc={r['E_calc']:.1f} "
                  f"E_obs={r['E_obs']:.1f} resid={r['resid']:+.1f} "
                  f"purity={pur} top={top}")
        if not hits:
            print("no matching level")

    elif a.cmd == "splitting":
        ta, tb = a.terms.split(",")
        d, ea, eb = term_splitting(a.outg11, a.config, ta, tb)
        msg = f"{a.config}: {tb}-{ta} = {d:+.0f} cm^-1  (E({ta})={ea:.0f}, E({tb})={eb:.0f})"
        if a.observed:
            oa, ob = (float(x) for x in a.observed.split(","))
            msg += f"  | observed {tb}-{ta} = {ob-oa:+.0f}  -> off {d-(ob-oa):+.0f}"
        print(msg)

    elif a.cmd == "param":
        v, info = PL.param_value_cm1(a.ing11, a.outgine, a.config, a.param)
        msg = f"{a.config} {info['phys']} = {v:.1f} cm^-1  (slot {info['key'].split('|')[-1]})"
        if a.vs is not None:
            msg += f"  | ref {a.vs:.1f}  ratio {v/a.vs:.3f}"
        print(msg)


if __name__ == "__main__":
    main()
