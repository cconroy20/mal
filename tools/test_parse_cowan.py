#!/usr/bin/env python3
"""Regression tests for tools/parse_cowan.py.

Run standalone (no pytest needed):   python3 tools/test_parse_cowan.py
Or under pytest:                     pytest tools/test_parse_cowan.py

THE BUG THIS GUARDS AGAINST (2026-06-30): parse_compositions read RCG's
EIGENVALUES block but terminated on the first BLANK LINE. RCG prints eigenvalues
11-per-line separated by blank lines, so the parser kept only the first ~11 of
each J-block and silently dropped the rest -- for the full Mg I basis, 192 of
1105 levels parsed and 64 of 122 configs seen. Every high-lying doubly-excited
config (3p3d, 3p4s, 3d.nl) vanished, which is exactly the data the fit needs.

THE INVARIANT (self-contained, needs no external reference): within each
EIGENVALUES(J=) block of an OUTG11, the number of eigenvalues MUST equal the
number of CONFIG. NO. entries -- both are the J-block matrix dimension. If the
eigenvalue read stops early, the counts diverge. This catches the regression at
its source, independent of any downstream fit or basis size.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import parse_cowan as P


def _block_counts(path):
    """For each EIGENVALUES(J=) block: (Jlabel, n_eigenvalues, n_config_no)."""
    lines = open(path, errors="replace").readlines()
    n = len(lines)
    i = 0
    out = []
    while i < n:
        m = re.search(r"EIGENVALUES\s*\(J=\s*([-\d.]+)\)", lines[i])
        if not m:
            i += 1
            continue
        j = i + 1
        evs = []
        while j < n and not lines[j].strip():
            j += 1
        while j < n:
            s = lines[j]
            if any(t in s for t in ("CONFIG. NO.", "EIGENVECTORS", "G-VALUES")) \
                    or re.search(r"EIGENVALUES\s*\(J=", s):
                break
            evs += P._NUM.findall(s)
            j += 1
        cfgnos = []
        if j < n and "CONFIG. NO." in lines[j]:
            k = j + 1
            while k < n:
                s = lines[k]
                if any(t in s for t in ("G-VALUES", "EIGENVECTORS")) \
                        or re.search(r"EIGENVALUES\s*\(J=", s):
                    break
                cfgnos += re.findall(r"\b\d+\b", s)
                k += 1
        if cfgnos:
            out.append((m.group(1), len(evs), len(cfgnos)))
        i = j
    return out


# reference OUTG11s that must always parse fully
FULL = os.path.join(ROOT, "work", "mg1_full", "OUTG11.abinitio")   # 122 configs
NINE = os.path.join(ROOT, "work", "mg1", "OUTG11")                 # 9 configs


def test_eigenvalue_count_matches_config_no():
    """Every J-block: n_eigenvalues == n_CONFIG.NO (the matrix dimension). This is
    the direct guard on the blank-line-truncation bug."""
    for path in (FULL, NINE):
        if not os.path.exists(path):
            continue
        for Jlab, nev, ncfg in _block_counts(path):
            assert nev == ncfg, (
                f"{os.path.basename(path)} J={Jlab}: {nev} eigenvalues but "
                f"{ncfg} CONFIG.NO entries -- eigenvalue read truncated early "
                f"(the blank-line bug). Levels are being silently dropped.")


def test_full_basis_parses_all_configs():
    """The 122-config Mg I run must expose all 122 configs and its full level set
    (was 64 configs / 192 levels under the bug). Guards the downstream symptom."""
    if not os.path.exists(FULL):
        return
    comp = P.parse_compositions(FULL)
    total = sum(len(v) for v in comp.values())
    cfgs = {L["config"].replace("Mg I", "").strip()
            for levs in comp.values() for L in levs}
    assert len(cfgs) >= 120, f"only {len(cfgs)} configs parsed (expected 122)"
    assert total >= 1000, f"only {total} levels parsed (expected ~1105)"
    # the doubly-excited valence configs the bug dropped must be present
    flat = {c.replace(" ", "").replace(".", "") for c in cfgs}
    for want in ("3p3d", "3p4s"):
        assert any(c.startswith(want) for c in flat), f"{want} missing from parse"


def test_nine_config_unchanged():
    """Small bases (single eigenvalue row, no blank-line wrap) are unaffected: the
    9-config Mg I deck still parses to exactly 30 levels across 9 configs."""
    if not os.path.exists(NINE):
        return
    comp = P.parse_compositions(NINE)
    total = sum(len(v) for v in comp.values())
    cfgs = {L["config"].replace("Mg I", "").strip()
            for levs in comp.values() for L in levs}
    assert total == 30, f"9-config parse changed: {total} levels (expected 30)"
    assert len(cfgs) == 9, f"9-config parse changed: {len(cfgs)} configs (expected 9)"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {t.__name__}\n      {e}")
    print(f"\n{len(tests)-failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
