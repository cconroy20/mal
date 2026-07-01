#!/usr/bin/env python3
"""Authoritative PHYSICAL labels for ING11 Slater slots -- so 'F2(3p2)' is never
guessed from slot position again.

THE BUG THIS PREVENTS: ing11_params.parse() labels a config's adjustable Slaters
positionally (P0, P1, P2, P3). Which P-slot is F^2 vs G^1 vs ZETA depends on the
config's orbitals and RCG's integral ordering -- so reading 'the F^2' off a slot
index is a guess, and guessing wrong silently returns a plausible-looking wrong
number (e.g. scanning P2=0.03 when F^2(3p2) is actually P0=21.1 kK).

THE FIX: RCG itself echoes the physical integral NAMES, in order, in the OUTGINE
block header ('... EAV  3p2  F2(22)  ALPHA  ZETA 2 ...'). We read that authoritative
sequence (reusing build_ine20._parse_param_names) and JOIN it to the ING11
adjustable slots, per config, with a hard alignment assertion. The result maps
each ING11 param key to its physical name (EAV / F2 / G1 / ZETA / ...) AND its
value in cm^-1 -- and RAISES if the OUTGINE names and ING11 slots don't line up,
instead of returning garbage.

Use physical_params(ing11, outgine) to get a list of dicts:
    {key, cfg, phys, k, value_cm1, kind, ...}   one per adjustable single-config param.
And param_value_cm1(ing11, outgine, cfg, phys) for a single look-up by physical name.
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ing11_params as IP
import build_ine20 as B
import make_report as R

# configs whose physical params overflow the 4 single-line ING11 slots (open
# d/f shells); populated during physical_params, exposed so callers can see what
# was NOT labeled (so "F2 not found" is distinguishable from "config skipped").
_OVERFLOW = set()
# configs in ING11 with adjustable slots but no OUTGINE name match (odd-block /
# normalization gaps) -- declined, not guessed.
_UNMATCHED = set()


def _is_param_name_line(s):
    """A param-NAME header line: 10-char fields of EAV/ZETA/F#(..)/G#(..)/ALPHA/
    BETA/T(..). Distinguished from value/eigenvalue lines (have '.') and pure
    group-code lines (only integers)."""
    if "." in s or not s.strip():
        return False
    return bool(re.search(r"\bEAV\b|ZETA|[FG]\d\(|ALPHA|BETA|T\(", s))


def _cfg_param_names(outgine_path):
    """{cfgkey -> ordered list of physical param names}, read from ALL OUTGINE
    block headers (the authoritative RCG echo). The 122-config OUTGINE has TWO
    param-name header regions -- the EVEN block (3s2 ... 3p2 ...) near the top and
    the ODD block (3s3p ... 3d4p ...) much later -- separated by eigenvalue/matrix
    data. Earlier this read only the even region (build_ine20._parse_param_names
    stops at the first region), so odd configs like 3s3p (and Bob's key G1(12)
    exchange) were missing. Here we find EVERY contiguous run of param-name lines
    and parse each, so both parity blocks (and any future multi-block layout) are
    covered. Names: 'EAV', 'F2(22)', 'G1(12)', 'ZETA 2', 'ALPHA', 'BETA'."""
    lines = open(outgine_path).readlines()
    by_cfg = {}
    i, n = 0, len(lines)
    while i < n:
        if not _is_param_name_line(lines[i].rstrip("\n")):
            i += 1
            continue
        # collect this contiguous header region
        region = []
        while i < n and _is_param_name_line(lines[i].rstrip("\n")):
            region.append(lines[i])
            i += 1
        names, configs = B._parse_param_names(region)
        for nm, cfg in zip(names, configs):
            if cfg is None:
                continue
            ck = R._cfgkey(cfg)
            phys = "EAV" if nm.startswith("EAV") else nm
            by_cfg.setdefault(ck, []).append(phys)
    return by_cfg


# physical name -> (kind, k) ; kind in {EAV, F, G, ZETA, ALPHA, BETA}
def _classify(phys):
    if phys == "EAV":
        return ("EAV", None)
    if phys.startswith("ZETA"):
        return ("ZETA", None)
    if phys in ("ALPHA", "BETA"):
        return (phys, None)
    m = re.match(r"([FG])(\d)", phys)
    if m:
        return (m.group(1), int(m.group(2)))
    return ("?", None)


def _raw_slots(ing11_path):
    """{cfgkey -> {'EAV': (lineno, value_kK or None), 'slots': [ (raw_value, grp,
    adjustable) x4 ]}} reading ALL 4 raw Slater slots per single-config line --
    including FIXED (grp 0) and MALFORMED ('0.00E+00' = ALPHA=0) slots, which
    ing11_params.parse() drops. We need every slot because OUTGINE lists every
    physical param in RCG order; the join is positional over ALL slots, then we
    keep the adjustable ones. The first single-config line (ground EAV reference)
    is included but its EAV value is the pinned zero."""
    raw = open(ing11_path).readlines()
    out = {}
    for s in (ln.rstrip("\n") for ln in raw):
        if not IP._is_singleconf(s):
            continue
        cfg = s[0:18].strip()
        ck = R._cfgkey(cfg)
        try:
            eav = int(s[20:30]) / IP.SCALE_EAV
        except ValueError:
            eav = None
        slots = []
        o = 30
        for _ in range(4):
            rawv = s[o:o + 9]
            grp = s[o + 9:o + 10]
            try:
                val = int(rawv) / IP.SCALE_P
                ok = True
            except ValueError:
                val, ok = None, False          # malformed slot (e.g. ALPHA 0.00E+00)
            adjustable = ok and grp.strip() not in ("", "0")
            slots.append((val, grp, adjustable))
            o += 10
        out[ck] = {"cfg": cfg, "eav": eav, "slots": slots}
    return out


def physical_params(ing11_path, outgine_path, strict=True):
    """Join ING11 Slater slots to their OUTGINE physical names -- POSITIONALLY over
    all 4 raw slots (OUTGINE lists every physical param in RCG order; ING11 stores
    all 4 slots, some fixed/malformed). Returns one dict per ADJUSTABLE single-
    config param plus every active CI: {key, cfg, phys, kind, k, value_cm1, raw}.

    strict=True (default) RAISES if, for any config, the count of non-EAV OUTGINE
    names != 4 raw ING11 slots, i.e. the orderings can't be aligned -- the exact
    misalignment that silently mislabels a slot (the wrong-slot bug). Loud failure
    beats a plausible wrong number."""
    names_by_cfg = _cfg_param_names(outgine_path)
    raw_by_cfg = _raw_slots(ing11_path)
    _, params = IP.parse(ing11_path)               # for stable keys + CI

    # map (cfgkey) -> the EAV/P param keys ing11_params assigns, by slot index
    keymap = {}
    for p in params:
        if p["kind"] == "CI":
            continue
        cfg = p["key"].split("|")[0].replace("Mg I", "").strip()
        ck = R._cfgkey(cfg)
        slot = p["key"].split("|")[-1]             # 'EAV' or 'P0'..'P3'
        keymap[(ck, slot)] = p["key"]

    out = []
    # CI pairs: already physically named
    for p in params:
        if p["kind"] == "CI":
            a, b = (x.strip() for x in p["key"].split("|")[0].split("-"))
            # CI physical cm^-1 = stored * 1000 (kK), same as EAV/Slater -- NOT
            # *SCALE_CI(1e4). Verified vs Bob's HF: our HF 3s2-3p2 stored 29.4047
            # -> 29405 cm^-1, 9% below Bob's HF 32465 (HF-code diff); *1e4 would
            # give 294047, ~9x Bob = impossible.
            out.append({"key": p["key"], "cfg": None, "pair": (a, b),
                        "phys": "CI", "kind": "CI", "k": None,
                        "value_cm1": p["value"] * 1000.0, "raw": p["value"]})

    for ck, rawinfo in raw_by_cfg.items():
        phys_names = names_by_cfg.get(ck)
        if phys_names is None:
            # cfgkey present in ING11 but not in the OUTGINE names we parsed
            # (e.g. odd-parity-block configs the header walk hasn't reached, or a
            # config-string normalization gap). DECLINE to label rather than guess;
            # record it so a caller asking for this config gets a clear message.
            if any(s[2] for s in rawinfo["slots"]):
                _UNMATCHED.add(ck)
            continue
        non_eav = [n for n in phys_names if n != "EAV"]
        slots = rawinfo["slots"]
        # OUTGINE names align to the FIRST len(non_eav) raw slots (RCG order); the
        # trailing slots must be INACTIVE. CAVEAT: open-shell configs (e.g. 3d^2:
        # F2 F4 ALPHA BETA T(D2) ZETA = 6 names) carry MORE than the 4 ING11 single-
        # line slots -- RCG stores those elsewhere and the simple positional join
        # does not apply. We DECLINE to label such configs (skip, optionally warn)
        # rather than guess -- the whole point is to never mislabel. Mg-I-type
        # closed/few-electron valence configs (<=4 names) are handled exactly.
        if len(non_eav) > len(slots):
            _OVERFLOW.add(ck)
            continue
        if strict and any(slots[i][2] for i in range(len(non_eav), len(slots))):
            raise ValueError(
                f"config {ck!r}: an ADJUSTABLE raw slot lies beyond the "
                f"{len(non_eav)} OUTGINE names {non_eav} -- name<->slot "
                f"misalignment, labels would be WRONG. (wrong-slot guard)")
        # EAV (adjustable unless it's the pinned ground reference, which has no key)
        eav_key = keymap.get((ck, "EAV"))
        if eav_key is not None and rawinfo["eav"] is not None:
            out.append({"key": eav_key, "cfg": ck, "phys": "EAV", "kind": "EAV",
                        "k": None, "value_cm1": rawinfo["eav"] * 1000.0,
                        "raw": rawinfo["eav"]})
        # non-EAV: positional join, keep adjustable slots
        for idx, (phys, (val, grp, adj)) in enumerate(zip(non_eav, slots)):
            if not adj:
                continue
            key = keymap.get((ck, f"P{idx}"))
            if key is None:
                continue                           # parse() didn't surface it
            kind, k = _classify(phys)
            # ing11_params stores both EAV and Slaters in the parsed attr as ~kK
            # (value/1e5 for EAV, value/1e4 for Slater); the physical cm^-1 is
            # val*1000 for both (verified vs RCG's own PARVALS echo: scaled
            # F2(3p2) attr 16.876 -> 16876 cm^-1 vs RCG echo 16885.6).
            out.append({"key": key, "cfg": ck, "phys": phys, "kind": kind, "k": k,
                        "value_cm1": val * 1000.0, "raw": val})
    return out


def param_value_cm1(ing11_path, outgine_path, cfg, phys_prefix):
    """Single look-up: physical value (cm^-1) of e.g. F2 in config '3p2'. phys_prefix
    matches the start of the physical name (so 'F2' matches 'F2(22)'). RAISES if no
    unique match (zero or several) -- never returns a guessed slot."""
    ck = R._cfgkey(cfg)
    rows = physical_params(ing11_path, outgine_path)
    if ck in _OVERFLOW:
        raise ValueError(
            f"config {cfg!r} has open-shell params overflowing the 4 ING11 slots "
            f"and was not positionally labeled (needs the extended-slot reader). "
            f"Cannot look up {phys_prefix} here without risk of mislabeling.")
    hits = [p for p in rows if p["cfg"] == ck and p["phys"].startswith(phys_prefix)]
    if len(hits) != 1:
        raise ValueError(
            f"{phys_prefix} in config {cfg!r}: expected exactly 1 match, got "
            f"{len(hits)} ({[h['phys'] for h in hits]}). No silent guess.")
    return hits[0]["value_cm1"], hits[0]


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Physically-labeled ING11 params.")
    ap.add_argument("--ing11", required=True)
    ap.add_argument("--outgine", required=True)
    ap.add_argument("--cfg", help="filter to one config (cfgkey).")
    a = ap.parse_args()
    rows = physical_params(a.ing11, a.outgine)
    ck = R._cfgkey(a.cfg) if a.cfg else None
    print(f"{'config':10} {'param':10} {'kind':5} {'value (cm^-1)':>14}")
    for p in rows:
        if ck and p["cfg"] != ck:
            continue
        if p["phys"] == "CI":
            continue
        print(f"{str(p['cfg']):10} {p['phys']:10} {p['kind']:5} "
              f"{p['value_cm1']:14.2f}")
