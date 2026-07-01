#!/usr/bin/env python3
"""Read and write the adjustable radial parameters in an ING11 deck.

ING11 stores two kinds of parameter lines:

  * single-configuration lines (EAV, ZETA, Slater F^k/G^k) -- integer-packed:
        Mg I   3s3p        3   2563750      4112   2015474        00        00hf...
        [0:18]=label [18:20]=count  then 4*(F9.4 value*1e4 + I1 group code) from
        col 30, with EAV as F10.5 (value*1e5) in [20:30].
  * CI pair lines (R^k interaction integrals) -- decimal A-format:
         3s2     - 3p2     1  29.40475   0.00005 ...
        value is F8.4 with a trailing 1-digit group code glued on (29.4047 + '5'),
        in one of 5 slots; the active slot is set by k.

A parameter is ADJUSTABLE iff its group code is nonzero. parse() returns the
ordered list of adjustable parameters (each with a stable key, current value,
and the location needed to write it back); write() puts a value vector back,
preserving every fixed field and all formatting. A read->write round-trip with
the same values reproduces the file byte-for-byte (see _selftest).
"""
import re

SCALE_EAV = 1e5
SCALE_P = 1e4
SCALE_CI = 1e4


# The single-config / CI marker in cols [68:80] is 'hf' for non-relativistic RCN
# (irel=0) and 'hr' for relativistic (irel=1) -- accept BOTH. (Without 'hr' every
# irel=1 ING11 line was silently rejected, so _raw_slots/param_labels/scale_hf all
# returned 0 params on a relativistic deck.)
def _hf_marker(line):
    seg = line[68:80]
    return "hf" in seg or "hr" in seg


def _is_singleconf(line):
    if len(line) < 70 or not _hf_marker(line):
        return False
    if not re.match(r"\d", line[19:20]):
        return False
    if " - " in line[0:18]:
        return False
    return True


def _is_ci(line):
    return " - " in line[0:18] and _hf_marker(line)


def values_by_key(path):
    """Map {key -> value} for EVERY single-config and CI field in `path`,
    regardless of group code (unlike parse(), which returns only adjustable
    params). Used to look up a parameter's ab-initio value as a ridge-prior
    centre even when that field is held fixed (group 0) in the ab-initio deck."""
    with open(path) as f:
        raw = f.readlines()
    out = {}
    first = True
    for ln in raw:
        s = ln.rstrip("\n")
        if _is_singleconf(s):
            cfg = s[0:18].strip()
            if not first:
                try:
                    out[f"{cfg}|EAV"] = int(s[20:30]) / SCALE_EAV
                except ValueError:
                    pass
            first = False
            o = 30
            for slot in range(4):
                try:
                    out[f"{cfg}|P{slot}"] = int(s[o:o + 9]) / SCALE_P
                except ValueError:
                    pass
                o += 10
        elif _is_ci(s):
            pair = s[0:18].strip()
            o = 20
            for slot in range(5):
                m = re.match(r"\s*(-?\d+\.\d+)", s[o:o + 10][:-1])
                if m:
                    out[f"{pair}|CI{slot}"] = float(m.group(1))
                o += 10
    return out


def parse(path):
    """Return (lines, params). `lines` is the raw file (list of str, no newline
    stripping issues). `params` is an ordered list of dicts:
        {key, kind, value, lineno, col, width, group}
    kind in {'EAV','P','CI'}; only ADJUSTABLE params (group code != 0) are
    included for EAV/P. CI params are included when their active slot is set.
    EAV of the FIRST single-config line (the ground config, pinned to 0) is
    excluded -- it is the energy zero, never adjusted."""
    with open(path) as f:
        raw = f.readlines()
    params = []
    first_singleconf = True
    for i, ln in enumerate(raw):
        s = ln.rstrip("\n")
        if _is_singleconf(s):
            cfg = s[0:18].strip()
            # EAV in [20:30], F10.5 packed (value*1e5). Adjustable unless this is
            # the ground config (pinned energy zero).
            if not first_singleconf:
                params.append({"key": f"{cfg}|EAV", "kind": "EAV",
                               "value": int(s[20:30]) / SCALE_EAV,
                               "lineno": i, "col": 20, "width": 10, "group": ""})
            first_singleconf = False
            # 4 * (F9.4 value + I1 group) from col 30
            o = 30
            for slot in range(4):
                raw_v = s[o:o + 9]
                grp = s[o + 9:o + 10]
                try:
                    v = int(raw_v) / SCALE_P
                except ValueError:
                    v = None
                if v is not None and grp.strip() and grp != "0":
                    params.append({"key": f"{cfg}|P{slot}", "kind": "P",
                                   "value": v, "lineno": i, "col": o,
                                   "width": 9, "group": grp})
                o += 10
        elif _is_ci(s):
            pair = s[0:18].strip()
            # decimal slots after the k index; find the nonzero (active) one.
            # slots start at col 20, each 10 wide: F8.4 value + 1-digit group.
            o = 20
            for slot in range(5):
                field = s[o:o + 10]
                # 10-char field = F8.4-ish value + 1-digit group code (last char).
                # The displayed number includes the group digit glued on, e.g.
                # '   1.83485' = value 1.8348, group '5'. Split off the last char.
                grp = field[-1]
                m = re.match(r"\s*(-?\d+\.\d+)", field[:-1])
                if m and grp.strip():
                    v = float(m.group(1))
                    # active CI value: not the 0.0000 placeholder
                    if abs(v) > 1e-3:
                        params.append({"key": f"{pair}|CI{slot}", "kind": "CI",
                                       "value": v, "lineno": i, "col": o,
                                       "width": 10, "group": grp})
                o += 10
    return raw, params


def _fmt_eav(v):
    return f"{round(v * SCALE_EAV):10d}"


def _fmt_p(v, grp):
    return f"{round(v * SCALE_P):9d}" + grp


def _fmt_ci(v, grp):
    # value (F.4) with the 1-digit group code glued on, right-justified to 10,
    # matching ING11's '   1.83485' = 1.8348 + group '5'.
    return f"{v:.4f}{grp}".rjust(10)


def write(raw_lines, params, values, out_path):
    """Write `values` (parallel to `params`) into a copy of raw_lines, preserving
    every fixed field, and save to out_path."""
    out = list(raw_lines)
    # group writes per line so multiple fields on one line compose correctly
    by_line = {}
    for p, v in zip(params, values):
        by_line.setdefault(p["lineno"], []).append((p, v))
    for lineno, edits in by_line.items():
        s = out[lineno].rstrip("\n")
        for p, v in edits:
            c = p["col"]
            if p["kind"] == "EAV":
                field = _fmt_eav(v)        # 10 chars, at [20:30]
            elif p["kind"] == "P":
                field = _fmt_p(v, p["group"])   # 9-digit value + 1 group = 10
            else:  # CI
                field = _fmt_ci(v, p["group"])  # F8.4 + group, 10 chars
            assert len(field) == 10, (p, repr(field))
            s = s[:c] + field + s[c + 10:]
        out[lineno] = s + "\n"
    with open(out_path, "w") as f:
        f.writelines(out)


def _selftest(path):
    """Read params, write them back unchanged -> file must be byte-identical."""
    import tempfile, os, filecmp
    raw, params = parse(path)
    vals = [p["value"] for p in params]
    fd, tmp = tempfile.mkstemp()
    os.close(fd)
    write(raw, params, vals, tmp)
    ok = filecmp.cmp(path, tmp, shallow=False)
    if not ok:
        # show first differing line
        a = open(path).readlines()
        b = open(tmp).readlines()
        for i, (x, y) in enumerate(zip(a, b)):
            if x != y:
                print(f"DIFF line {i}:\n  orig: {x!r}\n  new:  {y!r}")
                break
    os.remove(tmp)
    print(f"round-trip {'OK' if ok else 'FAILED'}  ({len(params)} params)")
    for p in params:
        print(f"  {p['key']:18} kind={p['kind']:3} val={p['value']:12.5f}")
    return ok


if __name__ == "__main__":
    import sys
    _selftest(sys.argv[1])
