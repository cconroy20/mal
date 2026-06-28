#!/usr/bin/env python3
"""
Substitute RCE's fitted parameters into an ING11 deck so a final RCG run yields
FITTED oscillator strengths (closing the semi-empirical loop).

RCE writes, at the end of PARVALS, a "PARAMETER VALUES FOR RCG INPUT" section
with the converged parameters per configuration and per CI pair, e.g.

    Mg I   3s3p            25.6396   0.0407   21.9152
    Mg I   3s4p            47.9659   0.0066    1.1021
     3s3p   - 3s4p         0.0000    8.4132

ING11 stores the same numbers as fixed-width 10-column integer fields equal to
(value * 1e5):

    Mg I   3s3p        3   2122740      2812   2788514        00        00hf85998585
    [0:10][10:20]=label+count  [20:30][30:40][40:50]=values   ...  [70:80]=scale code

and the CI (R^k) lines as:

     3s3p    - 3s4p    2   0.00005   1.93535 ...   (decimal, A-format)
(the single-config EAV/zeta/F/G lines are the integer-packed ones we rewrite).

This tool rewrites the single-configuration parameter lines (EAV, zeta, Slater
F^k/G^k) from the PARVALS "FOR RCG INPUT" block, matched by configuration name.
CI pair lines are left as-is (RCE kept those fixed in our setup).

Validation: substituting the *ab initio* PARVALS values back must reproduce the
original ING11 values bit-for-bit (round-trip), and re-running RCG must give the
same gf — see tests in work/run_mg1.sh.

Usage:
    python3 subst_fitted_params.py --parvals PARVALS --ing11 ING11 --out ING11.fit
"""
import argparse
import re

SCALE = 1e5


def read_parvals_rcg(path):
    """Return dict configname -> list of fitted parameter values (floats) from
    the 'PARAMETER VALUES FOR RCG INPUT' section. Only the single-configuration
    lines (label like 'Mg I   3s3p') are captured; CI pair lines (starting with
    a space + 'conf - conf') are skipped."""
    out = {}
    in_sec = False
    with open(path) as f:
        for ln in f:
            if "PARAMETER VALUES FOR RCG INPUT" in ln:
                in_sec = True
                continue
            # any other 'PARAMETER VALUES' header (cycle blocks) ends the section
            if "PARAMETER VALUES" in ln:
                in_sec = False
                continue
            if not in_sec:
                continue
            s = ln.rstrip("\n")
            if not s.strip():
                continue
            # CI pair lines look like '  3s3p   - 3s4p   ...' (have a ' - ')
            if " - " in s:
                continue
            # single-config line must START with a config label (letter), e.g.
            # 'Mg I   3s3p   <v1> <v2> ...'; numeric-only lines are skipped.
            if not re.match(r"\s*[A-Za-z]", s):
                continue
            m = re.match(r"(.{20})\s+(.*)$", s)
            if not m:
                continue
            label = m.group(1).strip()
            vals = [float(x) for x in re.findall(r"-?\d+\.\d+", m.group(2))]
            if vals:
                out[label] = vals
    return out


def _cfg_of_ing11_line(line):
    """Config label of an ING11 single-config parameter line, or None.
    Such lines have the param count at col 19 and 'hf' scaling code near col 70."""
    if len(line) < 70 or "hf" not in line[68:80]:
        return None
    # must have an integer param-count in cols 18:20
    if not re.match(r"\d", line[19:20]):
        return None
    label = line[0:18].strip()        # cols 0-17 are the config label; 18-19 = count
    # exclude the CI pair lines (they contain ' - ' in the label region)
    if " - " in line[0:18]:
        return None
    return label


def _fmt_eav(v):
    """EAV is read by RCG as F10.5 from a 10-char field written WITHOUT a decimal
    point, i.e. the integer round(v*1e5) right-justified in 10 cols."""
    return f"{round(v * 1e5):10d}"


def _fmt_param(v, group):
    """Non-EAV params (zeta, F^k, G^k) are F9.4 (9-char digit field, value*1e4)
    followed by a 1-digit I1 GROUP CODE. Preserve the original group code."""
    return f"{round(v * 1e4):9d}" + (group if group else "0")


def subst(parvals_path, ing11_path, out_path):
    """Rewrite the single-config parameter lines per RCG's actual format
    (2A6,A7,1X, F10.5 for EAV, then 4*(F9.4 value + I1 group code), A2, 4I2).
    Group codes and all trailing fields are preserved from the original line."""
    fitted = read_parvals_rcg(parvals_path)
    n_sub = 0
    with open(ing11_path) as f:
        lines = f.readlines()
    out = []
    for ln in lines:
        s = ln.rstrip("\n")
        cfg = _cfg_of_ing11_line(s)
        if cfg and cfg in fitted:
            vals = fitted[cfg]
            # EAV: cols [20:30]; params: 4*(F9.4 + I1) starting col 30
            new = s[0:20] + _fmt_eav(vals[0])
            o = 30
            for i in range(4):
                grp = s[o + 9:o + 10]                  # keep original group code
                v = vals[i + 1] if (i + 1) < len(vals) else None
                if v is None:
                    new += s[o:o + 10]                 # leave field untouched
                else:
                    new += _fmt_param(v, grp)
                o += 10
            new += s[o:]                               # A2 + trailing, unchanged
            out.append(new + "\n")
            n_sub += 1
        else:
            out.append(ln)
    with open(out_path, "w") as f:
        f.writelines(out)
    return n_sub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parvals", required=True)
    ap.add_argument("--ing11", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    n = subst(a.parvals, a.ing11, a.out)
    print(f"substituted {n} config parameter lines -> {a.out}")


if __name__ == "__main__":
    main()
