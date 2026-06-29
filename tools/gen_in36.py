#!/usr/bin/env python3
"""Generate an RCN in36 configuration deck for Mg I from a list of config labels
(matching Bob Kurucz's b1200{e,o}.com config set). Reproduces the exact column
format of the hand-written work/mg1/in36 (verified byte-for-byte on the original
9 configs).

A config label is Bob's compact form: a sequence of <n><l>[<occ>] orbital
tokens, e.g. '3s2' (3s^2), '3d4s' (3d 4s), '3p2' (3p^2), '3s11d' (3s 11d),
'3d9i' (3d 9i). For Mg I every config sits on the closed 2p^6 core, so the RCN
occupation string is '2p6' + the valence orbitals.

in36 config line layout (cols, 0-based):
  [0:10]  '   12    0'   Z=12, charge 0
  [10:15] 'Mg I '        spectrum name
  [15:30] '  ' + label   (label left-justified in 13 cols)
  [30:]   '    ' + each orbital token left-justified in a 5-char field
          (core 2p6 first, then valence)

Usage:
  tools/gen_in36.py --out work/mg1_full/in36          # default = Bob's full set
  tools/gen_in36.py --configs 3s2,3s3p,3p2 --out ...
"""
import argparse
import re

HEADER = ("200-90 0 2  01.  4.0    5.E-08    1.E-11"
          "-2 00090 0 1.0  0.65  0.0 1.00   -6")
CORE = "2p6"

# Bob Kurucz's Mg I configuration set (b1200e.com + b1200o.com), even then odd.
BOB_EVEN = [
    "3s2", "3s4s", "3s5s", "3s6s", "3s7s", "3s8s", "3s9s", "3s10s", "3s11s",
    "3d4s", "3d5s", "3d6s", "3d7s", "3d8s", "3d9s", "3d10s", "3d11s",
    "3s3d", "3s4d", "3s5d", "3s6d", "3s7d", "3s8d", "3s9d", "3s10d", "3s11d",
    "3d2", "3d4d", "3d5d", "3d6d", "3d7d", "3d8d", "3d9d", "3d10d", "3d11d",
    "3s5g", "3s6g", "3s7g", "3s8g", "3s9g", "3d5g", "3d6g", "3d7g", "3d8g",
    "3d9g", "3s7i", "3s8i", "3s9i", "3d7i", "3d8i", "3d9i",
    "3p2", "3p4p", "3p5p", "3p6p", "3p7p", "3p8p", "3p9p", "3p10p", "3p11p",
    "3p12p",
]
BOB_ODD = [
    "3s3p", "3s4p", "3s5p", "3s6p", "3s7p", "3s8p", "3s9p", "3s10p", "3s11p",
    "3s12p", "3d4p", "3d5p", "3d6p", "3d7p", "3d8p", "3d9p", "3d10p", "3d11p",
    "3d12p", "3s4f", "3s5f", "3s6f", "3s7f", "3s8f", "3s9f", "3s10f", "3s11f",
    "3d4f", "3d5f", "3d6f", "3d7f", "3d8f", "3d9f", "3d10f", "3d11f",
    "3s6h", "3s7h", "3s8h", "3s9h", "3d6h", "3d7h", "3d8h", "3d9h",
    "3p4s", "3p5s", "3p6s", "3p7s", "3p8s", "3p9s", "3p10s", "3p11s", "3p12s",
    "3p3d", "3p4d", "3p5d", "3p6d", "3p7d", "3p8d", "3p9d", "3p10d", "3p11d",
]

# n, l, optional occupation. n is the maximal digit run BEFORE an orbital letter.
# An occupation digit only counts when it is the LAST char or is followed by
# another orbital letter -- NOT when followed by a digit (that digit belongs to
# the next orbital's n, e.g. '3s11d' = 3s + 11d, and '3s4s' = 3s + 4s).
_ORB = re.compile(r"(\d+)([spdfghiklm])(\d(?![spdfghiklm\d]))?")


def orbital_tokens(label):
    """Valence orbital tokens for a config label, e.g. '3d4s'->['3d','4s'],
    '3p2'->['3p2'], '3s11d'->['3s','11d']. Occupation digit kept only when >1."""
    toks = []
    for n, l, occ in _ORB.findall(label):
        toks.append(f"{n}{l}{occ}" if occ and occ != "1" else f"{n}{l}")
    if not toks:
        raise ValueError(f"unparseable config: {label!r}")
    return toks


def config_line(label):
    occ = "    " + "".join(f"{t:<5}" for t in [CORE] + orbital_tokens(label))
    return f"   12    0Mg I   {label:<13}"[:30] + occ


def build(configs):
    return "\n".join([HEADER] + [config_line(c) for c in configs] + ["   -1"]) \
        + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--configs", default=None,
                    help="comma-separated config labels; default = Bob's full set")
    a = ap.parse_args()
    configs = ([c.strip() for c in a.configs.split(",")] if a.configs
              else BOB_EVEN + BOB_ODD)
    import os
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        f.write(build(configs))
    print(f"wrote {a.out}  ({len(configs)} configs)")


if __name__ == "__main__":
    main()
