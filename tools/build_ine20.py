#!/usr/bin/env python3
"""
Build an RCE input deck (INE20) from RCG's OUTGINE template by substituting the
OBSERVED (NIST) energy levels in place of the theoretical term-value
placeholders, so RCE fits the radial parameters to experiment.

How OUTGINE is structured (per parity block):
  ... header / config names / parameter names ...
  for each J (ascending):
      <line of term values T(M) for that J>   (Format 7F10.5, kK = 1000 cm^-1)
      <line of flags NF(M)>                    (Format 7I10)
  ... group codes / denominators / parameter values / trailer / -1 ...

We replace each J-block's T-values with observed levels and set the flag NF(M)
= +M to include a level in the fit, or -M to exclude it (no NIST match). The
matching strategy: within each (parity-block, J), sort computed eigenvalues and
NIST levels of that J & parity ascending and pair them in order (RCE itself can
reorder, but order-matching the close-lying ab initio values is robust for the
low levels we target). Unmatched computed slots are excluded (flag -M); we never
invent observed data.

Energies are in 1000 cm^-1 (kK) throughout OUTGINE.

Usage:
    python3 build_ine20.py --outgine work/mg1/OUTGINE \
        --nist data/nist/MgI_levels.tsv --out work/mg1/INE20
"""
import argparse
import re

KK = 1000.0


def load_nist_by_Jparity(path):
    """Return dict (parity, Jstr) -> sorted list of observed E in kK."""
    table = {}
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            J, level, parity, pred = p[2], p[3], p[4], p[5]
            try:
                Jf = float(J)
                E = float(level)
            except ValueError:
                continue
            if pred == "1":           # skip NIST-flagged predicted/uncertain
                continue
            key = (parity, "%g" % Jf)
            table.setdefault(key, []).append(E / KK)
    for k in table:
        table[k].sort()
    return table


FLOAT7 = re.compile(r"(-?\d+\.\d+|-?\d+\.|\d*\.\d+)")


def is_value_line(line):
    """A T-value line: only floats (7F10.5), no letters."""
    s = line.rstrip("\n")
    if not s.strip():
        return False
    if re.search(r"[A-Za-z]", s):
        return False
    nums = FLOAT7.findall(s)
    # value lines have floats with decimals; flag lines are pure integers
    return len(nums) > 0 and "." in s


def is_flag_line(line):
    s = line.rstrip("\n")
    if not s.strip() or "." in s:
        return False
    return bool(re.fullmatch(r"[\s\d-]+", s)) and bool(re.search(r"\d", s))


def _is_groupcode_line(line):
    """The parameter-group-code lines are integer lines whose entries are large
    in magnitude (|code| >= 100). Per-level flag lines are small sequential
    integers (1,2,3,...). Use the magnitude to tell them apart."""
    ints = re.findall(r"-?\d+", line)
    if not ints:
        return False
    return any(abs(int(x)) >= 100 for x in ints)


def _is_level_flag_line(line):
    """A genuine per-level flag line: small integers, each |v| <= 50 (level
    indices), e.g. '   1   2   3' or with negatives for excluded levels."""
    if not is_flag_line(line):
        return False
    ints = [int(x) for x in re.findall(r"-?\d+", line)]
    return bool(ints) and all(abs(v) <= 50 for v in ints)


def fmt_values(vals):
    return "".join(f"{v:10.4f}" for v in vals)


def fmt_flags(flags):
    return "".join(f"{f:10d}" for f in flags)


def build(outgine_path, nist_table, out_path):
    with open(outgine_path) as f:
        lines = f.readlines()
    _build_focused(lines, nist_table, out_path)


def _infer_block_parity(header_block_lines):
    """Infer parity from the config-names line of a block."""
    text = " ".join(header_block_lines)
    cfgs = re.findall(r"(\d[spdfg]\d?)(?![spdfg])", text)
    lmap = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}
    if not cfgs:
        return "e"
    # parity of the first listed configuration (all configs in a block share it)
    # take its two valence orbitals
    lsum = 0
    for c in cfgs[:2]:
        mm = re.match(r"\d([spdfg])(\d?)", c)
        l = lmap[mm.group(1)]; occ = int(mm.group(2)) if mm.group(2) else 1
        lsum += l * occ
    return "o" if lsum % 2 else "e"


def _build_focused(lines, nist_table, out_path):
    """Walk blocks; within each, replace consecutive (value-line, flag-line)
    pairs (one per J, ascending from the block's minimum J) with observed
    levels for that (parity, J)."""
    out = list(lines)
    n = len(lines)
    i = 0
    block_index = 0

    def is_block_header(k):
        # A genuine block header (7 ints) is immediately followed by a line of
        # the form 'N 0 1 1  0.000  1  3  1.0000' (the config-set descriptor).
        if not re.match(r"\s*\d+\s+\d+\s+\d+\s+-?\d+\s+\d+\s+\d+\s+\d+\s*$",
                        lines[k]):
            return False
        if k + 1 >= n:
            return False
        return bool(re.match(r"\s*\d+\s+0\s+1\s+1\s+[\d.]+", lines[k + 1]))

    while i < n:
        # block boundary
        if is_block_header(i):
            # gather block header lines (until first value/flag pair) for parity
            j = i + 1
            hdr = []
            while j < n and not (is_value_line(lines[j]) and
                                 j + 1 < n and is_flag_line(lines[j + 1])):
                hdr.append(lines[j]); j += 1
            parity = _infer_block_parity(hdr)
            # now consume value/flag pairs == J-blocks, ascending J
            # integer J for even-electron neutral (Mg I): J = 0,1,2,...
            # STOP when we reach the parameter-group-code section: those are
            # flag-like lines whose entries are large (|code| >= 100), unlike the
            # small per-level flags (1,2,3,...). The T-value section ends there;
            # everything after (group codes, denominators, parameter values,
            # trailer, -1) must be left untouched.
            Jval = 0.0
            while (j + 1 < n and is_value_line(lines[j])
                   and _is_level_flag_line(lines[j + 1])):
                ncomp = len(FLOAT7.findall(lines[j]))
                obs = list(nist_table.get((parity, "%g" % Jval), []))
                newvals, newflags = [], []
                for m in range(ncomp):
                    if m < len(obs):
                        newvals.append(obs[m]); newflags.append(m + 1)
                    else:
                        # no observed level: keep computed placeholder, exclude
                        comp = float(FLOAT7.findall(lines[j])[m])
                        newvals.append(comp); newflags.append(-(m + 1))
                out[j] = fmt_values(newvals) + "\n"
                out[j + 1] = fmt_flags(newflags) + "\n"
                j += 2
                Jval += 1.0
            # --- free the physical parameters: rewrite the group-code block ---
            # The group codes immediately follow the J-blocks (lines of +-100..).
            # Parameter NAMES live in the header lines (after the config-names
            # line). We free single-configuration physical parameters (EAV,
            # ZETA, F^k, G^k) by setting their code 100 -> 0, and keep the
            # configuration-interaction parameters (names like '120D1212') fixed.
            pnames = _parse_param_names(hdr)
            gc_start = j
            gc_lines = []
            while j < n and _is_groupcode_line(lines[j]):
                gc_lines.append(j); j += 1
            if gc_lines and pnames:
                codes = []
                for gj in gc_lines:
                    codes.extend(int(x) for x in re.findall(r"-?\d+", lines[gj]))
                newcodes = []
                for L, code in enumerate(codes):
                    name = pnames[L] if L < len(pnames) else ""
                    if (_is_physical_param(name)
                            and not _is_eav_reference(name, L,
                                                      block_index == 0)):
                        newcodes.append(0)        # free
                    else:
                        newcodes.append(code)     # keep (fixed / CI)
                # re-emit group-code lines, 7 per line (Format 7I10)
                for k, gj in enumerate(gc_lines):
                    chunk = newcodes[k * 7:(k + 1) * 7]
                    out[gj] = "".join(f"{c:10d}" for c in chunk) + "\n"
            block_index += 1
            i = j
        else:
            i += 1
    with open(out_path, "w") as f:
        f.writelines(out)


def _parse_param_names(hdr_lines):
    """Parameter names are A6,A4 (10-char) fields packed 7 per line, appearing
    after the config-names line in the block header. Return ordered list."""
    # the config-names line has the species + config tokens; the parameter-name
    # lines are those containing tokens like 'EAV', 'ZETA', 'G1(12)', '120D1212'.
    names = []
    started = False
    for ln in hdr_lines:
        if re.search(r"\bEAV\b|ZETA|^\s*[FG]\d", ln) or started:
            if re.search(r"EAV|ZETA|[FG]\d|\dD\d|\dE\d", ln):
                started = True
                # split into 10-char fields
                s = ln.rstrip("\n")
                for c in range(0, len(s), 10):
                    field = s[c:c + 10].strip()
                    if field:
                        names.append(field)
            elif started:
                break
    return names


def _is_physical_param(name):
    """Single-configuration physical parameter: EAV, ZETA, Slater F^k/G^k."""
    return bool(re.match(r"(EAV|ZETA|[FG]\d)", name))


def _is_eav_reference(name, L, is_first_block):
    """Hold ONLY the global ground-state EAV (first EAV of the first/even block)
    fixed as the energy zero. Every other block's EAV must be free, or its
    levels cannot reach the observed energies."""
    return bool(is_first_block and name.startswith("EAV") and L == 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outgine", required=True)
    ap.add_argument("--nist", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    nist = load_nist_by_Jparity(a.nist)
    build(a.outgine, nist, a.out)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
