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
import os
import re

KK = 1000.0


def _termkey(term):
    """multiplicity+L core, ignoring any parent: '(2S) 3P'->'3P', '3P'->'3P'."""
    t = re.sub(r"\([^)]*\)", " ", term)
    ms = re.findall(r"\d[A-Z]", t)
    return ms[-1] if ms else term.strip()


# n (one or more digits) + l + optional occupation; the occupation digit only
# counts when NOT followed by another orbital letter or digit, so '3s11d' parses
# as 3s + 11d (not 3s^1 + 1d) and '3p2' as 3p^2.
_ORB = re.compile(r"\d+[spdfghik](?:\d(?![spdfghik\d]))?")


def _cfgkey(config):
    """Reduce a config to its outer valence orbital(s), occupation-normalized,
    matching the scheme in make_report so NIST and computed labels agree."""
    toks = _ORB.findall(config)
    norm = []
    for t in toks:
        m = re.match(r"(\d+[spdfghik])(\d?)", t)
        nl, occ = m.group(1), m.group(2)
        norm.append(nl + (occ if occ and occ != "1" else ""))
    if norm and re.fullmatch(r"\d+[spdfghik]2", norm[-1]):
        return norm[-1]
    return ".".join(norm[-2:]) if norm else config.strip()


def load_nist_levels(path, max_E=None):
    """Return dict (parity, Jstr) -> list of dicts {E_kK, cfgk, tk} for observed
    (non-predicted) levels, keeping config and term so we can match by identity.
    `max_E` (cm^-1) drops levels above an energy ceiling -- use the ionization
    limit to fit only the BOUND spectrum (Bob fits a curated low set; our NIST
    cache also carries hundreds of high-Rydberg/autoionizing levels that, fit
    over an incomplete basis, make RCE thrash)."""
    table = {}
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            config, term, J, level, parity, pred = p[0], p[1], p[2], p[3], p[4], p[5]
            try:
                Jf = float(J); E = float(level)
            except ValueError:
                continue
            if pred == "1" or (max_E is not None and E > max_E):
                continue
            table.setdefault((parity, "%g" % Jf), []).append(
                {"E": E / KK, "cfgk": _cfgkey(config), "tk": _termkey(term)})
    return table


def load_computed_terms(outg11_path):
    """From OUTG11, get computed levels per (parity, Jstr) in ENERGY ORDER, each
    tagged with its ROBUST identity (cfgk, tk) = the DOMINANT eigenvector basis
    state (config, term). This is reliable even for strongly-mixed / reordered
    levels, unlike the per-level config in the ENERGY-MATRIX header (which labels
    every level by the block's first config). Returns
    dict (parity,J) -> [{E, cfgk, tk}, ...] energy-ordered."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "parse_cowan", os.path.join(here, "parse_cowan.py"))
    pc = importlib.util.module_from_spec(spec); spec.loader.exec_module(pc)
    comp = pc.parse_compositions(outg11_path)
    out = {}
    for key, levs in comp.items():
        out[key] = [{"E": d["E_calc"] / KK, "cfgk": _cfgkey(d["config"]),
                     "tk": _termkey(d["term"])} for d in levs]
    return out


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
    """A parameter-group-code line: a PURE-INTEGER line (no letters/parens/dots)
    whose every entry is a TIE-GROUP CODE -- 0 or a multiple of 100 (+100 = an
    adjustable parameter, -100 = fixed; larger multiples tie params together).
    This is the exact, general RCG/OUTGINE convention (verified on the 9- and
    122-config Mg I decks), and it cleanly excludes the three lookalikes that
    broke earlier heuristics: float value/eigenvalue lines (have '.'),
    parameter-NAME lines containing CI names like '120D1114' (have letters),
    and sequential per-level index lines '99 100 101 102' (not all mult-of-100)."""
    s = line.rstrip("\n")
    if not s.strip() or re.search(r"[A-Za-z().]", s):
        return False
    ints = re.findall(r"-?\d+", s)
    if not ints:
        return False
    return all(int(x) % 100 == 0 for x in ints)


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


def build(outgine_path, nist_levels, computed_terms, out_path, free_ci_with=(),
          freeze=(), ruleset=None):
    with open(outgine_path) as f:
        lines = f.readlines()
    return _build_focused(lines, nist_levels, computed_terms, out_path,
                          free_ci_with, freeze, ruleset)


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


def _block_orbitals(hdr):
    """Ordered l-letters of the block's orbitals, from the config-names header
    (used to map F^k/G^k orbital indices (ij) to l). Best-effort."""
    text = " ".join(hdr)
    # the first config-names line lists configs; take distinct orbitals in order
    orbs = []
    for nl in re.findall(r"\d([spdfghik])", text):
        if nl not in orbs:
            orbs.append(nl)
    return orbs


def _build_focused(lines, nist_levels, computed_terms, out_path, free_ci_with=(),
                   freeze=(), ruleset=None):
    """Walk blocks; within each, replace consecutive (value-line, flag-line)
    pairs (one per J, ascending from the block's minimum J) with observed
    levels. Each computed J-slot (energy-ordered) is matched to the observed
    NIST level of the SAME (config, term) identity -- not by energy rank, which
    swaps near-degenerate same-J terms (e.g. 3s3d 1D vs 3D)."""
    out = list(lines)
    n = len(lines)
    i = 0
    block_index = 0
    freed_tally = __import__("collections").Counter()

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
                Jkey = (parity, "%g" % Jval)
                slots = computed_terms.get(Jkey, [])     # energy-ordered
                obs = list(nist_levels.get(Jkey, []))
                # Match by exact (config, term) IDENTITY -- both the computed
                # slot (dominant eigenvector basis state) and the NIST level carry
                # a reliable config+term now. Within an identity, pair in energy
                # order. Robust to mixing and energy reordering.
                obs_by_id = {}
                for o in sorted(obs, key=lambda d: d["E"]):
                    obs_by_id.setdefault((o["cfgk"], o["tk"]), []).append(o)
                seen_id = {}
                newvals, newflags = [], []
                for m in range(ncomp):
                    comp = float(FLOAT7.findall(lines[j])[m])
                    match = None
                    if m < len(slots):
                        ident = (slots[m]["cfgk"], slots[m]["tk"])
                        rank = seen_id.get(ident, 0)
                        seen_id[ident] = rank + 1
                        cand = obs_by_id.get(ident, [])
                        if rank < len(cand):
                            match = cand[rank]
                    if match is not None:
                        newvals.append(match["E"]); newflags.append(m + 1)
                    else:
                        newvals.append(comp); newflags.append(-(m + 1))  # exclude
                out[j] = fmt_values(newvals) + "\n"
                out[j + 1] = fmt_flags(newflags) + "\n"
                j += 2
                Jval += 1.0
            # --- free the physical parameters: rewrite the group-code block ---
            # The group-code run (lines of 0/+-100..) is the next CONTIGUOUS run
            # of group-code lines; on large bases it is NOT immediately after the
            # J-blocks (more J value/index lines intervene when a J has many
            # levels), so scan forward to it rather than assuming adjacency. We
            # must not cross into the NEXT parity block, so stop at its header.
            pnames, pcfgs = _parse_param_names(hdr)
            scan = j
            while scan < n and not _is_groupcode_line(lines[scan]) \
                    and not is_block_header(scan):
                scan += 1
            gc_lines = []
            while scan < n and _is_groupcode_line(lines[scan]):
                gc_lines.append(scan); scan += 1
            j = scan
            if gc_lines and pnames:
                codes = []
                for gj in gc_lines:
                    codes.extend(int(x) for x in re.findall(r"-?\d+", lines[gj]))
                # INVARIANT: the group-code run must align 1:1 with the parameter-
                # name list (name[L] <-> code[L]). If these ever diverge, the
                # free/freeze decisions would be applied to the wrong parameters
                # -- a silent, dangerous failure. Refuse rather than mis-fit.
                if len(codes) != len(pnames):
                    raise ValueError(
                        f"OUTGINE block {block_index}: {len(pnames)} parameter "
                        f"names but {len(codes)} group codes -- name<->code "
                        f"misalignment; the parser would free the wrong params. "
                        f"(Check _parse_param_names / _is_groupcode_line.)")
                borbs = _block_orbitals(hdr) if ruleset else None
                newcodes = []
                for L, code in enumerate(codes):
                    name = pnames[L] if L < len(pnames) else ""
                    cfg = pcfgs[L] if L < len(pcfgs) else None
                    is_ref = _is_eav_reference(name, L, block_index == 0)
                    if ruleset is not None:
                        # ruleset-driven (notes/bob_fit_ruleset.md): free by family
                        # + shell, never the ground-EAV reference, honor --freeze.
                        free = (not is_ref
                                and not _is_frozen(name, cfg, freeze)
                                and _ruleset_free(name, cfg,
                                                  ruleset["observed_cfgs"],
                                                  ruleset["open_ls"],
                                                  ruleset["is_ion"], borbs))
                        newcodes.append(0 if free else code)
                        if free:
                            freed_tally[_param_family(name)] += 1
                    elif (_is_physical_param(name) and not is_ref
                            and not _is_frozen(name, cfg, freeze)):
                        newcodes.append(0)        # free physical param
                    elif _is_free_ci(name, free_ci_with):
                        newcodes.append(0)        # free selected CI integral
                    else:
                        newcodes.append(code)     # keep (fixed)
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
    return dict(freed_tally)


def _parse_param_names(hdr_lines):
    """Parameter names are A6,A4 (10-char) fields packed 7 per line, appearing
    after the config-names line in the block header. Return (names, configs):
    two parallel ordered lists, where configs[i] is the configuration that owns
    parameter names[i]. Each 'EAV  <cfg>' field opens a new config; subsequent
    non-EAV params (ZETA, F^k, G^k) belong to that config until the next EAV."""
    names, configs = [], []
    cur_cfg = None
    started = False
    for ln in hdr_lines:
        if re.search(r"\bEAV\b|ZETA|^\s*[FG]\d", ln) or started:
            if re.search(r"EAV|ZETA|[FG]\d|\dD\d|\dE\d", ln):
                started = True
                # split into 10-char fields
                s = ln.rstrip("\n")
                for c in range(0, len(s), 10):
                    field = s[c:c + 10].strip()
                    if not field:
                        continue
                    m = re.match(r"EAV\s+(\S+)", field)
                    if m:
                        cur_cfg = m.group(1)
                    names.append(field)
                    configs.append(cur_cfg)
            elif started:
                break
    return names, configs


def _is_physical_param(name):
    """Single-configuration physical parameter: EAV, ZETA, Slater F^k/G^k."""
    return bool(re.match(r"(EAV|ZETA|[FG]\d)", name))


# -------------------------- ruleset-driven free/freeze --------------------------
# Implements notes/bob_fit_ruleset.md (derived from Bob's fits of 10 species):
# free only the integrals the observed level structure directly constrains, by
# parameter family, shell-aware. Everything else stays FIXED (group code kept) and
# is held at scaled HF by the prior. Off-diagonal CI (R^k) is NEVER freed here.
_FG = re.compile(r"([FG])(\d)\((\d)(\d)\)")


def _param_family(name):
    """Coarse family of a single-config parameter name, for free-tally reporting."""
    nm = name.strip()
    if nm.startswith("EAV"):
        return "EAV"
    if nm.startswith("ZETA"):
        return "ZETA"
    if nm.startswith("ALPHA") or nm.startswith("BETA"):
        return "ALPHA/BETA"
    if re.match(r"F\d", nm):
        return "F^k"
    if re.match(r"G\d", nm):
        return "G^k"
    return "other"


def _orbital_l(config, idx, block_orbitals):
    """l-letter of the idx-th (1-based) orbital of `config`, using the block's
    orbital list when available; falls back to the config's own orbitals."""
    if block_orbitals and 1 <= idx <= len(block_orbitals):
        return block_orbitals[idx - 1]
    orbs = re.findall(r"\d([spdfghik])", config or "")
    return orbs[idx - 1] if 1 <= idx <= len(orbs) else None


def _ruleset_free(name, cfg, observed_cfgs, open_ls, is_ion, block_orbitals):
    """Decide free (True) vs freeze (False) for a single-config physical param,
    per the ruleset. `observed_cfgs` = set of _cfgkey strings with >=1 observed
    level; `open_ls` = set of open-shell l-letters in the ground config (e.g.
    {'p'} for p^2, {'d'} for d^n); `is_ion` = charge>0."""
    obs = _cfgkey(cfg) in observed_cfgs if cfg else False
    nm = name.strip()
    # Rule 1: EAV of any observed config (the global ground EAV is handled by the
    # caller's _is_eav_reference pin).
    if nm.startswith("EAV"):
        return obs
    # Rule 3: exchange G^k -- free where an observed multiplet splitting exists.
    if re.match(r"G\d", nm):
        return obs
    # Rule 4: spin-orbit ZETA -- free where fine structure is observed (the
    # ion-dependent precision floor is a refinement handled by the prior width).
    if nm.startswith("ZETA"):
        return obs
    # Rule 5: direct F^k -- only for the open/same shell (intra-shell pair).
    if re.match(r"F\d", nm):
        m = _FG.match(nm)
        if not m or not obs:
            return False
        i, j = int(m.group(3)), int(m.group(4))
        if i == j:                                   # same-orbital (e.g. 3p^2)
            li = _orbital_l(cfg, i, block_orbitals)
            return li in open_ls if open_ls else True
        return False                                 # cross-shell direct: freeze
    # Rule 6: Trees ALPHA (open shell, always) / BETA (ions only, sparingly).
    if nm.startswith("ALPHA"):
        return bool(open_ls) and obs
    if nm.startswith("BETA"):
        return bool(open_ls) and obs and is_ion
    return False


def _ci_configs(name):
    """For a CI parameter name like '161D1122' or '120D1113', return the pair of
    coupled configuration indices. The name is <c1><c2><k><D|E><....>: first two
    digits are the 1-based config numbers, third is the multipole k. So
    '161D1122' -> configs (1,6), '120D1113' -> (1,2). None if not a CI name."""
    m = re.match(r"(\d)(\d)\d[DE]\d+$", name)
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2)))


def _is_free_ci(name, free_ci_pairs):
    """Free a configuration-interaction integral only if it couples one of the
    specific config PAIRS in `free_ci_pairs` (a set of frozensets of 1-based
    config indices). Targeting individual pairs (e.g. {1,6} = 3s^2-3p^2) avoids
    the instability of freeing many correlated CI integrals at once."""
    pair = _ci_configs(name)
    if pair is None or not free_ci_pairs:
        return False
    return frozenset(pair) in free_ci_pairs


def _is_frozen(name, cfg, freeze):
    """Hold a single-configuration physical parameter at its ab-initio value
    instead of letting the energy fit move it. `freeze` is a set of (config,
    param-prefix) tuples, e.g. ('3s4d','G2') or ('3s5p','G1'). Use this for
    UNDER-CONSTRAINED Rydberg exchange/Slater integrals that the levels fit can
    drive to unphysical values (e.g. G->0), wrecking the gf even as energies
    improve. Matching is by config substring + param name prefix; an empty
    config in the tuple freezes that param prefix in EVERY config."""
    if not freeze or not cfg:
        return False
    nm = name.strip()
    for fc, fp in freeze:
        if nm.startswith(fp) and (not fc or fc == cfg):
            return True
    return False


def _is_eav_reference(name, L, is_first_block):
    """Hold ONLY the global ground-state EAV (first EAV of the first/even block)
    fixed as the energy zero. Every other block's EAV must be free, or its
    levels cannot reach the observed energies."""
    return bool(is_first_block and name.startswith("EAV") and L == 0)


def _ground_open_shells(nist_path):
    """Open-shell l-letters of the lowest observed config (partially-filled l).
    e.g. Mg I 3s^2 -> set() (closed); Fe II 3d^7 -> {'d'}; Si I 3s^2 3p^2 -> {'p'}.
    Reads the raw config string (col 0) of the lowest-energy observed level."""
    full = {"s": 2, "p": 6, "d": 10, "f": 14}
    best_E, best_cfg = None, None
    with open(nist_path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 6 or p[5] == "1":
                continue
            try:
                E = float(p[3])
            except ValueError:
                continue
            if best_E is None or E < best_E:
                best_E, best_cfg = E, p[0]
    opens = set()
    for n, l, occ in re.findall(r"(\d)([spdfghik])(\d*)", best_cfg or ""):
        o = int(occ) if occ else 1
        if l in full and 0 < o < full[l]:
            opens.add(l)
    return opens, best_cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outgine", required=True)
    ap.add_argument("--nist", required=True)
    ap.add_argument("--outg11", required=True,
                    help="ab initio OUTG11, for per-slot term identities.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--free-ci-pairs", default="",
                    help="comma-separated config PAIRS (1-based, 'i-j') whose CI "
                         "integral should be freed, e.g. '1-6' for 3s^2-3p^2. "
                         "Target specific pairs to avoid over-freeing.")
    ap.add_argument("--freeze-params", default="",
                    help="comma-separated 'config:PARAM' to HOLD at ab-initio "
                         "value (not freed by the energy fit), e.g. "
                         "'3s4d:G2,3s5p:G1'. Use for under-constrained Rydberg "
                         "exchange integrals the fit drives unphysical, which "
                         "wrecks gf. Empty config ':G1' freezes that param in "
                         "all configs.")
    ap.add_argument("--ruleset", action="store_true",
                    help="select free params by the physics ruleset "
                         "(notes/bob_fit_ruleset.md): EAV(observed) + G^k + ZETA "
                         "where observed + same-shell F^k + Trees ALPHA/BETA for "
                         "open shells; never free off-diagonal CI. Shell-aware.")
    ap.add_argument("--ion", action="store_true",
                    help="this species is an ION (charge>0); enables Trees BETA "
                         "freeing per the ruleset. Default: neutral.")
    ap.add_argument("--max-energy", type=float, default=None,
                    help="fit only observed levels below this energy (cm^-1), "
                         "e.g. the ionization limit, to exclude high-Rydberg / "
                         "autoionizing levels that make the fit thrash.")
    a = ap.parse_args()
    free_ci = set()
    for tok in a.free_ci_pairs.split(","):
        tok = tok.strip()
        if "-" in tok:
            i, j = tok.split("-"); free_ci.add(frozenset((int(i), int(j))))
    freeze = set()
    for tok in a.freeze_params.split(","):
        tok = tok.strip()
        if ":" in tok:
            cfg, prm = tok.split(":", 1)
            freeze.add((cfg.strip(), prm.strip()))
    nist = load_nist_levels(a.nist, max_E=a.max_energy)
    cterms = load_computed_terms(a.outg11)
    if a.max_energy is not None:
        print(f"  fitting {sum(len(v) for v in nist.values())} observed levels "
              f"below {a.max_energy:.0f} cm^-1")
    ruleset = None
    if a.ruleset:
        observed = {d["cfgk"] for v in nist.values() for d in v}
        open_ls, gcfg = _ground_open_shells(a.nist)
        ruleset = {"observed_cfgs": observed, "open_ls": open_ls,
                   "is_ion": a.ion}
        print(f"  ruleset: ground={gcfg!r} open_shells={open_ls or 'closed'} "
              f"ion={a.ion}  observed_configs={len(observed)}")
    tally = build(a.outgine, nist, cterms, a.out, free_ci, freeze, ruleset)
    print(f"wrote {a.out}  (free CI with configs {free_ci or 'none'}; "
          f"frozen params {freeze or 'none'}; "
          f"ruleset={'on' if ruleset else 'off'})")
    if tally:
        print(f"  freed by family: {tally}  (total {sum(tally.values())})")


if __name__ == "__main__":
    main()
