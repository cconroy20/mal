#!/usr/bin/env python3
"""
Parse Cowan RCG output (OUTG11) into tidy tables for diagnostics.

Two products:
  1. levels : computed energy levels (config, term, J, parity, E_calc in cm^-1),
              extracted from the per-J "EIGENVALUES" blocks.
  2. lines  : the electric-dipole transitions (lower/upper E, J, term, wavelength,
              log gf, GA), from the "ELEC DIP SPECTRUM" table.

This is the parsing layer for the per-species diagnostic report. When RCE is
wired up, an analogous reader will pull the *fitted* levels (and the observed
levels merged in) from LEVELS1/2/3; the report compares the two.

Energies in OUTG11 are in units of 1000 cm^-1; we convert to cm^-1.

Usage (as a library):
    from parse_cowan import parse_outg11
    levels, lines = parse_outg11("work/mg1/OUTG11")
"""
import re

KK = 1000.0  # OUTG11 energy unit is 1000 cm^-1


# orbital token n-l-occ: the occupation digit only counts if it is NOT itself
# the principal quantum number of a following orbital (i.e. not followed by an
# orbital letter). Using a lookahead keeps '3p' in '3s3p' from being read as
# occupation 3 of 3s.
_NL_OCC = re.compile(r"([1-9])([spdfghik])(\d(?![spdfghik]))?")


def _parity_from_config(config):
    """Parity = (-1)^(sum of l*occ) over the configuration string, e.g.
    '3s3p' -> odd, '3s2' -> even, '3s3d' -> even. Returns 'e' or 'o'."""
    lmap = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5, "i": 6}
    lsum = 0
    for n, l, occ in _NL_OCC.findall(config):
        o = int(occ) if occ else 1
        lsum += lmap[l] * o
    return "o" if (lsum % 2) else "e"


def parse_outg11(path):
    """Return (levels, lines).

    levels: list of dicts {config, term, J, parity, E_calc}  (E in cm^-1)
    lines : list of dicts {E_low, J_low, term_low, E_up, J_up, term_up,
                           lambda_A, loggf, GA}
    """
    with open(path, errors="replace") as f:
        text = f.read()
    return _parse_levels(text), _parse_lines(text)


# A basis row inside an EIGENVECTORS block:
#   "  <rownum>: <config>  (<parent>)<term>  <k>  c1 c2 ... cM"
# <rownum> = position in the printed list (can repeat across J-blocks),
# <k>      = the basis-state index, STABLE across the column-chunks of one block,
# c1..cM   = this basis state's coefficients in the eigenvectors of THIS chunk.
_BASIS_ROW = re.compile(
    r"^\s*\d+:\s+(\S+)\s+\((.*?)\)\s*(\d[A-Z]\*?)\s+(\d+)\s+(.*\S)\s*$")
# A column-chunk header: leading spaces, a lone chunk number, then config tokens.
_CHUNK_HDR = re.compile(r"^\s+\d+\s+[1-9][spdfghiklm]")
_NUM = re.compile(r"-?\d+\.\d+")


def parse_compositions(path):
    """Parse OUTG11's EIGENVALUES + EIGENVECTORS(LS) blocks and assign each
    computed level a ROBUST identity = the dominant basis state (config, term),
    i.e. the basis row with the largest |coefficient| in that eigenvalue's column.

    Returns dict (parity, Jstr) -> list (energy-ordered) of dicts
        {E_calc (cm^-1), config, term, parity}.

    The eigenvector matrix is printed in COLUMN-CHUNKS (~11 eigenvectors wide):
        EIGENVECTORS ( LS COUPLING)
          <chunk#>   <cfg> <cfg> ...        <- chunk header (column configs)
                     (parent)term ...       <- column parent-terms
          R: cfg (parent)term  k  c1 .. c11 <- one row per basis state, THIS chunk
          ...
          <chunk#>   <cfg> ...              <- next chunk: eigenvectors 12..22
          R: cfg (parent)term  k  c12 .. c22
          ...
        PURITY=...                          <- end of LS block
    Each basis state (keyed by k) recurs once per chunk; we concatenate its
    per-chunk coefficients into its full eigenvector row, then take, for every
    eigenvalue column, the basis state with the largest |coefficient|. Handles a
    single chunk (small bases) and many chunks (large bases) identically.
    """
    with open(path, errors="replace") as f:
        lines = f.readlines()
    out = {}
    i, n = 0, len(lines)
    while i < n:
        m = re.search(r"EIGENVALUES\s*\(J=\s*([-\d.]+)\)", lines[i])
        if not m:
            i += 1
            continue
        Jval = "%g" % float(m.group(1))
        # --- eigenvalues: numeric lines until the eigenvector/other section ---
        # RCG prints the eigenvalues 11-per-line with a BLANK LINE between each row
        # (see OUTG11: "EIGENVALUES (J=..)" then rows of 11 floats separated by
        # blanks, ended by "CONFIG. NO."). The block therefore must NOT terminate on
        # a blank line -- doing so stopped after the FIRST 11 values and silently
        # dropped the rest of the J-block (for J=2 that was 11 of 113 levels kept,
        # losing every high-lying doubly-excited level: 3p3d, 3d.nl, ...). Only the
        # explicit section headers (or the next EIGENVALUES block) end it.
        evs = []
        j = i + 1
        while j < n and not lines[j].strip():
            j += 1
        while j < n:
            s = lines[j]
            if any(t in s for t in ("CONFIG. NO.", "EIGENVECTORS", "G-VALUES")) \
                    or re.search(r"EIGENVALUES\s*\(J=", s):
                break
            nums = _NUM.findall(s)
            if nums:
                evs.extend(float(x) for x in nums)
            j += 1
        # --- locate the LS eigenvector block ---
        while j < n and not ("EIGENVECTORS" in lines[j] and "LS" in lines[j]):
            if re.search(r"EIGENVALUES\s*\(J=", lines[j]):
                break
            j += 1
        if j >= n or "EIGENVECTORS" not in lines[j]:
            i = j
            continue
        # --- read column-chunks, concatenating each basis state's coefficients ---
        # coeffs[k] = full list of this basis state's eigenvector coefficients
        # (across all chunks); label[k] = (config, term). col_base tracks how many
        # eigenvector columns precede the current chunk.
        coeffs, label = {}, {}
        col_base = 0
        chunk_width = 0          # columns seen in the current chunk
        k = j + 1
        while k < n:
            s = lines[k].rstrip("\n")
            if "PURITY" in s or re.search(r"EIGENVALUES\s*\(J=", s) \
                    or ("EIGENVECTORS" in s and "JJ" in s):
                break
            mr = _BASIS_ROW.match(s)
            if mr:
                cfg, term, bk = mr.group(1), mr.group(3), int(mr.group(4))
                cs = [float(x) for x in _NUM.findall(mr.group(5))]
                row = coeffs.setdefault(bk, [])
                # pad to the chunk's starting column, then append this chunk's cs
                if len(row) < col_base:
                    row.extend([0.0] * (col_base - len(row)))
                row.extend(cs)
                label[bk] = (cfg, term)
                chunk_width = max(chunk_width, len(cs))
            elif _CHUNK_HDR.match(s):
                # new chunk: advance the column base by the previous chunk width
                col_base += chunk_width
                chunk_width = 0
            k += 1
        # --- dominant basis state per eigenvalue column ---
        levs = []
        for col, E in enumerate(evs):
            best, bestc = None, -1.0
            for bk, cs in coeffs.items():
                if col < len(cs) and abs(cs[col]) > bestc:
                    bestc = abs(cs[col]); best = label[bk]
            if best:
                levs.append({"E_calc": E * KK, "config": best[0],
                             "term": best[1],
                             "parity": _parity_from_config(best[0])})
        if levs:
            out.setdefault((levs[0]["parity"], Jval), []).extend(levs)
        i = k
    for key in out:
        out[key].sort(key=lambda d: d["E_calc"])
    return out


def identify_lines(path, etol=2.0):
    """Attach ROBUST upper/lower level identities to each E1 line by matching
    the line's E_low/E_up (within etol cm^-1, same J) to the eigenvector-
    composition level table from parse_compositions(). This replaces the
    transition table's own config/term labels (the OUTG11 line block labels are
    derived from the block header, which is unreliable) with the dominant-
    eigenvector identity — the same identity the RCE level fit is built on.

    Returns the parse_outg11 line dicts, each augmented with:
        config_low, term_id_low, config_up, term_id_up   (or None if unmatched)
    Lines whose endpoints can't both be identified keep None and should be
    dropped from gf comparisons rather than matched by wavelength.
    """
    comp = parse_compositions(path)
    _, lines = parse_outg11(path)

    def _Jstr(J):
        try:
            return "%g" % float(J)
        except (ValueError, TypeError):
            return None

    def _find(E, J):
        levs = comp.get(("e", _Jstr(J))) or []
        levs = levs + (comp.get(("o", _Jstr(J))) or [])
        best, bd = None, etol
        for L in levs:
            d = abs(L["E_calc"] - E)
            if d <= bd:
                bd, best = d, L
        return best

    for d in lines:
        lo = _find(d["E_low"], d["J_low"])
        up = _find(d["E_up"], d["J_up"])
        d["config_low"] = lo["config"] if lo else None
        d["term_id_low"] = lo["term"] if lo else None
        d["config_up"] = up["config"] if up else None
        d["term_id_up"] = up["term"] if up else None
    return lines


def parse_levels1(path):
    """Parse RCE's LEVELS1 output -> list of dicts with the OBSERVED and FITTED
    energies together:
        {config, term, J, parity, E_obs, E_fit}   (energies in cm^-1)
    LEVELS1 row format (energies in kK):
        <E_obs> <E_fit> <residual> <gJ> <config> J= <J> <pct>% <conf> (<parent>) <term> ...
    The first numeric is observed, the second fitted. The dominant component's
    term (after the leading percentage) gives the label.
    """
    out = []
    with open(path, errors="replace") as f:
        for ln in f:
            s = ln.rstrip("\n")
            if not s.strip() or "J=" not in s:
                continue
            # leading: E_obs  E_fit  residual  gJ  config
            m = re.match(r"\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+"
                         r"(-?\d+\.\d+)\s+(\S+)\s+J=\s*([\d.]+)", s)
            if not m:
                continue
            # dominant component's term. LEVELS1 contains BOTH an LS-coupling
            # block (parent like '(2S) 3P') and a JJ-coupling block (parent like
            # '(2S 1/2) 1/2'). Keep only LS rows: the term token is <mult><L>
            # with NO slash in the parenthetical.
            mc = re.search(r"\d+%\s+\S+\s+\(([^)]*)\)\s*(\d[A-Z]\*?)", s)
            if not mc or "/" in mc.group(1):
                continue                       # skip JJ-coupling rows
            term = mc.group(2)
            E_obs = float(m.group(1)) * KK
            E_fit = float(m.group(2)) * KK
            config = m.group(5)
            J = m.group(6)
            parity = _parity_from_config(config)
            out.append({"config": config, "term": term, "J": J,
                        "parity": parity, "E_obs": E_obs, "E_fit": E_fit})
    return out


def _parse_levels(text):
    """Extract computed levels from the EIGENVALUES / EIGENVECTORS blocks.

    Each block looks like:
        0  EIGENVALUES      (J= 1.0)
                                       16.566   35.170
           CONFIG. NO.
                                         1        1
        ...
           EIGENVECTORS   (    LS COUPLING)
                            1          3s3p     3s3p
                                      (2S) 3P  (2S) 1P
    We read the eigenvalues for each J, then the term labels from the LS
    eigenvector header, and the config from the ENERGY MATRIX header above.
    """
    levels = []
    lines = text.splitlines()
    cur_config = None
    cur_J = None
    cur_evs = None  # pending eigenvalues for the current J, awaiting term labels
    for i, ln in enumerate(lines):
        m_em = re.search(r"ENERGY MATRIX.*COUPLING\)\s+J=\s*([-\d.]+)\s+(.*?)\s+CONFIG", ln)
        if m_em:
            cur_config = m_em.group(2).strip()
        m_ev = re.search(r"EIGENVALUES\s*\(J=\s*([-\d.]+)\)", ln)
        if m_ev:
            cur_J = m_ev.group(1)
            # gather eigenvalues from following non-blank numeric lines
            evs = []
            j = i + 1
            while j < len(lines):
                s = lines[j]
                if any(t in s for t in ("CONFIG. NO.", "EIGENVECTORS",
                                        "G-VALUES")):
                    break
                nums = re.findall(r"-?\d+\.\d+", s)
                if nums:
                    evs.extend(float(x) for x in nums)
                elif s.strip() == "" and evs:
                    break
                j += 1
            cur_evs = evs
            continue
        # the FIRST 'EIGENVECTORS (LS COUPLING)' after an EIGENVALUES block
        # carries the term labels matching cur_evs
        if cur_evs is not None and "EIGENVECTORS" in ln and "LS" in ln:
            blob = " ".join(lines[i + 1:i + 7])
            toks = re.findall(r"\(\s*\w+\s*\)\s*\d[A-Z]\*?", blob)
            toks = [re.sub(r"\s+", " ", t).strip() for t in toks]
            parity = _parity_from_config(cur_config or "")
            for k, ev in enumerate(cur_evs):
                levels.append({
                    "config": cur_config or "?",
                    "term": toks[k] if k < len(toks) else "?",
                    "J": cur_J,
                    "parity": parity,
                    "E_calc": ev * KK,
                })
            cur_evs = None
    return levels


def _parse_lines(text):
    """Parse the ELEC DIP SPECTRUM transition table."""
    out = []
    in_tab = False
    upper = None
    for ln in text.splitlines():
        if "ELEC DIP SPECTRUM" in ln and "ENERGIES" in ln:
            in_tab = True
            continue
        if not in_tab:
            continue
        # upper-level marker:
        #  " * * *  <E> <J> <idx> (parent)<term>  <#> <species> <config>  * * *"
        m_up = re.search(
            r"\* \* \*\s+(-?[\d.]+)\s+([\d.]+)\s+\d+\s+(\(.*?\)\s*\d[A-Z]\*?)"
            r"\s+\d+\s+\S.*?\s+(\S+)\s+\* \* \*", ln)
        if m_up:
            upper = {"E": float(m_up.group(1)) * KK, "J": m_up.group(2),
                     "term": re.sub(r"\s+", " ", m_up.group(3)).strip(),
                     "config": m_up.group(4).strip()}
            continue
        # transition rows: "<idx> <E_low> <J_low> <N> (parent)<term_low> ... lam loggf GA ..."
        m_row = re.match(
            r"\s*\d+\s+(-?[\d.]+)\s+([\d.]+)\s+\d+\s+(\(.*?\)\s*\d[A-Z]\*?)", ln)
        toks = ln.split()
        if m_row and len(toks) >= 8 and upper is not None:
            gi = next((i for i, t in enumerate(toks)
                       if re.match(r"^[\d.]+E[+-]\d+$", t)), None)
            if gi is None or gi < 3:
                continue
            try:
                lam = float(toks[gi - 2]); loggf = float(toks[gi - 1])
                GA = float(toks[gi])
                E_low = float(m_row.group(1)) * KK
                J_low = m_row.group(2)
                term_low = re.sub(r"\s+", " ", m_row.group(3)).strip()
            except (ValueError, IndexError):
                continue
            if 0 < lam < 1e6 and -15 < loggf < 5:
                out.append({"E_low": E_low, "J_low": J_low,
                            "term_low": term_low,
                            "E_up": upper["E"], "J_up": upper["J"],
                            "term_up": upper["term"],
                            "config_up": upper.get("config", "?"),
                            "lambda_A": lam, "loggf": loggf, "GA": GA})
    return out


if __name__ == "__main__":
    import sys
    lv, ln = parse_outg11(sys.argv[1])
    print(f"levels: {len(lv)}")
    for d in lv[:20]:
        print(f"  {d['config']:14s} {d['term']:10s} J={d['J']:4s} "
              f"par={d['parity']}  E={d['E_calc']:12.3f}")
    print(f"lines: {len(ln)}")
    for d in ln[:10]:
        print(f"  lam={d['lambda_A']:10.3f}  loggf={d['loggf']:+.3f}  "
              f"up={d['term_up']}")
    # eigenvector-composition identities (the robust labels); the lowest few are
    # an easy correctness spot-check -- the ground state must be physical, NOT
    # mislabelled (a multi-chunk parsing bug used to call it e.g. '3p4p').
    comp = parse_compositions(sys.argv[1])
    flat = sorted((L["E_calc"], L["config"], L["term"], J)
                  for (p, J), v in comp.items() for L in v)
    print(f"compositions: {len(flat)} levels; lowest 8 by energy:")
    for E, c, t, J in flat[:8]:
        print(f"  E={E:11.1f}  J={J:4s}  {c:10s} {t}")
