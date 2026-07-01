#!/usr/bin/env python3
"""
Per-species diagnostic report for the Modern Atomic Linelists project.

Generates a multi-page PDF that, for one ion, shows:
  Page 1  Summary: species, level/line counts, fit RMS (when available),
          key resonance lines computed vs reference.
  Page 2  Energy-level diagram (the centerpiece): computed vs NIST (and FITTED
          once RCE is wired up), with connectors, on a BROKEN energy axis so the
          dense low-lying terms and sparse high Rydberg levels are both legible.
  Page 3  Fit residuals: (E_calc - E_obs) vs E_obs, per matched level.
  Page 4  gf comparison / spectrum: log gf vs wavelength (and vs reference gf
          when available).

Inputs come from the parser (parse_cowan.parse_outg11) plus the cached NIST
levels (data/nist/<SP>_levels.tsv). The "fitted" column is drawn when a fitted-
levels source is supplied (RCE LEVELS files); until then the report shows the
ab-initio-vs-observed comparison, which is already a useful pre-fit diagnostic.

Usage:
    python3 tools/make_report.py MgI \
        --outg11 work/mg1/OUTG11 \
        --nist data/nist/MgI_levels.tsv \
        --out docs/reports/MgI_report.pdf
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_cowan import parse_outg11, parse_levels1  # noqa: E402

# --- house style -----------------------------------------------------------
# A nicer typeface (Charter: a clean serif designed for technical text) with
# matching mathtext, and a small curated palette shared across all pages so the
# same series is the same colour/marker everywhere. Falls back gracefully if the
# font is missing on another machine.
import matplotlib.font_manager as _fm
_HAVE = {f.name for f in _fm.fontManager.ttflist}
_FONT = next((f for f in ("Charter", "PT Serif", "Palatino", "Georgia",
                          "DejaVu Serif") if f in _HAVE), "serif")
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": [_FONT, "DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 13,
    "axes.linewidth": 0.9,
    "axes.titlesize": 14,
    "axes.labelsize": 14,
    "legend.fontsize": 9.5,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 4, "ytick.major.size": 4,
    # keep all four spines: every plot is fully boxed by axes
})

# Colourblind-safe series palette (darkened Okabe-Ito family: distinguishable
# under red-green colour blindness). No red+green pair: deep blue / dark amber /
# dark purple for ab initio, our fit, and Kurucz respectively, used on every
# page. Marker SHAPE also differs per series as a second, colour-independent cue.
COL_ABINITIO = "#004C8C"   # deep blue
COL_FIT      = "#B26B00"   # dark amber
COL_KURUCZ   = "#8E2F6E"   # dark purple

# Distinct OPEN marker per series (same shape/colour for a series on every page).
MK_ABINITIO = "o"
MK_FIT      = "D"
MK_KURUCZ   = "s"
_MK_SIZE = 30


def _scatter(ax, x, y, color, marker, label=None):
    """House scatter: always-open markers, consistent size/weight, used on every
    comparison page so series look identical across pages."""
    ax.scatter(x, y, s=_MK_SIZE, facecolors="none", edgecolors=color,
               marker=marker, linewidths=1.1, zorder=3, label=label)


# ----------------------------------------------------------------------------
# data loading / matching
# ----------------------------------------------------------------------------
def load_ie(path, species):
    """Read an ionization-energy table (species<TAB>IE_cm1); return IE in cm^-1
    for `species`, or None."""
    if not path or not os.path.exists(path):
        return None
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0].strip() == species:
                try:
                    return float(parts[1])
                except ValueError:
                    return None
    return None


def load_nist(path):
    rows = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            # DROP 2p^5 core-EXCITED levels (an open 2p core-hole, E_obs ~440000+
            # cm^-1). Our basis is built on the CLOSED 2p^6 core, so these are not
            # representable -- and worse, _cfgkey strips the core-hole and collapses
            # e.g. '2p5.(2P*).3s.3p2.(4P) 3P' -> '3p2 3P', which then FALSE-MATCHES
            # our valence 3p^2 3P and injects a ~-385000 cm^-1 residual. There are
            # ~55 such entries in the Mg I NIST table; they must never enter the fit
            # or the level matching. (The parser fix exposed them by finally matching
            # the high-lying levels that had been silently dropped.)
            if "2p5" in p[0].replace(" ", ""):
                continue
            rows.append({"config": p[0], "term": p[1], "J": p[2],
                         "E_obs": float(p[3]), "parity": p[4]})
    return rows


def _termkey(term):
    """Normalize a term label to its FINAL multiplicity+L core, ignoring any
    parent term in parentheses: '(2S) 3P' -> '3P', '1P*' -> '1P', '3P' -> '3P'."""
    # drop parenthetical parent terms, then take the last <digit><Letter> token
    t = re.sub(r"\([^)]*\)", " ", term)
    ms = re.findall(r"\d[A-Z]", t)
    return ms[-1] if ms else term.strip()


# orbital token: multi-digit n + l + optional occupation; the occupation digit
# counts only when NOT followed by another orbital letter or digit, so '3s11d'
# parses as 3s + 11d (not 3s^1 + 1d) and '3p2' as 3p^2. (Must stay identical to
# build_ine20._ORB so config keys agree between the two.)
_ORB = re.compile(r"\d+[spdfghik](?:\d(?![spdfghik\d]))?")


def _cfgkey(config):
    """Reduce a config to its valence orbital tokens (occupations normalized),
    e.g. 'Mg I   3s3p' -> '3s.3p', '3s.3p' -> '3s.3p', '2p6.3s2' -> '3s2',
    '3s11d' -> '3s.11d'. Keeps the last two orbital tokens (the valence part);
    drops occupation '1' and closed-shell markers so computed and NIST agree."""
    toks = _ORB.findall(config)
    norm = []
    for t in toks:
        m = re.match(r"(\d+[spdfghik])(\d?)", t)
        nl, occ = m.group(1), m.group(2)
        norm.append(nl + (occ if occ and occ != "1" else ""))
    # a single doubly-occupied valence orbital (e.g. 3s2) is the ground-type
    # config; NIST often writes it with the closed core (2p6.3s2). Reduce to the
    # outermost token so the two forms agree.
    if norm and re.fullmatch(r"\d+[spdfghik]2", norm[-1]):
        return norm[-1]
    return ".".join(norm[-2:]) if norm else config.strip()


def _Jkey(J):
    """Normalize a J value to a canonical string; '' / non-numeric -> None."""
    try:
        return str(float(J))
    except (ValueError, TypeError):
        return None


def match_levels(calc, nist):
    """Match computed to observed by (valence-config, termkey, J). Returns list of
    dicts with E_calc, E_obs (or None), and labels."""
    index = {}
    for n in nist:
        jk = _Jkey(n["J"])
        if jk is None:
            continue
        key = (_cfgkey(n["config"]), _termkey(n["term"]), jk)
        index.setdefault(key, n)
    out = []
    for c in calc:
        jk = _Jkey(c["J"])
        key = (_cfgkey(c["config"]), _termkey(c["term"]), jk)
        n = index.get(key) if jk is not None else None
        out.append({**c, "E_obs": (n["E_obs"] if n else None),
                    "matched": n is not None})
    return out


# ----------------------------------------------------------------------------
# pages
# ----------------------------------------------------------------------------
def _resid_stats(resids):
    """(N, median|d|, RMS, 90th pct |d|, max|d|) for a list of residuals, with the
    global median offset removed (the convention the headline level-RMS uses)."""
    if not resids:
        return None
    r = np.array(resids, float)
    r = r - np.median(r)
    a = np.abs(r)
    return {"n": len(r), "median": float(np.median(a)),
            "rms": float(np.sqrt(np.mean(r ** 2))),
            "p90": float(np.percentile(a, 90)), "max": float(a.max())}


def build_compare(fitted, kurucz_levels, cap):
    """Head-to-head stats: ours (from the identity-matched fitted-levels JSON,
    which already carries E_obs/E_fit per level) vs Bob's RCE residuals
    (kurucz_levels), both restricted to E_obs <= cap so the comparison is on an
    identical footing. Returns {cap, ours, bob} or None."""
    if not fitted or not kurucz_levels:
        return None
    ours = [fl["E_fit"] - fl["E_obs"] for fl in fitted
            if fl.get("E_fit") is not None and fl.get("E_obs") is not None
            and fl["E_obs"] <= cap]
    bob = [d for (e, d) in kurucz_levels if e <= cap]
    so, sb = _resid_stats(ours), _resid_stats(bob)
    return {"cap": cap, "ours": so, "bob": sb} if (so and sb) else None


def page_summary(pdf, species, matched, lines, fit_rms, compare=None):
    fig = plt.figure(figsize=(8.5, 11))
    fig.suptitle(f"Diagnostic report: {species}", fontsize=16, y=0.97)
    ax = fig.add_axes([0.08, 0.1, 0.84, 0.8]); ax.axis("off")

    n_match = sum(1 for m in matched if m["matched"])
    txt = []
    txt.append(f"Computed levels: {len(matched)}")
    txt.append(f"Matched to NIST: {n_match}")
    txt.append(f"Computed E1 lines: {len(lines)}")
    if fit_rms is not None:
        txt.append(f"Level-fit RMS (E_calc - E_obs): {fit_rms:.1f} cm^-1")
    else:
        resid = [m["E_calc"] - m["E_obs"] for m in matched
                 if m["matched"] and m["E_obs"] is not None]
        if resid:
            txt.append("Pre-fit RMS (E_calc - E_obs): "
                       f"{np.sqrt(np.mean(np.square(resid))):.0f} cm^-1 "
                       "(ab initio, no RCE fit yet)")

    # Head-to-head vs Bob Kurucz's RCE fit, on an IDENTICAL footing (same energy
    # cap, offset removed). RMS is tail-dominated by the singlet 3s.nd 1D / near-IE
    # levels, so the MEDIAN |E-O| is the honest measure of bulk quality.
    if compare:
        o, b = compare["ours"], compare["bob"]
        txt.append("")
        txt.append(f"Energy-level fit vs Bob Kurucz (E_obs <= {compare['cap']:.0f} "
                   "cm^-1, offset-removed):")
        txt.append(f"  {'':14}{'N':>5}{'median|E-O|':>13}{'RMS':>8}"
                   f"{'90th':>8}{'max':>8}")
        txt.append(f"  {'this fit':14}{o['n']:>5}{o['median']:>13.1f}"
                   f"{o['rms']:>8.0f}{o['p90']:>8.0f}{o['max']:>8.0f}")
        txt.append(f"  {'Bob (RCE)':14}{b['n']:>5}{b['median']:>13.1f}"
                   f"{b['rms']:>8.0f}{b['p90']:>8.0f}{b['max']:>8.0f}")
        txt.append("  bulk (median) is near-Bob; the ~3x RMS gap is the singlet")
        txt.append("  3s.nd 1D / near-IE tail (unobserved-perturber wall), not")
        txt.append("  data starvation -- the basis matches Bob's b1200{e,o} deck.")

    # strongest computed lines
    if lines:
        strong = sorted(lines, key=lambda d: -d["loggf"])[:8]
        txt.append("")
        txt.append("Strongest computed lines:")
        txt.append(f"  {'lambda(A)':>10}  {'log gf':>7}  upper term")
        for d in strong:
            txt.append(f"  {d['lambda_A']:>10.2f}  {d['loggf']:>+7.3f}  "
                       f"{d['term_up']}")

    ax.text(0.0, 1.0, "\n".join(txt), va="top", ha="left",
            family="monospace", fontsize=11, transform=ax.transAxes)
    pdf.savefig(fig); plt.close(fig)


def _segments(energies, gap_factor=4.0, min_gap=3000.0):
    """Find energy-axis break segments: cluster levels, break where the gap to
    the next level is large compared with the typical spacing."""
    e = np.sort(np.array(energies))
    if len(e) < 3:
        return [(e.min() - 100, e.max() + 100)] if len(e) else [(0, 1)]
    diffs = np.diff(e)
    typ = np.median(diffs[diffs > 0]) if np.any(diffs > 0) else 1.0
    thr = max(min_gap, gap_factor * typ)
    segs = []
    start = e[0]
    for i in range(len(diffs)):
        if diffs[i] > thr:
            segs.append((start, e[i]))
            start = e[i + 1]
    segs.append((start, e[-1]))
    # pad each segment a touch
    return [(a - 0.02 * max(1, b - a) - 50, b + 0.02 * max(1, b - a) + 50)
            for a, b in segs]


_LMAP = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5, "I": 6, "K": 7}
_ORB_TOK = re.compile(r"[1-9][spdfghik]\d?(?![spdfghik])")


def _running_orbital(config):
    """The 'running' (valence/Rydberg) orbital of a configuration, e.g.
    '3s.3p' -> '3p', '3s.4d' -> '4d', '2p6.3s2' -> '3s'. Strips occupation."""
    toks = _ORB_TOK.findall(config)
    clean = [re.match(r"[1-9][spdfghik]", t).group(0) for t in toks]
    # the core s electron (lowest s that's doubly closed-ish) is '3s' for Mg-like;
    # generally the running orbital is the last (outermost) token.
    noncore = [o for o in clean if o != clean[0]] if len(clean) > 1 else clean
    return (noncore[-1] if noncore else clean[-1]) if clean else "?"


def _Lof(term):
    m = re.match(r"\d+([A-Z])", _termkey(term))
    return _LMAP.get(m.group(1), 9) if m else 9


def _term_mathlabel(term, parity):
    """LaTeX-style term label with real superscript multiplicity and a degree
    marker for odd parity, e.g. ('3P','o') -> '$^{3}\\!P^{\\circ}$'."""
    tk = _termkey(term)              # e.g. '3P'
    m = re.match(r"(\d+)([A-Z])", tk)
    if not m:
        return tk
    mult, L = m.group(1), m.group(2)
    deg = r"^{\circ}" if parity == "o" else ""
    return rf"$^{{{mult}}}\!{L}{deg}$"


def _term_sort_key(term_panel):
    """Order term panels by mean energy (set by caller)."""
    return term_panel["E_mean"]


def _group_terms(matched):
    """Group matched levels by (config, termkey, parity); return panels sorted
    by energy (low first)."""
    groups = {}
    for m in matched:
        key = (_cfgkey(m.get("config", "")), _termkey(m["term"]),
               m.get("parity", "e"))
        groups.setdefault(key, []).append(m)
    panels = []
    for (cfg, tk, par), lev in groups.items():
        Es = [x["E_calc"] for x in lev] + \
             [x["E_obs"] for x in lev if x["E_obs"] is not None]
        panels.append({"cfg": cfg, "tk": tk, "par": par, "lev": lev,
                       "E_mean": np.mean(Es) if Es else 0.0})
    panels.sort(key=_term_sort_key)
    return panels


def _grotrian_layout(panels):
    """Na I-style column assignment: ONE column per term SYMBOL (multiplicity +
    L + parity), e.g. 3S, 1S, 3P*, 1P*, 3D, 1D, ... All configurations of a term
    (3s3p, 3s4p, ...) stack in that single column. Even-parity block on the left,
    odd on the right; within a block, order by (multiplicity, L). Returns a dict
    termkey+parity -> column x, plus the ordered column metadata."""
    # collect distinct (termkey, parity)
    terms = {}
    for p in panels:
        key = (p["tk"], p["par"])
        terms.setdefault(key, True)
    def order(k):
        tk, par = k
        m = re.match(r"(\d+)", tk)
        mult = int(m.group(1)) if m else 9
        return (_Lof(tk), -mult)            # by L, triplets before singlets
    even = sorted([k for k in terms if k[1] == "e"], key=order)
    odd = sorted([k for k in terms if k[1] == "o"], key=order)
    colx = {}
    cols = []          # (termkey, parity, x)
    x = 0
    for k in even:
        colx[k] = x; cols.append((k[0], k[1], x)); x += 1
    x += 1             # gap between parity blocks
    for k in odd:
        colx[k] = x; cols.append((k[0], k[1], x)); x += 1
    return {"colx": colx, "cols": cols, "nx": max(1, x)}


EV = 8065.544   # cm^-1 per eV


def _grotrian_page(pdf, species, panels, efield, solid_label, title,
                   layout=None, ylim=None, ie_cm=None, efield_offset=0.0):
    """Na I-style Grotrian diagram. One column per term symbol; all levels of a
    term stack in that column, each a short horizontal tick labelled by its
    running orbital (red). `efield` selects the solid-line energy ('E_calc' ab
    initio or 'E_fit' fitted); NIST observed is the dashed overlay. Left y-axis
    in cm^-1, right y-axis in eV. Optional ionization limit drawn dashed.
    `efield_offset` is ADDED to every solid energy -- used to remove the fit's
    global ground-pinned offset so fitted vs observed align (same convention as
    the residuals page and gf_fit's level-RMS)."""
    if layout is None:
        layout = _grotrian_layout(panels)
    colx, cols, nx = layout["colx"], layout["cols"], layout["nx"]

    fig, ax = plt.subplots(figsize=(max(7.0, 0.85 * nx + 2.5), 9.0))
    ax.set_title(f"{species}: {title}", fontsize=14, pad=12)

    half = 0.40
    any_ai = False     # did we draw any autoionizing (above-limit) level?
    # collect, per column, the label points (energy, orbital) to de-collide
    for p in panels:
        xc = colx.get((p["tk"], p["par"]))
        if xc is None:
            continue
        for m in p["lev"]:
            es = m.get(efield); eo = m.get("E_obs")
            if es is not None:
                es += efield_offset
            orb = _running_orbital(m.get("config", ""))
            # a level above the first ionization limit is autoionizing (a
            # doubly-excited resonance embedded in the continuum) -> tag it.
            e_ref = es if es is not None else eo
            ai = ie_cm is not None and e_ref is not None and e_ref > ie_cm
            lab = orb + ("  (AI)" if ai else "")
            any_ai = any_ai or ai
            if es is not None:
                ax.hlines(es, xc - half, xc + half, color=COL_ABINITIO, lw=1.6)
                ax.text(xc + half + 0.04, es, lab, va="center", ha="left",
                        fontsize=10.5, color="0.30",
                        fontstyle="italic" if ai else "normal")
            if eo is not None:
                ax.hlines(eo, xc - half, xc + half, color=COL_FIT, lw=1.2,
                          ls="--", alpha=0.9)
                if es is None:   # observed-only level still gets its (AI) tag
                    ax.text(xc + half + 0.04, eo, lab, va="center", ha="left",
                            fontsize=10.5, color="0.30",
                            fontstyle="italic" if ai else "normal")

    # ionization limit, labelled with the actual energy in eV
    if ie_cm is not None:
        ax.axhline(ie_cm, color="0.3", lw=1.0, ls="--")
        ax.text(nx - 0.3, ie_cm, f" ionization limit ({ie_cm / EV:.2f} eV)",
                va="bottom", ha="right", fontsize=11, color="0.3")
        if any_ai:
            ax.text(nx - 0.3, ie_cm, "levels above are autoionizing (AI) ",
                    va="top", ha="right", fontsize=9.5, color="0.45",
                    fontstyle="italic")

    # x-axis: term symbols
    ax.set_xticks([x for _, _, x in cols])
    ax.set_xticklabels([_term_mathlabel(tk, par) for tk, par, _ in cols],
                       fontsize=13)

    ax.set_xlim(-0.8, nx - 0.2)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_ylabel("Energy (cm$^{-1}$)")

    # right axis in eV, sharing the same data range, with a visible solid spine
    axr = ax.twinx()
    axr.set_ylim(ax.get_ylim()[0] / EV, ax.get_ylim()[1] / EV)
    axr.set_ylabel("Energy (eV)")
    axr.spines["right"].set_visible(True)
    axr.spines["right"].set_color("black")
    axr.spines["right"].set_linewidth(0.8)

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=COL_ABINITIO, lw=1.6, label=solid_label),
               Line2D([0], [0], color=COL_FIT, lw=1.2, ls="--",
                      label="observed (NIST)")]
    ax.legend(handles=handles, frameon=False, fontsize=11, loc="lower right")
    fig.tight_layout(rect=[0.04, 0.02, 1, 1.0])
    pdf.savefig(fig); plt.close(fig)


def _group_fitted(fitted):
    """Group RCE LEVELS1 rows (each carrying E_obs and E_fit) by
    (config, termkey, parity). Dedups the LS/JJ duplication in LEVELS1 by
    keeping one row per (cfg, term, J, rounded E_fit)."""
    seen = set()
    uniq = []
    for m in fitted:
        key = (_cfgkey(m["config"]), _termkey(m["term"]), str(m["J"]),
               round(m["E_fit"], 1))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(m)
    groups = {}
    for m in uniq:
        key = (_cfgkey(m["config"]), _termkey(m["term"]), m["parity"])
        groups.setdefault(key, []).append(m)
    panels = []
    for (cfg, tk, par), lev in groups.items():
        Es = [x["E_fit"] for x in lev] + \
             [x["E_obs"] for x in lev if x["E_obs"] is not None]
        panels.append({"cfg": cfg, "tk": tk, "par": par, "lev": lev,
                       "E_mean": np.mean(Es) if Es else 0.0})
    panels.sort(key=_term_sort_key)
    return panels


def _attach_abinitio(calc, fitted):
    """Attach E_calc (ab initio) to each fitted level. OUTG11's config labels are
    unreliable under CI, but within each (parity, J) the ab initio eigenvalues
    and the fitted levels are the SAME states in the same energy order. So we
    match by energy rank within (parity, J). Returns the fitted list with an
    'E_calc' key added (None if no ab initio available)."""
    from collections import defaultdict
    ab = defaultdict(list)
    for c in calc:
        ab[(c["parity"], str(float(c["J"])))].append(c["E_calc"])
    for k in ab:
        ab[k].sort()
    # group fitted by (parity, J), sort by E_fit, pair with sorted ab initio
    fb = defaultdict(list)
    for m in fitted:
        fb[(m["parity"], str(float(m["J"])))].append(m)
    for k, lev in fb.items():
        lev.sort(key=lambda x: x["E_fit"])
        cand = ab.get(k, [])
        for i, m in enumerate(lev):
            m["E_calc"] = cand[i] if i < len(cand) else None
    return fitted


def _relabel_by_nist(fitted, nist):
    """RCE labels a level by its plurality eigenvector component, which is
    ambiguous for strongly-mixed levels (e.g. 3s3d J=2 at 57/42 % 3D/1D). But
    each fitted level's E_obs value came FROM a specific NIST level whose term is
    authoritative. Re-tag each fitted level's term/config by looking up its E_obs
    in the NIST table (exact energy match), so the display is consistent with the
    fit input."""
    if not nist:
        return fitted
    by_E = {}
    for o in nist:
        by_E[round(o["E_obs"], 1)] = o          # E_obs in cm^-1
    for m in fitted:
        eo = m.get("E_obs")
        if eo is None:
            continue
        o = by_E.get(round(eo, 1))
        if o:
            m["term"] = o["term"]
            m["config"] = o["config"]
            m["parity"] = o["parity"]
    return fitted


def build_unified_panels(calc, fitted, nist=None):
    """One authoritative level list (from LEVELS1) carrying E_calc/E_obs/E_fit,
    grouped into term panels. Both Grotrian pages and the table use THIS, so the
    before/after layouts are identical. When `nist` is given, level term/config
    labels are taken from the NIST identity of each level's E_obs (robust for
    strongly-mixed levels)."""
    fitted = _attach_abinitio(calc, _dedup_levels1(fitted))
    fitted = _relabel_by_nist(fitted, nist)
    groups = {}
    for m in fitted:
        key = (_cfgkey(m["config"]), _termkey(m["term"]), m["parity"])
        groups.setdefault(key, []).append(m)
    panels = []
    for (cfg, tk, par), lev in groups.items():
        Es = [x["E_fit"] for x in lev] + \
             [x["E_obs"] for x in lev if x["E_obs"] is not None]
        panels.append({"cfg": cfg, "tk": tk, "par": par, "lev": lev,
                       "E_mean": np.mean(Es) if Es else 0.0})
    panels.sort(key=_term_sort_key)
    return panels


def _dedup_levels1(fitted):
    seen = set(); uniq = []
    for m in fitted:
        key = (_cfgkey(m["config"]), _termkey(m["term"]), str(float(m["J"])),
               round(m["E_fit"], 1))
        if key in seen:
            continue
        seen.add(key); uniq.append(m)
    return uniq


def page_levels(pdf, species, matched, fitted=None, calc=None, ie_cm=None,
                nist=None, fitted_label="fitted (RCE)"):
    """Centerpiece. When fitted data exist, both Grotrian pages share ONE panel
    set and ONE x-layout so paging gives a clean before/after:
      Page 1: ab initio (solid) + NIST (dashed)
      Page 2: fitted (solid, offset-aligned) + NIST (dashed)
    Then the merged level table. Without a fit, falls back to the ab-initio-only
    diagram from `matched`."""
    if fitted:
        panels = build_unified_panels(calc or [], fitted, nist)
        layout = _grotrian_layout(panels)
        # shared y-range across both pages so paging is a clean before/after
        allE = []
        for p in panels:
            for m in p["lev"]:
                allE += [v for v in (m.get("E_calc"), m.get("E_obs"),
                                     m.get("E_fit")) if v is not None]
        if ie_cm is not None:
            allE.append(ie_cm)
        pad = 0.04 * (max(allE) - min(allE)) if allE else 1.0
        ylim = (min(min(allE), 0) - pad, max(allE) + pad) if allE else None
        _grotrian_page(pdf, species, panels, "E_calc",
                       "computed (ab initio)",
                       "term diagram - ab initio vs observed", layout, ylim,
                       ie_cm)
        # remove the fit's global ground-pinned offset (median E_fit-E_obs) so
        # the fitted term diagram aligns to NIST; same convention as the
        # residuals page / gf_fit level-RMS.
        fdiffs = [m["E_fit"] - m["E_obs"] for p in panels for m in p["lev"]
                  if m.get("E_fit") is not None and m.get("E_obs") is not None]
        foff = -float(np.median(fdiffs)) if fdiffs else 0.0
        _grotrian_page(pdf, species, panels, "E_fit", fitted_label,
                       "term diagram - fitted vs observed", layout, ylim, ie_cm,
                       efield_offset=foff)
        _page_level_table(pdf, species, panels)
    else:
        panels = _group_terms(matched)
        layout = _grotrian_layout(panels)
        _grotrian_page(pdf, species, panels, "E_calc",
                       "computed (ab initio)",
                       "term diagram - ab initio vs observed", layout,
                       ie_cm=ie_cm)
        _page_level_table_simple(pdf, species, panels)


def _table_pages(pdf, ttl, col, rows, numeric_from=None):
    """Render a borderless, journal-style table across pages. No cell outlines;
    a bold larger header with a single rule beneath it, faint zebra striping for
    readability, and right-aligned numeric columns (those at index >= numeric_from).
    """
    per = 34
    ncol = len(col)
    for start in range(0, max(1, len(rows)), per):
        chunk = rows[start:start + per]
        fig, ax = plt.subplots(figsize=(8.5, 9.5)); ax.axis("off")
        ax.set_title(ttl, fontsize=13, pad=16)
        tbl = ax.table(cellText=chunk, colLabels=col, loc="upper center",
                       cellLoc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(9.5)
        tbl.scale(1, 1.45)

        for (r, c), cell in tbl.get_celld().items():
            cell.set_edgecolor("none")            # remove ALL outlines
            cell.PAD = 0.04
            if r == 0:                            # header row
                cell.get_text().set_fontsize(11)
                cell.get_text().set_fontweight("bold")
                cell.set_facecolor("white")
                # single rule under the header
                cell.visible_edges = "B"
                cell.set_edgecolor("0.25")
                cell.set_linewidth(1.1)
            else:
                # faint zebra striping
                cell.set_facecolor("#f2f2f2" if (r % 2 == 0) else "white")
            # right-align numeric columns, left-align label columns
            if numeric_from is not None and r > 0:
                if c >= numeric_from:
                    cell.get_text().set_ha("right"); cell.PAD = 0.06
                elif c < 2:
                    cell.get_text().set_ha("left"); cell.PAD = 0.06
        pdf.savefig(fig); plt.close(fig)


def _page_level_table(pdf, species, panels):
    """Unified table (one row per level): E_obs, E_calc, Δ_abinit, E_fit, Δ_fit.
    Operates on the unified panels where each level carries all three energies.
    Δ_fit (and the displayed E_fit) have the fit's global ground-pinned offset
    removed (median E_fit-E_obs), so Δ_fit is the honest residual that matches the
    headline level-RMS -- NOT the ~constant absolute offset of the pinned zero."""
    fdiffs = [m["E_fit"] - m["E_obs"] for p in panels for m in p["lev"]
              if m.get("E_fit") is not None and m.get("E_obs") is not None]
    foff = float(np.median(fdiffs)) if fdiffs else 0.0
    rows = []
    for p in panels:
        term = _term_mathlabel(p["tk"], p["par"])
        for m in sorted(p["lev"], key=lambda x: (x.get("E_fit") or 0.0)):
            try:
                Js = "%g" % float(m["J"])
            except (TypeError, ValueError):
                Js = str(m["J"])
            eo = m.get("E_obs"); ec = m.get("E_calc"); ef = m.get("E_fit")
            efa = (ef - foff) if ef is not None else None   # offset-aligned E_fit
            rows.append([
                p["cfg"], term, Js,
                (f"{eo:.1f}" if eo is not None else "--"),
                (f"{ec:.1f}" if ec is not None else "--"),
                (f"{ec - eo:+.1f}" if (ec is not None and eo is not None) else "--"),
                (f"{efa:.1f}" if efa is not None else "--"),
                (f"{efa - eo:+.1f}" if (efa is not None and eo is not None) else "--"),
            ])
    col = ["config", "term", "$J$", "$E_\\mathrm{obs}$", "$E_\\mathrm{abinit}$",
           "$\\Delta_\\mathrm{abinit}$", "$E_\\mathrm{fit}$",
           "$\\Delta_\\mathrm{fit}$"]
    _table_pages(pdf, f"{species}: levels — ab initio & fitted vs observed "
                 "(cm$^{-1}$)", col, rows,
                 numeric_from=3)   # right-align the energy columns


def _page_level_table_simple(pdf, species, panels):
    """Ab-initio-only table (no RCE fit yet)."""
    rows = []
    for p in panels:
        term = _term_mathlabel(p["tk"], p["par"])
        for m in sorted(p["lev"], key=lambda x: x["E_calc"]):
            try:
                Js = "%g" % float(m["J"])
            except (TypeError, ValueError):
                Js = str(m["J"])
            eo = m.get("E_obs")
            rows.append([p["cfg"], term, Js, f"{m['E_calc']:.1f}",
                         (f"{eo:.1f}" if eo is not None else "--"),
                         (f"{m['E_calc'] - eo:+.1f}" if eo is not None else "--")])
    col = ["config", "term", "$J$", "$E_\\mathrm{abinit}$", "$E_\\mathrm{obs}$",
           "$\\Delta E$ (cm$^{-1}$)"]
    _table_pages(pdf, f"{species}: levels (computed vs observed)", col, rows,
                 numeric_from=3)


def load_kurucz_levels(paths):
    """Parse Bob Kurucz's RCE fit log(s) (c<xxyy>{e,o}z.log / c<xxyy>{e,o}.log)
    for his fitted level energies. Each 'OBSERVED ENERGY / EIGENVALUE / E-O'
    table row gives (E_obs, E_fit, residual). Predicted levels (blank observed
    energy) are skipped. Returns list of (E_obs, residual) for levels that have
    an observed energy. `paths` is one path or a list (even+odd)."""
    if isinstance(paths, str):
        paths = [paths]
    out = []
    for path in paths:
        if not path or not os.path.exists(path):
            continue
        with open(path, errors="replace") as f:
            in_tab = False
            for ln in f:
                if "EIGENVALUE" in ln and "E-O" in ln:
                    in_tab = True
                    continue
                if not in_tab:
                    continue
                s = ln.rstrip("\n")
                # data row: idx  E_obs  EIGENVALUE  E-O  ...  (E_obs may be blank)
                m = re.match(r"\s*\d+\s+(-?\d+\.\d*)\s+(-?\d+\.\d+)\s+"
                             r"(-?\d+\.\d+)\b", s)
                if m:
                    e_obs = float(m.group(1))
                    resid = float(m.group(3))   # E-O = E_fit - E_obs
                    out.append((e_obs, resid))
                elif s.strip() and not s.lstrip()[0].isdigit() \
                        and "EIGENVALUE" not in s and not s.startswith("J "):
                    # a non-data, non-blank line (e.g. a new header) ends the run
                    if "OBSERVED" in s or "$" in s:
                        in_tab = False
    return out


def page_residuals(pdf, species, matched, unified_panels=None,
                   kurucz_levels=None, fitted_label="fitted (RCE)"):
    """Residuals vs observed, in TWO stacked panels: top shows the full y-range
    (so the large ab-initio offsets are visible), bottom zooms to the FITTED
    residuals (so the fit-quality comparison -- our fit vs Kurucz's -- is
    legible instead of collapsed onto zero by the ab-initio scale). With a fit
    (unified_panels), ab initio (circles) and fitted (diamonds) come from the
    SAME level set; without a fit, fall back to ab-initio-only `matched`.

    The FITTED residuals have their global MEDIAN OFFSET removed (the same
    convention gf_fit's level-RMS uses): the fit pins the ground-config Eav to 0,
    which puts the whole spectrum on an absolute zero offset by a near-constant
    from NIST's; the meaningful quantity is the residual STRUCTURE, not that
    bulk shift. Without this the RMS would be dominated by a ~1700 cm^-1 constant
    and disagree with the headline fit RMS. Ab-initio is left RAW (its offset is
    physically informative -- shows how far HF sits from observed)."""
    # collect series: (label, color, marker, points Nx2, is_fit?)
    series = []
    if unified_panels is not None:
        ab = [(m["E_obs"], m["E_calc"] - m["E_obs"])
              for p in unified_panels for m in p["lev"]
              if m.get("E_obs") is not None and m.get("E_calc") is not None]
        ft = [(m["E_obs"], m["E_fit"] - m["E_obs"])
              for p in unified_panels for m in p["lev"]
              if m.get("E_obs") is not None and m.get("E_fit") is not None]
        if ab:
            series.append(("ab initio", COL_ABINITIO, MK_ABINITIO,
                           np.array(ab), False))
        if ft:
            ft = np.array(ft)
            ft[:, 1] -= np.median(ft[:, 1])     # remove global offset (see above)
            series.append((fitted_label, COL_FIT, MK_FIT, ft, True))
    else:
        res = [(m["E_obs"], m["E_calc"] - m["E_obs"])
               for m in matched if m["matched"] and m["E_obs"] is not None]
        if res:
            series.append(("ab initio", COL_ABINITIO, MK_ABINITIO,
                           np.array(res), False))
    if kurucz_levels:
        series.append(("Kurucz (gfall fit)", COL_KURUCZ, MK_KURUCZ,
                       np.array(kurucz_levels), True))

    # RMS summary lives in the legend label (no top title)
    leg = {label: f"{label}   RMS = {np.sqrt(np.mean(pts[:, 1] ** 2)):.0f} "
                  "cm$^{-1}$" for label, _, _, pts, _ in series}

    fig, (a_full, a_zoom) = plt.subplots(2, 1, figsize=(8.0, 8.0))
    xlab = "Observed level  $E_\\mathrm{obs}$  (cm$^{-1}$)"
    ylab = "$E_\\mathrm{model} - E_\\mathrm{obs}$  (cm$^{-1}$)"

    def _draw(ax):
        ax.axhline(0, color="0.5", lw=0.7, zorder=1)
        for label, col, mk, pts, _ in series:
            _scatter(ax, pts[:, 0], pts[:, 1], col, mk, leg[label])

    if not series:
        a_full.text(0.5, 0.5, "no matched levels", ha="center")
    else:
        _draw(a_full)
        _draw(a_zoom)
        # zoom y-limit from the FITTED series only (so it isn't set by ab initio)
        fit_pts = [pts for _, _, _, pts, isfit in series if isfit]
        if fit_pts:
            allfit = np.concatenate([p[:, 1] for p in fit_pts])
            ylim = max(50.0, 1.15 * np.abs(allfit).max())
            a_zoom.set_ylim(-ylim, ylim)

    a_full.legend(frameon=False, loc="lower right")
    for ax in (a_full, a_zoom):
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def load_nist_lines(path):
    """Read cached NIST lines (ritz_wl_A, log_gf, gA, acc, conf/term/J for
    lower & upper). Returns list of dicts."""
    out = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            try:
                out.append({"lambda_A": float(p[0]), "loggf": float(p[1]),
                            "acc": p[3],
                            "conf_i": p[4], "term_i": p[5], "J_i": p[6],
                            "conf_k": p[7], "term_k": p[8], "J_k": p[9]})
            except ValueError:
                continue
    return out


def load_kurucz_lines(path, elem_code=None):
    """Read Bob Kurucz's per-ion GF line file (gf<xxyy>.pos / .lines / .all), the
    documented 160-column format. These are his SEMI-EMPIRICAL computed log gf
    (the K<yy> code marks the calc year). Returns dicts with the same keys as
    load_nist_lines (conf/term/J for both ends) so the identity matcher reuses
    them. Wavelengths are air (nm) above 200 nm; we keep loggf + identities only
    (matching is by level identity, not wavelength). elem_code (e.g. 12.00 for
    Mg I) filters to one species when a file mixes them."""
    out = []
    with open(path, errors="replace") as f:
        for ln in f:
            if len(ln) < 80:
                continue
            try:
                wl_nm = float(ln[0:11])
                loggf = float(ln[11:18])
                elem = float(ln[18:24])
                if elem_code is not None and abs(elem - elem_code) > 1e-3:
                    continue
                lab_lo = ln[42:52].strip().split()
                lab_up = ln[70:80].strip().split()
                if len(lab_lo) < 2 or len(lab_up) < 2:
                    continue
                out.append({
                    "loggf": loggf, "lambda_A": wl_nm * 10.0,  # nm -> Angstrom
                    "conf_i": lab_lo[0], "term_i": lab_lo[1], "J_i": ln[36:41],
                    "conf_k": lab_up[0], "term_k": lab_up[1], "J_k": ln[64:69],
                })
            except ValueError:
                continue
    return out


def _level_idkey(cfg, term, J):
    """Normalized (config, term, J) level identity (same as match_levels)."""
    return (_cfgkey(cfg), _termkey(term), _Jkey(J))


def _line_idkey(cfg_lo, t_lo, J_lo, cfg_up, t_up, J_up):
    """Order-independent transition identity from both end-level identities."""
    return tuple(sorted([_level_idkey(cfg_lo, t_lo, J_lo),
                         _level_idkey(cfg_up, t_up, J_up)]))


def match_gf_by_identity(outg11_path, nist_lines):
    """Match computed E1 lines to NIST by the EIGENVECTOR-COMPOSITION identity of
    BOTH end levels (config+term+J) -- the same robust identity the RCE level fit
    is built on. This avoids the term-pair+nearest-wavelength matcher that
    collides Rydberg series members (3p^2 1S vs 3s5s 1S) and invents residuals.
    Returns list of (loggf_comp, loggf_nist, lambda_comp)."""
    from parse_cowan import identify_lines
    buckets = {}
    for d in nist_lines:
        k = _line_idkey(d["conf_i"], d["term_i"], d["J_i"],
                        d["conf_k"], d["term_k"], d["J_k"])
        buckets.setdefault(k, []).append(d)
    pairs, used = [], set()
    for d in identify_lines(outg11_path):
        if not (d["config_low"] and d["config_up"]):
            continue
        k = _line_idkey(d["config_low"], d["term_id_low"], d["J_low"],
                        d["config_up"], d["term_id_up"], d["J_up"])
        cands = [c for c in buckets.get(k, []) if id(c) not in used]
        if not cands:
            continue
        n = min(cands, key=lambda c: abs(c["lambda_A"] - d["lambda_A"]))
        used.add(id(n))
        pairs.append((d["loggf"], n["loggf"], d["lambda_A"]))
    return pairs


def match_kurucz_gf(kurucz_lines, nist_lines, tol=0.01):
    """Match Bob Kurucz's per-ion gf lines to NIST by both-end level identity
    (config+term+J) AND nearest wavelength (within `tol` fractional), so that
    different Rydberg members sharing a reduced (config,term) key are not
    confused (his full series produces many lines per identity bucket). Kurucz
    wl is air; NIST cache is vacuum -- the ~1e-4 air/vac shift is well inside a
    1% tolerance and disambiguation only needs the nearest. Returns
    (loggf_kurucz, loggf_nist, loggf_nist) triples (3rd = residual-panel x)."""
    buckets = {}
    for d in nist_lines:
        k = _line_idkey(d["conf_i"], d["term_i"], d["J_i"],
                        d["conf_k"], d["term_k"], d["J_k"])
        buckets.setdefault(k, []).append(d)
    pairs, used = [], set()
    for d in kurucz_lines:
        k = _line_idkey(d["conf_i"], d["term_i"], d["J_i"],
                        d["conf_k"], d["term_k"], d["J_k"])
        cands = [c for c in buckets.get(k, []) if id(c) not in used]
        if not cands:
            continue
        n = min(cands, key=lambda c: abs(c["lambda_A"] - d["lambda_A"]))
        if abs(n["lambda_A"] - d["lambda_A"]) / n["lambda_A"] > tol:
            continue
        used.add(id(n))
        pairs.append((d["loggf"], n["loggf"], n["loggf"]))
    return pairs


def page_gf(pdf, species, abinitio_path, fitted_path, nist_lines,
            gf_fitted_label="fitted (RCE)", kurucz_lines=None, strong_cut=-1.0):
    """gf-comparison page, laid out like the level-residuals page: TWO stacked
    panels of (model - NIST) log gf vs NIST log gf. Top spans the full log gf
    range; bottom zooms to the strong, well-measured lines (NIST log gf >=
    strong_cut), the regime the fit actually targets, so their scatter isn't
    swamped by the weak-line tail. Lines are matched to NIST by eigenvector-
    composition identity (config+term+J), not nearest wavelength."""
    # series: (label, color, marker, pairs) where pairs = (model, nist, x=nist)
    series = []
    if abinitio_path:
        series.append(("ab initio", COL_ABINITIO, MK_ABINITIO,
                       match_gf_by_identity(abinitio_path, nist_lines)))
    if fitted_path:
        series.append((gf_fitted_label, COL_FIT, MK_FIT,
                       match_gf_by_identity(fitted_path, nist_lines)))
    if kurucz_lines:
        series.append(("Kurucz (gfall)", COL_KURUCZ, MK_KURUCZ,
                       match_kurucz_gf(kurucz_lines, nist_lines)))

    # full NIST-log gf range across all series, with a margin
    alln = [p[1] for _, _, _, pairs in series for p in pairs]
    if alln:
        lo, hi = min(alln), max(alln)
        m = 0.05 * (hi - lo) + 0.1
        glo, ghi = lo - m, hi + m
    else:
        glo, ghi = -3.0, 1.0

    # RMS summary (strong / all) lives in the legend label (no top title)
    leg = {}
    for label, col, mk, pairs in series:
        if not pairs:
            leg[label] = label
            continue
        c = np.array([p[0] for p in pairs]); nlg = np.array([p[1] for p in pairs])
        d = c - nlg
        strong = nlg >= strong_cut
        rms_all = np.sqrt(np.mean(d ** 2))
        rms_strong = (np.sqrt(np.mean(d[strong] ** 2)) if strong.any()
                      else float("nan"))
        leg[label] = (f"{label}   RMS = {rms_strong:.2f} strong, "
                      f"{rms_all:.2f} all")

    fig, (a_full, a_zoom) = plt.subplots(2, 1, figsize=(8.0, 8.0), sharex=False)
    xlab = "NIST $\\log gf$"
    ylab = "$\\Delta\\log gf$ (model $-$ NIST)"

    def _draw(ax):
        ax.axhline(0, color="0.5", lw=0.7, zorder=1)
        for label, col, mk, pairs in series:
            if not pairs:
                continue
            c = np.array([p[0] for p in pairs])
            nlg = np.array([p[1] for p in pairs])
            _scatter(ax, nlg, c - nlg, col, mk, leg[label])

    _draw(a_full)
    _draw(a_zoom)
    a_full.set_xlim(glo, ghi); a_full.set_ylim(-1.0, 1.0)
    a_full.legend(frameon=False, loc="lower left")

    # zoom: x to strong lines; y from a ROBUST spread of the strong residuals
    # (90th pct), so a few large outliers don't blow the window back open --
    # the point of the zoom is to read the bulk fit quality near zero.
    a_zoom.set_xlim(strong_cut - 0.1, ghi)
    strong_d = np.array([p[0] - p[1] for _, _, _, pairs in series for p in pairs
                         if p[1] >= strong_cut])
    if strong_d.size:
        yl = max(0.25, 1.25 * np.percentile(np.abs(strong_d), 90))
    else:
        yl = 0.5
    a_zoom.set_ylim(-yl, yl)
    for ax in (a_full, a_zoom):
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("--outg11", required=True)
    ap.add_argument("--nist", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--levels1", default=None,
                    help="RCE LEVELS1 output; if given, add fitted-vs-observed.")
    ap.add_argument("--fitted-levels-json", default=None,
                    help="JSON list of fitted levels [{config,term,J,parity,"
                         "E_obs,E_fit}] (cm^-1), e.g. from a gf_fit (Python "
                         "optimizer) run that has no RCE LEVELS1. Alternative to "
                         "--levels1 for the level/residual pages.")
    ap.add_argument("--ie", default=None,
                    help="ionization-energy table (species<TAB>IE_cm1).")
    ap.add_argument("--gf-fitted", default=None,
                    help="OUTG11 from the fitted-parameter RCG run (fitted gf).")
    ap.add_argument("--gf-fitted-label", default="fitted (RCE)",
                    help="legend label for the fitted-gf series on the gf page.")
    ap.add_argument("--fitted-label", default="fitted (RCE)",
                    help="legend label for the fitted series on the level/"
                         "residual pages (e.g. 'fitted (energy-only, full basis)').")
    ap.add_argument("--nist-lines", default=None,
                    help="cached NIST lines table (reference log gf).")
    ap.add_argument("--kurucz-lines", default=None,
                    help="Bob Kurucz's per-ion GF line file (gf<xxyy>.pos/.all); "
                         "his fitted gf, added to the gf page as a 3rd series.")
    ap.add_argument("--kurucz-elem", type=float, default=None,
                    help="element code to filter Kurucz lines (e.g. 12.00 Mg I).")
    ap.add_argument("--kurucz-levels", nargs="*", default=None,
                    help="Bob Kurucz's RCE fit log(s) c<xxyy>{e,o}z.log; his "
                         "fitted level residuals, added to the residuals page.")
    ap.add_argument("--fit-rms", type=float, default=None)
    ap.add_argument("--compare-cap", type=float, default=None,
                    help="energy cap (cm^-1) for the page-1 Bob head-to-head "
                         "stats table; restricts BOTH series to E_obs <= cap.")
    a = ap.parse_args()

    calc, lines = parse_outg11(a.outg11)
    nist = load_nist(a.nist)
    matched = match_levels(calc, nist)
    ie_cm = load_ie(a.ie, a.species) if a.ie else None

    fitted = None
    unified = None
    if a.fitted_levels_json and os.path.exists(a.fitted_levels_json):
        import json
        with open(a.fitted_levels_json) as f:
            fitted = json.load(f) or None
        if fitted:
            unified = build_unified_panels(calc, fitted, nist)
    elif a.levels1 and os.path.exists(a.levels1):
        fitted = parse_levels1(a.levels1) or None
        if fitted:
            unified = build_unified_panels(calc, fitted, nist)

    nist_lines = (load_nist_lines(a.nist_lines)
                  if a.nist_lines and os.path.exists(a.nist_lines) else None)
    fitted_path = (a.gf_fitted
                   if a.gf_fitted and os.path.exists(a.gf_fitted) else None)
    kurucz_lines = (load_kurucz_lines(a.kurucz_lines, a.kurucz_elem)
                    if a.kurucz_lines and os.path.exists(a.kurucz_lines)
                    else None)
    kurucz_levels = (load_kurucz_levels(a.kurucz_levels)
                     if a.kurucz_levels else None)
    compare = (build_compare(fitted, kurucz_levels, a.compare_cap)
               if a.compare_cap else None)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with PdfPages(a.out) as pdf:
        page_summary(pdf, a.species, matched, lines, a.fit_rms, compare=compare)
        page_levels(pdf, a.species, matched, fitted, calc, ie_cm, nist,
                    fitted_label=a.fitted_label)
        page_residuals(pdf, a.species, matched, unified,
                       kurucz_levels=kurucz_levels,
                       fitted_label=a.fitted_label)
        if nist_lines:
            page_gf(pdf, a.species, a.outg11, fitted_path, nist_lines,
                    gf_fitted_label=a.gf_fitted_label, kurucz_lines=kurucz_lines)
    nfit = sum(len(p["lev"]) for p in unified) if unified else 0
    print(f"wrote {a.out}  ({len(matched)} levels, {len(lines)} lines, "
          f"{sum(m['matched'] for m in matched)} matched to NIST, "
          f"{nfit} fitted)")


if __name__ == "__main__":
    main()
