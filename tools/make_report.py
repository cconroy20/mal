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
from parse_cowan import parse_outg11  # noqa: E402

plt.rcParams.update({
    "font.size": 10, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
})


# ----------------------------------------------------------------------------
# data loading / matching
# ----------------------------------------------------------------------------
def load_nist(path):
    rows = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 5:
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


# orbital token with non-greedy occupation (same lookahead trick as parity)
_ORB = re.compile(r"[1-9][spdfghik]\d?(?![spdfghik])")


def _cfgkey(config):
    """Reduce a config to its valence orbital tokens (occupations normalized),
    e.g. 'Mg I   3s3p' -> '3s.3p', '3s.3p' -> '3s.3p', '2p6.3s2' -> '3s2'.
    Keeps the last two orbital tokens (the valence part); drops occupation '1'
    and closed-shell markers so computed and NIST forms agree."""
    toks = _ORB.findall(config)
    norm = []
    for t in toks:
        m = re.match(r"([1-9][spdfghik])(\d?)", t)
        nl, occ = m.group(1), m.group(2)
        norm.append(nl + (occ if occ and occ != "1" else ""))
    # a single doubly-occupied valence orbital (e.g. 3s2) is the ground-type
    # config; NIST often writes it with the closed core (2p6.3s2). Reduce to the
    # outermost token so the two forms agree.
    if norm and re.fullmatch(r"[1-9][spdfghik]2", norm[-1]):
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
def page_summary(pdf, species, matched, lines, fit_rms):
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


def page_levels(pdf, species, matched):
    """The centerpiece: a Grotrian-style term diagram (the standard convention).
    x-axis = term, grouped by parity (even block | odd block), each ordered by
    energy; y-axis = energy. J-levels are short horizontal ticks at their
    energies; computed solid, NIST observed dashed. Fine structure that is
    unresolved at this ordinate is reported precisely on the following table
    page. (A 'fitted' series will be added once RCE output exists.)"""
    panels = _group_terms(matched)
    # split into even / odd parity blocks, each ordered by energy
    even = [p for p in panels if p["par"] == "e"]
    odd = [p for p in panels if p["par"] == "o"]
    even.sort(key=_term_sort_key)
    odd.sort(key=_term_sort_key)

    # assign an x-slot to each term: even block on the left, odd on the right,
    # with a gap between blocks.
    layout = []   # (panel, x)
    x = 0
    for p in even:
        layout.append((p, x)); x += 1
    x += 1        # gap between parity blocks
    for p in odd:
        layout.append((p, x)); x += 1
    nx = max(1, x)

    fig, ax = plt.subplots(figsize=(max(7.0, 0.9 * nx + 2), 9.0))
    fig.suptitle(f"{species}: term (Grotrian) diagram", fontsize=14, y=0.97)

    half = 0.36
    for p, xc in layout:
        for m in p["lev"]:
            ec = m["E_calc"]; eo = m["E_obs"]
            ax.hlines(ec, xc - half, xc + half, color="C0", lw=1.8)
            if eo is not None:
                ax.hlines(eo, xc - half, xc + half, color="C3", lw=1.4, ls="--")

    # term labels under the x-axis (with superscripts)
    ax.set_xticks([xc for _, xc in layout])
    ax.set_xticklabels(
        [rf"{_term_mathlabel(p['tk'], p['par'])}" for p, _ in layout],
        fontsize=9)
    # second-line config labels
    ymin, ymax = ax.get_ylim()
    for p, xc in layout:
        ax.annotate(p["cfg"], xy=(xc, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, -28), textcoords="offset points",
                    ha="center", va="top", fontsize=7, color="0.4")

    # parity block headers
    if even:
        xe = np.mean([xc for p, xc in layout if p["par"] == "e"])
        ax.text(xe, 1.01, "even parity", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9, style="italic")
    if odd:
        xo = np.mean([xc for p, xc in layout if p["par"] == "o"])
        ax.text(xo, 1.01, "odd parity", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9, style="italic")

    ax.set_xlim(-0.8, nx - 0.2)
    ax.set_ylabel("Energy (cm$^{-1}$)")
    ax.spines["bottom"].set_position(("outward", 0))

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color="C0", lw=1.8, label="computed (ab initio)"),
               Line2D([0], [0], color="C3", lw=1.4, ls="--", label="observed (NIST)")]
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="upper left")
    fig.text(0.5, 0.01,
             "Standard term diagram: columns are terms (grouped by parity), "
             "ticks are J-levels. Solid = computed, dashed = observed. "
             "Per-level energies and residuals are tabulated on the next page.",
             ha="center", fontsize=8, style="italic")
    fig.tight_layout(rect=[0.04, 0.05, 1, 0.94])
    pdf.savefig(fig); plt.close(fig)

    # ---- companion fine-structure table page ----
    _page_level_table(pdf, species, panels)


def _page_level_table(pdf, species, panels):
    """Tabulate every level: config, term, J, E_calc, E_obs, delta E."""
    rows = []
    for p in panels:
        for m in sorted(p["lev"], key=lambda x: x["E_calc"]):
            try:
                Js = "%g" % float(m["J"])
            except (TypeError, ValueError):
                Js = str(m["J"])
            eo = m["E_obs"]
            term = _term_mathlabel(p["tk"], p["par"])
            rows.append([p["cfg"], term, Js,
                         f"{m['E_calc']:.1f}",
                         (f"{eo:.1f}" if eo is not None else "--"),
                         (f"{m['E_calc'] - eo:+.1f}" if eo is not None else "--")])
    # paginate at ~40 rows/page
    per = 40
    for start in range(0, max(1, len(rows)), per):
        chunk = rows[start:start + per]
        fig, ax = plt.subplots(figsize=(7.5, 9.5)); ax.axis("off")
        fig.suptitle(f"{species}: levels (computed vs observed)", fontsize=13,
                     y=0.97)
        col = ["config", "term", "J", "E_calc", "E_obs", "ΔE (cm⁻¹)"]
        tbl = ax.table(cellText=chunk, colLabels=col, loc="upper center",
                       cellLoc="center")
        tbl.auto_set_font_size(False); tbl.set_fontsize(8.5)
        tbl.scale(1, 1.3)
        pdf.savefig(fig); plt.close(fig)


def page_residuals(pdf, species, matched):
    res = [(m["E_obs"], m["E_calc"] - m["E_obs"], m["term"])
           for m in matched if m["matched"] and m["E_obs"] is not None]
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    if res:
        eo = np.array([r[0] for r in res]); dd = np.array([r[1] for r in res])
        ax.axhline(0, color="k", lw=0.6)
        ax.scatter(eo, dd, s=28, color="C0", zorder=3)
        rms = np.sqrt(np.mean(dd**2))
        ax.set_title(f"{species}: level-fit residuals   "
                     f"(RMS = {rms:.0f} cm$^{{-1}}$, ab initio)", fontsize=11)
    else:
        ax.text(0.5, 0.5, "no matched levels", ha="center")
    ax.set_xlabel("Observed level  $E_\\mathrm{obs}$  (cm$^{-1}$)")
    ax.set_ylabel("$E_\\mathrm{calc} - E_\\mathrm{obs}$  (cm$^{-1}$)")
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_gf(pdf, species, lines):
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    if lines:
        lam = np.array([d["lambda_A"] for d in lines])
        gf = np.array([d["loggf"] for d in lines])
        ax.scatter(lam, gf, s=10, color="C0", alpha=0.6, edgecolors="none")
        ax.set_title(f"{species}: computed E1 spectrum "
                     f"({len(lines)} lines)", fontsize=11)
    else:
        ax.text(0.5, 0.5, "no lines", ha="center")
    ax.set_xlabel("Wavelength $\\lambda$  (Å)")
    ax.set_ylabel("$\\log gf$")
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("--outg11", required=True)
    ap.add_argument("--nist", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fit-rms", type=float, default=None)
    a = ap.parse_args()

    calc, lines = parse_outg11(a.outg11)
    nist = load_nist(a.nist)
    matched = match_levels(calc, nist)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with PdfPages(a.out) as pdf:
        page_summary(pdf, a.species, matched, lines, a.fit_rms)
        page_levels(pdf, a.species, matched)
        page_residuals(pdf, a.species, matched)
        # gf-vs-wavelength page deferred: only meaningful once RCE produces many
        # lines and we have a reference (Kurucz/lab) gf to compare against.
    print(f"wrote {a.out}  ({len(matched)} levels, {len(lines)} lines, "
          f"{sum(m['matched'] for m in matched)} matched to NIST)")


if __name__ == "__main__":
    main()
