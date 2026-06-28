#!/usr/bin/env python3
"""
Worked-example figures for docs/cowan_explainer.tex, built from a REAL run of the
Cowan chain on Sn7+ (configurations 4d^7 and 4d^6 5p), the shipped example.

Inputs (produced by running build/bin/{rcn,rcn2,rcg} in work/sn7plus/):
  - work/sn7plus/rwfn.dat : converged radial wavefunctions P_nl(r) (RCN, via the
                            physics-neutral dump patch; see notes/build_notes.md)
  - work/sn7plus/OUTG11   : energy levels, wavelengths, and log gf (RCG)

Outputs PDFs into docs/figs/:
  - sn7p_orbitals.pdf : real self-consistent P_nl(r) for the 4d^7 configuration
  - sn7p_spectrum.pdf : the computed electric-dipole spectrum (log gf vs wavelength)
"""
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.normpath(os.path.join(HERE, "..", "work", "sn7plus"))
FIGS = os.path.join(HERE, "figs")
os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({
    "font.size": 11, "axes.linewidth": 0.8, "figure.dpi": 150,
    "axes.spines.top": False, "axes.spines.right": False,
})


def read_rwfn(path):
    """Parse rwfn.dat -> list of dicts {conf,label,occ,r,P}."""
    orbitals = []
    cur = None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"#\s+(\S.*?)\s+label=\s*(\S+)\s+occ=\s*([\d.]+)", line)
                conf = line[2:].split("label=")[0].strip()
                lab = m.group(2)
                occ = float(m.group(3))
                cur = {"conf": conf, "label": lab, "occ": occ, "r": [], "P": []}
                orbitals.append(cur)
            elif line.strip():
                a, b = line.split()
                cur["r"].append(float(a)); cur["P"].append(float(b))
    for o in orbitals:
        o["r"] = np.array(o["r"]); o["P"] = np.array(o["P"])
    return orbitals


def fig_orbitals():
    orbs = read_rwfn(os.path.join(RUN, "rwfn.dat"))
    # take the 4d^7 configuration block (first 9 orbitals)
    sel = [o for o in orbs if o["conf"].endswith("4d7")]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for o in sel:
        ax.plot(o["r"], o["P"], lw=1.4, label=o["label"].strip())
    ax.axhline(0, color="k", lw=0.5, alpha=0.5)
    ax.set_xlabel(r"$r$  (atomic units, $a_0$)")
    ax.set_ylabel(r"$P_{n\ell}(r) = r\,R_{n\ell}(r)$")
    ax.set_xlim(0, 2.5)
    ax.set_title(r"Self-consistent radial wavefunctions, Sn$^{7+}$ $4d^7$ (RCN)",
                 fontsize=10, pad=10)
    ax.legend(frameon=False, ncol=3, columnspacing=0.9, fontsize=9,
              loc="lower right")
    fig.tight_layout()
    out = os.path.join(FIGS, "sn7p_orbitals.pdf")
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print("wrote", out, f"({len(sel)} orbitals)")


def read_outg11_spectrum(path):
    """Parse the ELEC DIP SPECTRUM table from OUTG11: rows with
    index, E, J, term, DELTA E, LAMBDA(A), LOG GF, GA, ... .
    Return arrays of wavelength (A) and log gf."""
    lam, loggf = [], []
    in_tab = False
    with open(path) as f:
        for line in f:
            if "LAMBDA(A)" in line and "LOG GF" in line:
                in_tab = True
                continue
            if not in_tab:
                continue
            # data rows: leading integer index, then floats; LOG GF is a signed
            # float typically in [-8, 1]. Be tolerant; require >=7 numeric tokens.
            toks = line.split()
            if len(toks) < 8:
                continue
            try:
                idx = int(toks[0])
            except ValueError:
                continue
            # find LAMBDA(A) and LOG GF by position: the documented columns are
            #   idx  E  J  <term...>  DELTA_E  LAMBDA  LOGGF  GA  CF/BR  group
            # term has variable width, so parse from the right: the last numeric
            # block is [..., LAMBDA, LOGGF, GA, CFBR, group]. GA is like 3.8E+09.
            nums = []
            for t in toks:
                try:
                    nums.append(float(t))
                except ValueError:
                    nums.append(None)
            # locate GA (scientific notation with E) -> LOGGF is two before it,
            # LAMBDA is three before it.
            gi = next((i for i, t in enumerate(toks)
                       if re.match(r"^[\d.]+E[+-]\d+$", t)), None)
            if gi is None or gi < 3:
                continue
            try:
                L = float(toks[gi - 2])
                G = float(toks[gi - 1])
            except ValueError:
                continue
            if 0.0 < L < 1e5 and -12 < G < 3:
                lam.append(L); loggf.append(G)
    return np.array(lam), np.array(loggf)


def fig_spectrum():
    lam, loggf = read_outg11_spectrum(os.path.join(RUN, "OUTG11"))
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    ax.scatter(lam, loggf, s=9, color="C0", alpha=0.6, edgecolors="none")
    ax.set_xlabel("Wavelength $\\lambda$  (Å)")
    ax.set_ylabel(r"$\log gf$")
    ax.set_title(r"Computed E1 spectrum of Sn$^{7+}$: "
                 r"$4d^7 \!-\! 4d^6 5p$ (RCG)", fontsize=10)
    ax.set_ylim(-9, 1)
    fig.tight_layout()
    out = os.path.join(FIGS, "sn7p_spectrum.pdf")
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print("wrote", out, f"({len(lam)} transitions)")


if __name__ == "__main__":
    fig_orbitals()
    fig_spectrum()
    print("done.")
