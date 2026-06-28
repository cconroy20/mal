#!/usr/bin/env python3
"""
Generate figures for docs/cowan_explainer.tex.

These are physically correct, ILLUSTRATIVE figures (hydrogenic orbitals and the
exact Slater kernel) used to build intuition for what the Cowan codes compute.
They are not output of a specific RCN run; when we drive a real Sn/Fe calculation
we will add figures of the actual self-consistent P_nl(r).

Outputs PDFs into docs/figs/.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import genlaguerre, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({
    "font.size": 11, "axes.linewidth": 0.8, "figure.dpi": 150,
    "axes.spines.top": False, "axes.spines.right": False,
})


def P_nl(n, l, r, Z=1.0):
    """Hydrogenic reduced radial function P_nl(r) = r * R_nl(r), in atomic units.
    Normalized so that integral_0^inf P_nl^2 dr = 1."""
    rho = 2.0 * Z * r / n
    norm = np.sqrt((2.0 * Z / n) ** 3 * factorial(n - l - 1) /
                   (2.0 * n * factorial(n + l)))
    L = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    R = norm * np.exp(-rho / 2.0) * rho ** l * L
    return r * R


def fig_orbitals():
    """Radial wavefunctions P_nl(r) for several orbitals."""
    r = np.linspace(1e-4, 30, 2000)
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    orbitals = [(1, 0, "1s"), (2, 0, "2s"), (2, 1, "2p"),
                (3, 1, "3p"), (3, 2, "3d")]
    for n, l, lab in orbitals:
        ax.plot(r, P_nl(n, l, r), label=lab, lw=1.6)
    ax.axhline(0, color="k", lw=0.5, alpha=0.5)
    ax.set_xlabel(r"$r$  (atomic units, $a_0$)")
    ax.set_ylabel(r"$P_{n\ell}(r) = r\,R_{n\ell}(r)$")
    ax.set_xlim(0, 25)
    ax.legend(frameon=False, ncol=5, columnspacing=1.0,
              loc="upper center", bbox_to_anchor=(0.5, 1.12))
    fig.tight_layout()
    out = os.path.join(FIGS, "orbitals.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def fig_slater_kernel():
    """The Slater-integral kernel r_<^k / r_>^(k+1) as a function of r1 for fixed r2,
    plus the integrand P_3d^2(r1) * kernel that F^k actually integrates."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.6))

    r1 = np.linspace(1e-3, 12, 1500)
    r2 = 3.0
    for k, c in zip([0, 2, 4], ["C0", "C1", "C2"]):
        rlo = np.minimum(r1, r2)
        rhi = np.maximum(r1, r2)
        ker = rlo ** k / rhi ** (k + 1)
        a1.plot(r1, ker, color=c, lw=1.6, label=fr"$k={k}$")
    a1.axvline(r2, color="k", ls=":", lw=0.9)
    a1.text(r2 + 0.2, a1.get_ylim()[1] * 0.85, r"$r_2$", fontsize=11)
    a1.set_xlabel(r"$r_1$  ($a_0$)")
    a1.set_ylabel(r"$r_<^{\,k}\,/\,r_>^{\,k+1}$")
    a1.set_title(r"Coulomb kernel (fixed $r_2$)", fontsize=10)
    a1.legend(frameon=False)

    # integrand of F^k(3d,3d): P_3d^2(r1) * kernel * P_3d^2(r2), shown vs r1 at fixed r2
    r = np.linspace(1e-3, 20, 1500)
    P3d = P_nl(3, 2, r, Z=1.0)
    rho2 = 5.0
    rlo = np.minimum(r, rho2); rhi = np.maximum(r, rho2)
    for k, c in zip([0, 2, 4], ["C0", "C1", "C2"]):
        integ = P3d ** 2 * (rlo ** k / rhi ** (k + 1))
        a2.plot(r, integ, color=c, lw=1.6, label=fr"$k={k}$")
    a2.set_xlabel(r"$r_1$  ($a_0$)")
    a2.set_ylabel(r"$P_{3d}^2(r_1)\, r_<^{\,k}/r_>^{\,k+1}$")
    a2.set_title(r"$F^k$ integrand (slice at fixed $r_2$)", fontsize=10)
    a2.set_xlim(0, 20)
    a2.legend(frameon=False)

    fig.tight_layout()
    out = os.path.join(FIGS, "slater_kernel.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def fig_semiempirical():
    """Schematic: ab initio (HF) levels are systematically too spread out;
    scaling Slater integrals by ~0.85 (the semi-empirical fit) brings the
    computed levels into agreement with observed levels."""
    rng = np.random.default_rng(7)
    base = np.array([0.0, 1.0, 1.45, 2.7, 3.1, 3.4, 4.6])      # "observed" pattern
    obs = base
    # HF: too spread (Slater integrals ~15% too large -> term splitting too large)
    hf = base * 1.15 + 0.05 * rng.standard_normal(base.size)
    # fitted: scale back ~0.86 + small residual -> close to observed
    fit = base * 1.00 + 0.03 * rng.standard_normal(base.size)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    cols = {"HF (ab initio)": (0, "C3", hf),
            "Fitted (RCE)": (1, "C0", fit),
            "Observed (NIST)": (2, "k", obs)}
    for lab, (x, c, vals) in cols.items():
        for v in vals:
            ax.hlines(v, x - 0.32, x + 0.32, color=c, lw=2)
        ax.text(x, max(hf) + 0.35, lab, ha="center", fontsize=10, color=c)
    # connect fitted->observed to show agreement
    for vf, vo in zip(fit, obs):
        ax.plot([1.32, 1.68], [vf, vo], color="0.6", lw=0.6, ls="--")
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(-0.3, max(hf) + 0.7)
    ax.set_xticks([])
    ax.set_ylabel(r"Energy (arb.\ units, schematic)")
    ax.set_title("Semi-empirical fit: scale radial parameters to match observation",
                 fontsize=10)
    fig.tight_layout()
    out = os.path.join(FIGS, "semiempirical.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    fig_orbitals()
    fig_slater_kernel()
    fig_semiempirical()
    print("done.")
