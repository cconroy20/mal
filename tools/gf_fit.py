#!/usr/bin/env python3
"""Tier-2 prototype: fit the Cowan radial parameters to BOTH observed energy
levels AND well-measured NIST gf values, using the existing RCG as a black-box
forward model wrapped in a Python optimizer.

Motivation: the standard RCE fit (energies only) was found to DEGRADE gf vs the
ab-initio RCG values (notes/mg1_gf_analysis.md), because fitting energies
reshapes eigenvectors in ways that hurt the dipole matrix elements. Here we put
gf into the objective so the fit can't trade away gf accuracy for a marginally
better energy.

Objective (minimized):
    chi2 = sum_i wE_i (E_calc_i - E_obs_i)^2            [levels, cm^-1]
         + lambda * sum_j wgf_j (loggf_calc_j - loggf_NIST_j)^2   [gf, dex]
  wE_i  : 1/sigma_E^2, sigma_E a flat energy tolerance (cm^-1)
  wgf_j : 1/sigma_gf^2 from NIST accuracy class; only log gf >= GF_MIN included
  lambda: relative weight of the gf term (tunable)

Forward model: write a trial parameter vector into ING11 via ing11_params, run
RCG (~8 ms for Mg I), parse levels + gf from OUTG11, match to NIST by
eigenvector-composition identity (make_report).

Usage:
    tools/gf_fit.py --run-dir work/mg1 --seed work/mg1/ING11.abinitio \
        --nist data/nist/MgI_levels.tsv --nist-lines data/nist/MgI_lines.tsv \
        --lambda 1.0 [--ridge 0.0] [--maxiter 4000]
"""
import argparse
import os
import subprocess
import sys

import numpy as np
from scipy.optimize import minimize, least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ing11_params as IP
import make_report as R
from parse_cowan import parse_outg11, parse_compositions, identify_lines


def _jstr(J):
    """Canonical J string matching parse_compositions keys ('%g' of float J)."""
    try:
        return "%g" % float(J)
    except (ValueError, TypeError):
        return None

# NIST accuracy class -> approx 1-sigma uncertainty in log gf (dex), ASD legend.
ACC_DEX = {"AAA": 0.013, "AA": 0.013, "A+": 0.013, "A": 0.022, "B+": 0.043,
           "B": 0.087, "C+": 0.13, "C": 0.22, "D+": 0.30, "D": 0.43, "E": 0.70}
GF_MIN = -1.0          # only fit gf for lines at least this strong
SIGMA_E = 50.0         # flat energy uncertainty (cm^-1) for the level term


def acc_sigma(a):
    return ACC_DEX.get(a.strip(), 0.5)


class Forward:
    """Black-box forward model: params -> RCG -> (levels, gf), with NIST targets
    precomputed so each evaluation just runs RCG and matches."""

    def __init__(self, run_dir, seed_ing11, nist_path, nist_lines_path):
        self.run_dir = os.path.abspath(run_dir)
        self.ing11 = os.path.join(self.run_dir, "ING11")
        self.outg11 = os.path.join(self.run_dir, "OUTG11")
        self.rcg = os.path.join(os.path.dirname(self.run_dir), "..",
                                "build", "bin", "rcg")
        self.rcg = os.path.abspath(self.rcg)
        self.raw, self.params = IP.parse(seed_ing11)
        self.seed = np.array([p["value"] for p in self.params])
        self.nist = R.load_nist(nist_path)
        self.nist_lines = R.load_nist_lines(nist_lines_path)
        # NIST level lookup by robust identity (cfgkey, termkey, Jkey)
        self.nist_lev = {}
        for nlv in self.nist:
            k = R._level_idkey(nlv["config"], nlv["term"], nlv["J"])
            self.nist_lev.setdefault(k, nlv["E_obs"])
        # ensure cowan.cfg points RCG at run_dir
        with open(os.path.join(self.run_dir, "cowan.cfg"), "w") as f:
            f.write(self.run_dir + "/\n")
        self.neval = 0
        self.setup_targets()      # pin residual structure from the seed

    def run(self, values):
        IP.write(self.raw, self.params, values, self.ing11)
        subprocess.run([self.rcg], cwd=self.run_dir, stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.neval += 1
        return None

    # ----- SLOT-based level access (stable under small parameter changes) -----
    # A level's "slot" = (parity, Jstr, energy-rank within that (parity,J) block).
    # Established once at the seed via eigenvector identity, then tracked by rank.
    # Ranks change only at true level crossings (rare, far from the optimum), so
    # the residual is smooth almost everywhere -> a clean finite-diff Jacobian.

    def _block_levels(self):
        """{(parity, Jstr) -> list of level dicts sorted by E_calc}."""
        comp = parse_compositions(self.outg11)
        return {bj: sorted(levs, key=lambda L: L["E_calc"])
                for bj, levs in comp.items()}

    def _level_at(self, blocks, slot):
        par, J, rank = slot
        levs = blocks.get((par, J))
        if levs is None or rank >= len(levs):
            return None
        return levs[rank]

    def _slot_of_energy(self, blocks, E, Jstr, parity=None):
        """Find the (parity,J,rank) slot whose level energy is closest to E,
        among blocks with this J (and this parity, if given). Constraining by
        parity is essential for line ends: the two ends of an E1 line have
        OPPOSITE parity, and without it a level in the other-parity block at a
        similar energy/J can be picked, collapsing two distinct lines onto one
        slot pair."""
        best = None
        for (par, J), levs in blocks.items():
            if J != Jstr or (parity is not None and par != parity):
                continue
            for rank, L in enumerate(levs):
                d = abs(L["E_calc"] - E)
                if best is None or d < best[0]:
                    best = (d, (par, J, rank))
        return best[1] if best else None

    @staticmethod
    def _parity_of(config):
        from parse_cowan import _parity_from_config
        return _parity_from_config(config)

    def energy_resid(self):
        """Vector of level residuals (E_calc - E_obs) over the fixed slots, with
        the global offset removed. Convenience for report() / metrics."""
        blocks = self._block_levels()
        out = []
        for slot, eobs in self.level_targets:
            L = self._level_at(blocks, slot)
            if L is not None:
                out.append(L["E_calc"] - eobs)
        out = np.array(out)
        if len(out):
            out = out - np.median(out)
        return out

    def gf_resid(self):
        """(residuals, sigmas) over the fixed line slots. Convenience for
        report() / metrics."""
        rows = self._current_line_rows()
        res, sig = [], []
        for slotpair, (nist, sigma) in self.line_targets.items():
            r = rows.get(slotpair)
            if r is not None:
                res.append(r - nist)
                sig.append(sigma)
        return np.array(res), np.array(sig)

    def _current_line_rows(self):
        """{(lower_slot, upper_slot) -> computed log gf} for the CURRENT OUTG11,
        attaching each computed line's ends to level slots by energy."""
        blocks = self._block_levels()
        out = {}
        for d in identify_lines(self.outg11):
            if not (d["config_low"] and d["config_up"]):
                continue
            lo = self._slot_of_energy(blocks, d["E_low"], _jstr(d["J_low"]),
                                      self._parity_of(d["config_low"]))
            up = self._slot_of_energy(blocks, d["E_up"], _jstr(d["J_up"]),
                                      self._parity_of(d["config_up"]))
            if lo and up:
                out.setdefault((lo, up), d["loggf"])
        return out

    def setup_targets(self):
        """Pin the residual structure from the seed. level_targets: list of
        (slot, E_obs). line_targets: {(lower_slot, upper_slot) -> (NIST loggf,
        sigma)} for strong matched lines. Both keyed by stable energy-rank slots
        so the residual vector keeps constant length/order across evaluations."""
        self.run(self.seed)
        blocks = self._block_levels()
        # level targets: each NIST-matched computed level -> its slot
        self.level_targets = []
        seen = set()
        for (par, J), levs in blocks.items():
            for rank, L in enumerate(levs):
                k = R._level_idkey(L["config"], L["term"], J)
                eobs = self.nist_lev.get(k)
                slot = (par, J, rank)
                if eobs is not None and k not in seen:
                    seen.add(k)
                    self.level_targets.append((slot, eobs))
        # line targets: strong matched lines -> {(lower_slot, upper_slot) ->
        # (NIST loggf, sigma)}. Each computed line (slot pair) is matched to its
        # NIST counterpart J-resolved, so multiplet components don't cross-match.
        self.line_targets = self._build_line_targets()

    def _build_line_targets(self):
        """For each computed line (keyed by its end slots), find the NIST line
        with the SAME level identities AND the same J pair, and store its log gf
        + sigma. J-resolved, so 3D J=1->3P J=1 and 3D J=2->3P J=2 stay distinct."""
        blocks = self._block_levels()
        # NIST lines indexed by (low-id, up-id, Jlo, Jup), order-independent
        nist_by_id = {}
        for nl in self.nist_lines:
            if nl["loggf"] < GF_MIN:
                continue
            a = (R._cfgkey(nl["conf_i"]), R._termkey(nl["term_i"]),
                 R._Jkey(nl["J_i"]))
            b = (R._cfgkey(nl["conf_k"]), R._termkey(nl["term_k"]),
                 R._Jkey(nl["J_k"]))
            nist_by_id[frozenset((a, b))] = (nl["loggf"], acc_sigma(nl["acc"]))
        out = {}
        for d in identify_lines(self.outg11):
            if not (d["config_low"] and d["config_up"]):
                continue
            a = (R._cfgkey(d["config_low"]), R._termkey(d["term_id_low"]),
                 R._Jkey(d["J_low"]))
            b = (R._cfgkey(d["config_up"]), R._termkey(d["term_id_up"]),
                 R._Jkey(d["J_up"]))
            hit = nist_by_id.get(frozenset((a, b)))
            if hit is None:
                continue
            lo = self._slot_of_energy(blocks, d["E_low"], _jstr(d["J_low"]),
                                      self._parity_of(d["config_low"]))
            up = self._slot_of_energy(blocks, d["E_up"], _jstr(d["J_up"]),
                                      self._parity_of(d["config_up"]))
            if lo and up:
                out[(lo, up)] = hit
        return out

    def resid_vector(self, x, scale, lam, ridge):
        """Stacked weighted residual vector for least_squares:
            [ (E_calc-E_obs)/SIGMA_E  per level slot,
              sqrt(lam)*(loggf-NIST)/sigma_gf  per line slot,
              sqrt(ridge)*(x-1)  (optional) ]
        Slots are tracked by energy rank, so the vector has constant length and
        smooth entries -> 0.5*||vec||^2 matches the scalar chi2."""
        self.run(x * scale)
        blocks = self._block_levels()
        line_rows = self._current_line_rows()
        # global energy offset from the level slots
        eres = []
        for slot, eobs in self.level_targets:
            L = self._level_at(blocks, slot)
            eres.append((L["E_calc"] - eobs) if L is not None else np.nan)
        eres = np.array(eres)
        good = ~np.isnan(eres)
        off = np.median(eres[good]) if good.any() else 0.0
        vec = []
        for v, ok in zip(eres, good):
            vec.append(((v - off) / SIGMA_E) if ok else 0.0)
        sl = np.sqrt(lam)
        for slotpair, (nist, sigma) in self.line_targets.items():
            r = line_rows.get(slotpair)
            vec.append(sl * (r - nist) / sigma if r is not None else 0.0)
        if ridge:
            vec.extend(np.sqrt(ridge) * (np.asarray(x) - 1.0))
        return np.asarray(vec)


def _scale_of(fwd):
    return np.where(np.abs(fwd.seed) > 1e-6, np.abs(fwd.seed), 1.0)


def make_objective(fwd, lam, ridge):
    """Scalar chi2 objective (Nelder-Mead path; also used by the lambda sweep)."""
    scale = _scale_of(fwd)

    def obj(x):
        values = x * scale          # x is in seed-relative units
        fwd.run(values)
        eres = fwd.energy_resid()
        chi2 = np.sum((eres / SIGMA_E) ** 2) if len(eres) else 0.0
        gres, gsig = fwd.gf_resid()
        if len(gres):
            chi2 += lam * np.sum((gres / gsig) ** 2)
        if ridge:
            chi2 += ridge * np.sum((x - 1.0) ** 2)   # pull toward seed
        return chi2

    return obj, scale


def fit_lm(fwd, lam, ridge, maxiter, presolve=True):
    """Levenberg-Marquardt fit on the stacked residual vector. Returns
    (best_values, result, scale, cov) where cov is the parameter covariance in
    SCALED units (x = value/scale), estimated from the Jacobian at the optimum
    (the basis for gf MC uncertainties). 0.5*||resid||^2 == the scalar chi2, so
    the optimum and lambda semantics match the Nelder-Mead path.

    LM is a LOCAL gradient method: from the bare RCE seed it under-converges
    (stops in a poorer basin) vs Nelder-Mead. So by default do a short NM
    PRE-SOLVE to land in the right basin, then LM to polish and -- critically --
    to produce the Jacobian/covariance. presolve=False uses LM alone (e.g. when
    the seed is already near the optimum)."""
    scale = _scale_of(fwd)
    x0 = np.ones_like(fwd.seed)
    if presolve:
        obj, _ = make_objective(fwd, lam, ridge)
        nm = minimize(obj, x0, method="Nelder-Mead",
                      options={"maxiter": maxiter, "maxfev": maxiter,
                               "xatol": 1e-4, "fatol": 1e-3, "adaptive": True})
        x0 = nm.x
    res = least_squares(
        lambda x: fwd.resid_vector(x, scale, lam, ridge), x0,
        method="lm", max_nfev=maxiter, xtol=1e-8, ftol=1e-8,
        diff_step=1e-3)              # finite-diff step in scaled units
    best = res.x * scale
    # --- covariance & its ERROR-MODEL assumptions (read before trusting sigmas) -
    # cov = (J^T J)^-1 * s^2, with s^2 = 2*cost/dof the EFFECTIVE per-residual
    # variance estimated FROM THE RESIDUAL SCATTER (reduced chi^2), NOT from
    # propagated measurement errors. This is the standard "let the fit's scatter
    # define the noise" scheme, appropriate here because the residuals are
    # dominated by MODEL inadequacy (missing configs/CI), not measurement noise.
    # The per-residual WEIGHTS that shape cov are: gf -> NIST accuracy class
    # (real); energy -> a flat SIGMA_E (a relative weight, NOT a true ~0.01 cm^-1
    # level error); ridge -> the prior strength. Consequences:
    #   * sigmas are CONDITIONAL ON THE MODEL; they do not include model error.
    #   * null-space (ridge-only) params get sigma ~ 1/sqrt(ridge): that is the
    #     PRIOR width, not data-derived. self.cov_reff records s for honesty.
    J = res.jac
    cov = None
    try:
        JTJ = J.T @ J
        cov = np.linalg.inv(JTJ)
        dof = max(len(res.fun) - len(res.x), 1)
        s2 = 2.0 * res.cost / dof            # res.cost = 0.5*||r||^2
        cov *= s2
        fwd.cov_reff = np.sqrt(s2)           # effective noise scale (reduced chi)
        fwd.cov_dof = dof
    except np.linalg.LinAlgError:
        pass
    return best, res, scale, cov


def report(fwd, values, label):
    fwd.run(values)
    eres = fwd.energy_resid()
    gres, gsig = fwd.gf_resid()
    erms = np.sqrt(np.mean(eres ** 2)) if len(eres) else float("nan")
    grms = np.sqrt(np.mean(gres ** 2)) if len(gres) else float("nan")
    wgrms = (np.sqrt(np.sum((gres / gsig) ** 2) / np.sum(1 / gsig ** 2))
             if len(gres) else float("nan"))
    print(f"{label:>16}: level RMS = {erms:7.1f} cm^-1   "
          f"strong-gf RMS = {grms:.3f}   acc-wtd gf RMS = {wgrms:.3f}   "
          f"(N_E={len(eres)}, N_gf={len(gres)})")
    return erms, grms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--seed", required=True, help="ING11 to seed from.")
    ap.add_argument("--nist", required=True)
    ap.add_argument("--nist-lines", required=True)
    ap.add_argument("--lambda", dest="lam", type=float, default=1.0)
    ap.add_argument("--ridge", type=float, default=0.0)
    ap.add_argument("--maxiter", type=int, default=4000)
    ap.add_argument("--method", choices=["lm", "nelder-mead"], default="lm",
                    help="lm = Levenberg-Marquardt on the residual vector "
                         "(fast, gives the Jacobian/covariance); nelder-mead = "
                         "the old derivative-free scalar minimizer.")
    ap.add_argument("--out", default=None, help="write fitted ING11 here.")
    a = ap.parse_args()

    fwd = Forward(a.run_dir, a.seed, a.nist, a.nist_lines)
    print(f"{len(fwd.params)} adjustable params; seeded from {a.seed}")
    report(fwd, fwd.seed, "seed")

    if a.method == "lm":
        best, res, scale, cov = fit_lm(fwd, a.lam, a.ridge, a.maxiter)
        print(f"\noptimizer (LM): {res.message}  "
              f"(nfev={fwd.neval}, cost={res.cost:.3f}, "
              f"||r||^2={2 * res.cost:.3f})")
        report(fwd, best, "gf-fit")
        if cov is not None:
            sig_scaled = np.sqrt(np.clip(np.diag(cov), 0, None))
            sig = sig_scaled * scale                       # value units
            reff = getattr(fwd, "cov_reff", float("nan"))
            # prior width per param (scaled units): ridge term is sqrt(ridge)*(x-1)
            # so a wholly-unconstrained param has scaled sigma ~ reff/sqrt(ridge).
            prior_w = (reff / np.sqrt(a.ridge)) if a.ridge else float("inf")
            print("\n  ERROR MODEL: sigmas are CONDITIONAL ON THE MODEL. The "
                  "per-residual\n  noise is the EFFECTIVE scatter (reduced "
                  f"chi = {reff:.2f}, dof={getattr(fwd,'cov_dof','?')}), not "
                  "propagated\n  measurement error; gf weights = NIST accuracy, "
                  "energy weight = flat\n  SIGMA_E (relative only), and ridge "
                  f"sets a prior of width ~{prior_w:.2g} (scaled).\n  Params at "
                  ">=0.7*prior width are PRIOR-DOMINATED (data uninformative).")
            print("\n  parameter 1-sigma (from Jacobian at optimum):")
            for p, b, s, ss in zip(fwd.params, best, sig, sig_scaled):
                flag = " <prior-dominated>" if (a.ridge
                        and ss >= 0.7 * prior_w) else ""
                print(f"    {p['key']:18} {b:12.5f} +/- {s:11.5f}{flag}")
    else:
        obj, scale = make_objective(fwd, a.lam, a.ridge)
        x0 = np.ones_like(fwd.seed)
        res = minimize(obj, x0, method="Nelder-Mead",
                       options={"maxiter": a.maxiter, "maxfev": a.maxiter,
                                "xatol": 1e-4, "fatol": 1e-3, "adaptive": True})
        best = res.x * scale
        print(f"\noptimizer (NM): {res.message}  "
              f"(nfev={fwd.neval}, chi2={res.fun:.3f})")
        report(fwd, best, "gf-fit")

    if a.out:
        IP.write(fwd.raw, fwd.params, best, a.out)
        print(f"wrote fitted parameters -> {a.out}")


if __name__ == "__main__":
    main()
