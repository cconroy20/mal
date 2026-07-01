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

Regularization: --ridge applies a per-parameter Gaussian PRIOR pulling each
radial parameter toward its ab-initio Hartree-Fock value (--abinitio), with a
per-kind width (PRIOR_SIGMA: loose EAV, tighter Slater, tightest CI). This is a
smooth, general form of Bob Kurucz's FIXEDHF discipline -- the data overrides the
prior where it has information, and the prior keeps the rest near HF, making the
otherwise rank-deficient fit well-posed (full-rank Jacobian -> usable covariance).
ridge~1 keeps Mg I gf quality (~0.08 dex) while curing the rank deficiency.

Usage:
    tools/gf_fit.py --run-dir work/mg1 --seed work/mg1/ING11.fit \
        --abinitio work/mg1/ING11.abinitio \
        --nist data/nist/MgI_levels.tsv --nist-lines data/nist/MgI_lines.tsv \
        --lambda 3 --ridge 1.0 [--free-kinds EAV,P] [--maxiter 6000]
"""
import argparse
import os
import re
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

# Per-kind ridge-PRIOR width (Bob-style discipline as a smooth prior): how far,
# in RELATIVE terms, the data is allowed to pull each parameter away from its
# ab-initio Hartree-Fock value. Loose on EAVs (centroids -- data should move
# them freely); tighter on single-config Slater F/G/zeta ('P'); tightest on CI
# interaction integrals (least trusted to fit, à la Bob's FIXEDHF for CI). These
# are fractional sigmas; sigma -> inf reproduces an unregularized free fit, small
# sigma approaches a hard freeze. Scaled by the global --ridge knob.
PRIOR_SIGMA = {"EAV": 1.0, "P": 0.30, "CI": 0.15}

# NIST accuracy class -> approx 1-sigma uncertainty in log gf (dex), ASD legend.
ACC_DEX = {"AAA": 0.013, "AA": 0.013, "A+": 0.013, "A": 0.022, "B+": 0.043,
           "B": 0.087, "C+": 0.13, "C": 0.22, "D+": 0.30, "D": 0.43, "E": 0.70}
GF_MIN = -1.0          # only fit gf for lines at least this strong
SIGMA_E = 50.0         # flat energy uncertainty (cm^-1) for the level term

# Above-cap levels are admitted as fit targets ONLY for these explicitly
# whitelisted valence-perturber configs (cfgkeys). 3p^2 is the one demonstrated
# win: its 1S@68275 fits to ~-13 cm^-1 and pins F^2(3p^2). The other above-cap
# valence configs the parser fix exposed (3p3d, 3p4s, 3p4d) look tempting but are
# DESTROYED when fit (under-constrained -> LM drags their EAV the wrong way; see
# setup_targets); they are reported from the ab-initio model, not fit. Extend this
# list per-ion only for a valence level shown to fit STABLY. (Kept as cfgkeys so a
# future 3d^2 perturber, if it fits cleanly, is a one-line addition.)
CAP_EXEMPT_CONFIGS = {"3p2"}

# Only RELEASE the ridge on zeta params with |zeta| >= this (cm^-1): the low-n,
# observable-fine-structure configs (3s3p ~28, 3p2 ~31). High-n Rydberg zetas
# (~1 cm^-1) are unconstrained noise that ill-conditions the fit if released.
ZETA_MIN_CM = 5.0


def _max_orbital_n(cfgkey):
    """Largest principal quantum number among a cfgkey's orbitals, e.g.
    '3p2'->3, '3s.10d'->10, '3d.4s'->4. cfgkey orbitals are <n><l>[occ] tokens
    joined by '.'. Used to tell a low-n valence config from a Rydberg member."""
    ns = [int(m) for m in re.findall(r"(\d+)[spdfghik]", cfgkey)]
    return max(ns) if ns else 0


def acc_sigma(a):
    return ACC_DEX.get(a.strip(), 0.5)


class Forward:
    """Black-box forward model: params -> RCG -> (levels, gf), with NIST targets
    precomputed so each evaluation just runs RCG and matches."""

    def __init__(self, run_dir, seed_ing11, nist_path, nist_lines_path,
                 free_kinds=None, abinitio_ing11=None, obs_configs=None,
                 free_ci_pairs=None, free_ci_sigma=None, max_energy=None,
                 outgine=None, zeta_sigma=None, zeta_configs=None,
                 fit_gf_scales=False, free_p_keys=None,
                 eav_sigma=None, eav_configs=None):
        # remember the construction kwargs so parallel-Jacobian worker PROCESSES
        # can rebuild an identical Forward in their own run_dir (run_dir overridden).
        self.init_kwargs = dict(
            seed_ing11=seed_ing11, nist_path=nist_path,
            nist_lines_path=nist_lines_path, free_kinds=free_kinds,
            abinitio_ing11=abinitio_ing11, obs_configs=obs_configs,
            free_ci_pairs=free_ci_pairs, free_ci_sigma=free_ci_sigma,
            max_energy=max_energy, outgine=outgine, zeta_sigma=zeta_sigma,
            zeta_configs=zeta_configs, fit_gf_scales=fit_gf_scales,
            free_p_keys=free_p_keys, eav_sigma=eav_sigma,
            eav_configs=eav_configs)
        self.run_dir = os.path.abspath(run_dir)
        # Optionally cap the fitted LEVELS at this observed energy (cm^-1), like
        # build_ine20 --max-energy: high-Rydberg / near-IE levels are unreliable
        # on an incomplete basis and their rank-based slots collide (degenerate
        # near-IE levels) -> bogus ~2000 cm^-1 residuals that would dominate a
        # pure-energy (gf-off) objective. None = no cap (every NIST-matched level).
        self.max_energy = max_energy
        self.ing11 = os.path.join(self.run_dir, "ING11")
        self.outg11 = os.path.join(self.run_dir, "OUTG11")
        self.rcg = os.path.join(os.path.dirname(self.run_dir), "..",
                                "build", "bin", "rcg")
        self.rcg = os.path.abspath(self.rcg)
        self.raw, allp = IP.parse(seed_ing11)
        # SELECTIVELY-FREED CI integrals: a set of frozensets of two config-name
        # strings (e.g. frozenset({"3s2","3p2"})). These specific interaction
        # integrals are freed regardless of free_kinds/obs_configs -- they couple
        # the singlet-system configs (3s2/3p2/3d2/3s.nd) whose FROZEN raw-HF CI
        # otherwise mis-sets the 1D/1P eigenvector mixing, costing ~0.5 dex on the
        # singlet gf (settled diagnosis). Freeing only the few DOMINANT low-l CI
        # avoids the divergence of freeing all correlated CI at once. CI keys carry
        # the readable pair, e.g. "3s2     - 3p2|CI0".
        free_ci_pairs = free_ci_pairs or set()

        def _ci_pair(p):
            if p["kind"] != "CI":
                return None
            a, b = (x.strip() for x in p["key"].split("|")[0].split("-"))
            return frozenset((a, b))

        # Bob-style discipline: optionally restrict the FREE set to certain
        # parameter kinds (e.g. {"EAV"} or {"EAV","P"}); the rest stay pinned at
        # their ab-initio seed value (IP.write only touches self.params). 'P' =
        # single-config Slater/zeta (F^k,G^k,zeta); 'CI' = interaction integrals.
        # HYBRID: free_p_keys force specific single-config Slater integrals (by
        # exact ING11 key) into the free set even when their kind ('P') is NOT in
        # free_kinds -- e.g. free the 3s.np G^1 exchange series individually while
        # the rest of the structure stays frozen+scaled (Bob's recipe). These same
        # keys are EXCLUDED from the G/F global scaling below (no double-control).
        free_p_keys = set(free_p_keys or ())
        if free_kinds:
            allp = [p for p in allp if p["kind"] in free_kinds
                    or _ci_pair(p) in free_ci_pairs
                    or p["key"] in free_p_keys]
        # ...and optionally drop params whose CONFIG has no observed level: the
        # data can't constrain those, so freezing them at HF cuts the free-param
        # count (Jacobian cost) with no loss -- on the full Mg I basis this is
        # 269 -> 170. obs_configs is a set of _cfgkey strings. Hand-picked free CI
        # pairs are exempt (they couple configs we deliberately want to fit).
        if obs_configs is not None:
            allp = [p for p in allp
                    if _ci_pair(p) in free_ci_pairs
                    or R._cfgkey(p["key"].split("|")[0].replace("Mg I", "")
                                 .strip()) in obs_configs]
        self.params = allp
        self.free_ci_pairs = free_ci_pairs
        self.seed = np.array([p["value"] for p in self.params])
        # Ridge-prior CENTRES = ab-initio (HF) values per free param, looked up
        # by key from the ab-initio ING11; fall back to the seed value if absent.
        # Per-kind PRIOR_SIGMA sets how far data may pull each from its centre.
        abv = IP.values_by_key(abinitio_ing11) if abinitio_ing11 else {}
        self.prior_centre = np.array([abv.get(p["key"], p["value"])
                                      for p in self.params])
        self.prior_sigma = np.array([PRIOR_SIGMA.get(p["kind"], 1.0)
                                     for p in self.params])
        # Optionally RELAX the ridge prior on the hand-picked free CI integrals:
        # the default CI sigma (0.15, tightest) is meant to PIN the ~5500 frozen
        # CI near scaled-HF, but the few CI we deliberately freed to reshape the
        # 1D/1P mixing should be allowed to move as far as the 9-config fit moved
        # them (~60%). free_ci_sigma overrides PRIOR_SIGMA['CI'] for ONLY the
        # free_ci_pairs params (None = leave at the default tight prior).
        if free_ci_sigma is not None and free_ci_pairs:
            for i, p in enumerate(self.params):
                if _ci_pair(p) in free_ci_pairs:
                    self.prior_sigma[i] = free_ci_sigma
        # ZETA POLISH: loosen the ridge on the spin-orbit (zeta) params so the
        # OBSERVED fine structure can pull them to the right value. zeta's HF (our
        # code) sits well below Bob's, so a screening scale can't reach the truth;
        # but the fine-structure data CAN (zeta=41 reproduces the observed 3s3p 3P
        # splitting). With the default sigma_P=0.30 and a centre below truth, the
        # ridge under-shoots (the move is >1 sigma vs a sub-SIGMA_E energy gain).
        # zeta_sigma raises the per-zeta sigma so that move is cheap.
        #
        # SURGICAL: release zeta ONLY for the named zeta_configs -- the few configs
        # with OBSERVED, resolved fine structure that actually constrains zeta
        # (e.g. {'3s.3p','3p2'}). Releasing zeta broadly (all 14, or all |zeta|>5
        # which pulls in the 3p.ns core-zeta series) adds near-unconstrained
        # directions that ILL-CONDITION the Jacobian -> LM thrashes (chi2 oscillated
        # 200-360, spiked to 3.6e5). zeta_configs are cfgkeys; default None = none.
        if zeta_sigma is not None and outgine is not None and zeta_configs:
            import param_labels as PL
            want = {R._cfgkey(c) for c in zeta_configs}
            zeta_keys = {p["key"]
                         for p in PL.physical_params(seed_ing11, outgine,
                                                     strict=False)
                         if p.get("kind") == "ZETA" and p.get("cfg") in want}
            for i, p in enumerate(self.params):
                if p["key"] in zeta_keys:
                    self.prior_sigma[i] = zeta_sigma
        # SURGICAL EAV RELEASE: loosen the ridge on the CENTROID (EAV) of named
        # configs so the data can pull it far from the scaled-HF centre. The default
        # EAV sigma (1.0 kK) pins the centroid near HF -- fine for low valence, but
        # for the doubly-excited configs the parser fix newly exposed (3p3d, 3p4s,
        # 3p4d) the scaled-HF centroid is ~7000-20000 cm^-1 (7-20 sigma) off, so the
        # ridge fights the fit and the level stays badly placed. Bob fits these
        # centroids freely (-> his 3p3d residuals <=3 cm^-1); this is the analogue.
        # SURGICAL by config (eav_configs): a broad EAV release re-loosens the whole
        # spectrum and de-anchors the well-placed bulk. eav_configs are cfgkeys.
        if eav_sigma is not None and eav_configs:
            want = {R._cfgkey(c) for c in eav_configs}
            for i, p in enumerate(self.params):
                if p["kind"] != "EAV":
                    continue
                cfg = R._cfgkey(p["key"].split("|")[0].replace("Mg I", "").strip())
                if cfg in want:
                    self.prior_sigma[i] = eav_sigma
        # TWO FITTED SCREENING SCALES (s_G, s_F): two GLOBAL multipliers that the
        # fit adjusts, applied at write-time to ALL FROZEN G (exchange) and F
        # (direct) single-config integrals -- relative to their AB-INITIO value.
        # This is Bob's two-factor screening (his 101 scaled FIXEDHF integrals use
        # G~0.6, F~0.8 on HIS HF), but FIT to our deck (our HF differs, so our
        # optimum is ~0.8/0.8 -- the scale isn't transferable, only the physical
        # target). 2 params capture the systematic HF screening with a strong,
        # well-conditioned signal (term splittings the EAV fit can't touch),
        # instead of dozens of ill-conditioned per-integral freedoms.
        self.fit_gf_scales = fit_gf_scales
        self._scale_slots = []      # (lineno, col, hf_value, group, 'G'|'F')
        self._n_scale = 0
        if fit_gf_scales and outgine is not None:
            import param_labels as PL
            # map each single-config Slater's ING11 location -> physical kind, via
            # the ab-initio deck (so the scale multiplies the AB-INITIO value).
            ab_raw, ab_params = IP.parse(abinitio_ing11)
            loc = {p["key"]: p for p in ab_params}      # key -> abinitio param (has lineno,col,value,group)
            for p in PL.physical_params(abinitio_ing11, outgine, strict=False):
                if (p.get("kind") in ("G", "F") and p["key"] in loc
                        and p["key"] not in free_p_keys):   # don't scale freed ones
                    ap = loc[p["key"]]
                    self._scale_slots.append(
                        (ap["lineno"], ap["col"], ap["value"], ap["group"],
                         p["kind"]))
            self._n_scale = 2       # [s_G, s_F] appended to the param vector
            # extend params/seed/prior with the two scale params (seed at 0.8)
            self.params = self.params + [
                {"key": "GLOBAL|sG", "kind": "SCALE", "value": 0.8},
                {"key": "GLOBAL|sF", "kind": "SCALE", "value": 0.8}]
            self.seed = np.append(self.seed, [0.8, 0.8])
            self.prior_centre = np.append(self.prior_centre, [0.8, 0.8])
            self.prior_sigma = np.append(self.prior_sigma, [0.5, 0.5])
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

    def clone_for_worker(self, worker_dir):
        """A lightweight copy that SHARES this Forward's params, targets, prior and
        scale (so its residual vector aligns column-for-column) but runs RCG in its
        OWN worker_dir (isolated ING11/OUTG11). Used by the parallel Jacobian: each
        worker perturbs one parameter and runs RCG concurrently without clobbering
        the master's or each other's files. Does NOT re-run setup_targets (targets
        are shared from the master, pinned once at the seed)."""
        import copy as _copy
        w = _copy.copy(self)                      # shallow: shares lists/dicts
        w.run_dir = os.path.abspath(worker_dir)
        w.ing11 = os.path.join(w.run_dir, "ING11")
        w.outg11 = os.path.join(w.run_dir, "OUTG11")
        with open(os.path.join(w.run_dir, "cowan.cfg"), "w") as f:
            f.write(w.run_dir + "/\n")
        w.neval = 0
        return w

    def run(self, values):
        values = np.asarray(values)
        raw = self.raw
        if self._n_scale:
            # last two entries are [s_G, s_F]; the rest map to ING11 slots
            s_G, s_F = float(values[-2]), float(values[-1])
            ing_params = self.params[:-self._n_scale]
            ing_values = values[:-self._n_scale]
            # apply the fitted scales to the FROZEN G/F slots (relative to ab-initio)
            raw = list(self.raw)
            by_line = {}
            for lineno, col, hf, grp, kind in self._scale_slots:
                by_line.setdefault(lineno, []).append((col, hf, grp, kind))
            for lineno, edits in by_line.items():
                s = raw[lineno].rstrip("\n")
                for col, hf, grp, kind in edits:
                    scaled = hf * (s_G if kind == "G" else s_F)
                    field = IP._fmt_p(scaled, grp)         # 9-digit value + group
                    s = s[:col] + field + s[col + 10:]
                raw[lineno] = s + "\n"
        else:
            ing_params, ing_values = self.params, values
        IP.write(raw, ing_params, ing_values, self.ing11)
        subprocess.run([self.rcg], cwd=self.run_dir, stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.neval += 1
        return None

    # ----- IDENTITY-based level access (RCE-style, re-resolved each eval) -----
    # Each computed level carries a ROBUST identity = its dominant eigenvector
    # basis state (config, term) + its block J. We match a target to the current
    # computed level of the SAME (cfgkey, termkey, Jstr) identity EVERY evaluation,
    # exactly as Cowan's RCE does internally -- rather than tracking a frozen
    # energy-RANK slot. Rank tracking manufactured ~2000 cm^-1 residuals when
    # near-degenerate high-n levels swapped order, and silently dropped levels it
    # couldn't track; identity matching is immune to both (a level keeps its
    # identity through crossings). The residual stays aligned by identity ->
    # constant length/order -> a clean finite-diff Jacobian, same as before.

    def _block_levels(self):
        """{(parity, Jstr) -> list of level dicts sorted by E_calc}."""
        comp = parse_compositions(self.outg11)
        return {bj: sorted(levs, key=lambda L: L["E_calc"])
                for bj, levs in comp.items()}

    def _levels_by_identity(self):
        """{(cfgkey, termkey, Jstr) -> level dict} for the CURRENT OUTG11, using
        each level's dominant-eigenvector (config, term) identity and its block J.
        If two computed levels share an identity (rare near-degeneracy), keep the
        one nearest the target observed energy when known, else the lowest-energy
        one -- deterministic either way."""
        comp = parse_compositions(self.outg11)
        by_id = {}
        for (par, J), levs in comp.items():
            for L in levs:
                k = R._level_idkey(L["config"], L["term"], J)
                if k[0] is None or k[1] is None or k[2] is None:
                    continue
                prev = by_id.get(k)
                if prev is None:
                    by_id[k] = L
                else:
                    # collision: prefer the one closer to this identity's E_obs
                    eobs = getattr(self, "_target_eobs", {}).get(k)
                    if eobs is not None:
                        if abs(L["E_calc"] - eobs) < abs(prev["E_calc"] - eobs):
                            by_id[k] = L
                    elif L["E_calc"] < prev["E_calc"]:
                        by_id[k] = L
        return by_id

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
        """Vector of level residuals (E_calc - E_obs) over the fixed IDENTITY
        targets, with the global offset removed. Convenience for report()/metrics.
        A target with no current same-identity computed level contributes nothing
        (it should be rare with identity matching, unlike rank slots)."""
        by_id = self._levels_by_identity()
        out = []
        for k, eobs in self.level_targets:
            L = by_id.get(k)
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
        (identity, E_obs) where identity=(cfgkey,termkey,Jstr) -- re-resolved to a
        computed level by identity each evaluation (RCE-style), NOT a frozen rank.
        line_targets: {(lower_slot, upper_slot) -> (NIST loggf, sigma)} for strong
        matched lines. The residual vector keeps constant length/order (one entry
        per target identity) across evaluations."""
        self.run(self.seed)
        # level targets: every computed level at the seed whose identity matches a
        # NIST level (below the energy cap). Identity matching means we no longer
        # lose targets to rank swaps, so this set is >= the old rank-slot set.
        #
        # CAP EXEMPTION: the max_energy cap drops unreliable HIGH-RYDBERG /
        # near-degenerate series members (whose rank-degenerate slots gave the old
        # ~14000 cm^-1 artifacts). ONE above-cap config is worth exempting: 3p^2. Its
        # 1S@68275 fits to ~-13 cm^-1 and pins the 1S-3P splitting -> F^2(3p^2) ->
        # the position of the (unobserved) 3p^2 1D perturber that mixes into the
        # 3s.nd 1D series. That level carries irreplaceable, WELL-CONSTRAINED signal.
        #
        # The OTHER above-cap valence configs the parser fix exposed (3p3d, 3p4s,
        # 3p4d) are the OPPOSITE: at ab-initio they sit ~1500 cm^-1 from observed
        # (good!), but as fit TARGETS they are DESTROYED -- weakly weighted (8 levels
        # vs 100) and under-constrained, the LM drags their EAV 8000-20000 the wrong
        # way (verified: mean|resid| 1477 ab-initio -> 10454 fitted; an EAV-ridge
        # release made it worse, not better -- these are flat, low-gradient
        # directions, exactly the unobserved-perturber pathology). So they must NOT
        # be fit targets. We report them from the (already-good) ab-initio model
        # instead. RULE: exempt an above-cap level only if its config is in the
        # explicit CAP_EXEMPT_CONFIGS whitelist of reliable valence perturbers.
        self.level_targets = []
        self._target_eobs = {}
        seen = set()
        comp = parse_compositions(self.outg11)
        for (par, J), levs in comp.items():
            for L in levs:
                k = R._level_idkey(L["config"], L["term"], J)
                eobs = self.nist_lev.get(k)
                if eobs is not None and k not in seen:
                    if (self.max_energy is not None and eobs > self.max_energy
                            and k[0] not in CAP_EXEMPT_CONFIGS):
                        continue          # above cap, not a whitelisted perturber
                    seen.add(k)
                    self.level_targets.append((k, eobs))
                    self._target_eobs[k] = eobs
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
            [ (E_calc-E_obs)/SIGMA_E              per level slot,
              sqrt(lam)*(loggf-NIST)/sigma_gf     per line slot,
              sqrt(ridge)*(theta-theta_HF)/(sigma_kind*|theta_HF|)  per param ]
        The last block is a per-parameter GAUSSIAN PRIOR pulling each radial
        parameter toward its ab-initio Hartree-Fock value with a per-kind width
        (PRIOR_SIGMA): a smooth, general form of Bob's FIXEDHF discipline (hard
        freeze = sigma->0; free fit = ridge=0). Levels are matched by IDENTITY
        each eval (not energy rank), so the vector has constant length and smooth
        entries -> 0.5*||vec||^2 is the regularized chi2."""
        self.run(x * scale)
        by_id = self._levels_by_identity()
        line_rows = self._current_line_rows()
        # global energy offset from the identity-matched levels
        eres = []
        for k, eobs in self.level_targets:
            L = by_id.get(k)
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
            theta = np.asarray(x) * scale
            denom = self.prior_sigma * np.maximum(np.abs(self.prior_centre), 1e-6)
            vec.extend(np.sqrt(ridge) * (theta - self.prior_centre) / denom)
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
            # same per-parameter HF prior as resid_vector (||.||^2 form)
            denom = fwd.prior_sigma * np.maximum(np.abs(fwd.prior_centre), 1e-6)
            chi2 += ridge * np.sum(((values - fwd.prior_centre) / denom) ** 2)
        return chi2

    return obj, scale


def fit_lm(fwd, lam, ridge, maxiter, presolve=True, progress=None,
           tol=1e-8, parallel_jac=None):
    """Levenberg-Marquardt fit on the stacked residual vector. Returns
    (best_values, result, scale, cov) where cov is the parameter covariance in
    SCALED units (x = value/scale), estimated from the Jacobian at the optimum
    (the basis for gf MC uncertainties). 0.5*||resid||^2 == the scalar chi2, so
    the optimum and lambda semantics match the Nelder-Mead path.

    LM is a LOCAL gradient method: from the bare RCE seed it under-converges
    (stops in a poorer basin) vs Nelder-Mead. So by default do a short NM
    PRE-SOLVE to land in the right basin, then LM to polish and -- critically --
    to produce the Jacobian/covariance. presolve=False uses LM alone (e.g. when
    the seed is already near the optimum).

    `tol` sets BOTH xtol and ftol. The default 1e-8 is publication-tight and
    makes the fit crawl for hundreds of extra evals once chi2 has plateaued (the
    last digits of a sub-cm^-1 move). For EXPLORATORY sessions pass a looser tol
    (e.g. 1e-5 / 1e-6): it stops the LM as soon as relative chi2/param change
    drops below tol, cutting wall-clock several-fold with no meaningful change to
    the answer. Tighten back to 1e-8 only for a final production fit."""
    scale = _scale_of(fwd)
    x0 = np.ones_like(fwd.seed)
    if presolve:
        obj, _ = make_objective(fwd, lam, ridge)
        nm = minimize(obj, x0, method="Nelder-Mead",
                      options={"maxiter": maxiter, "maxfev": maxiter,
                               "xatol": 1e-4, "fatol": 1e-3, "adaptive": True})
        x0 = nm.x
    # wrap the residual so we can emit progress (eval count + current chi2);
    # `progress(msg)` is any callable (e.g. append-to-file). Useful for the slow
    # full-basis fits where each eval is ~1s and a run is hundreds of evals.
    _last = {}
    def resid(x):
        v = fwd.resid_vector(x, scale, lam, ridge)
        if progress is not None:
            progress(f"eval {fwd.neval:5d}  chi2={float(v @ v):12.3f}")
        _last["x"] = np.array(x); _last["v"] = v
        return v
    # PARALLEL JACOBIAN: if a ParallelJac is supplied, hand scipy an analytic `jac`
    # that computes the n_param finite-difference columns CONCURRENTLY across worker
    # dirs (~Ncore x faster; bit-identical to scipy's serial 2-point differences).
    jac_arg = "2-point"
    if parallel_jac is not None:
        def jac(x, *a):
            # reuse the residual at x if scipy just evaluated it here
            f0 = _last["v"] if np.array_equal(_last.get("x"), x) else None
            return parallel_jac.jac(x, f0=f0)
        jac_arg = jac
    res = least_squares(
        resid, x0, method="lm", max_nfev=maxiter, xtol=tol, ftol=tol,
        jac=jac_arg, diff_step=1e-3)     # finite-diff step in scaled units
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
    ap.add_argument("--abinitio", default=None,
                    help="ab-initio ING11 for the ridge-prior centres (HF "
                         "values). Defaults to the seed if omitted.")
    ap.add_argument("--nist", required=True)
    ap.add_argument("--nist-lines", required=True)
    ap.add_argument("--free-kinds", default=None,
                    help="comma list restricting the free set, e.g. 'EAV' or "
                         "'EAV,P' (P=Slater/zeta, CI=interaction). Default: all.")
    ap.add_argument("--lambda", dest="lam", type=float, default=1.0)
    ap.add_argument("--ridge", type=float, default=0.0,
                    help="ridge-prior strength: pull each param toward its HF "
                         "value with per-kind width (PRIOR_SIGMA). 0 = free fit.")
    ap.add_argument("--maxiter", type=int, default=4000)
    ap.add_argument("--method", choices=["lm", "nelder-mead"], default="lm",
                    help="lm = Levenberg-Marquardt on the residual vector "
                         "(fast, gives the Jacobian/covariance); nelder-mead = "
                         "the old derivative-free scalar minimizer.")
    ap.add_argument("--out", default=None, help="write fitted ING11 here.")
    a = ap.parse_args()

    free_kinds = (set(k.strip() for k in a.free_kinds.split(","))
                  if a.free_kinds else None)
    fwd = Forward(a.run_dir, a.seed, a.nist, a.nist_lines,
                  free_kinds=free_kinds, abinitio_ing11=a.abinitio)
    print(f"{len(fwd.params)} adjustable params; seeded from {a.seed}"
          + (f"; ridge prior centred on {a.abinitio}" if a.abinitio else ""))
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
