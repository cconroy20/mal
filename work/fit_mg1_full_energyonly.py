#!/usr/bin/env python3
"""Full-basis (122-config) Mg I fit, BOB-STYLE PURE ENERGY (gf term OFF).

Bob fits OBSERVED ENERGY LEVELS only -- RCE has no transition data in its
objective. The gf-aware (lambda*gf) term is OUR Tier-2 addition, not his. This
driver turns it off (lambda=0) to (a) reproduce Bob's actual objective on the
full basis, and (b) test whether the gf term was helping or hurting, and whether
a pure-energy fit moves the 3p2 / singlet-1D CI that the gf term couldn't pull
(the energy residual on 3p2 is -1400 cm^-1 -- a strong, real signal).

Also applies --max-energy: high-Rydberg / near-IE levels are unreliable on this
basis and their rank-based slots collide -> bogus ~2000 cm^-1 residuals that
would DOMINATE a gf-off objective. Mg I IE = 61671.05 cm^-1.

Two variants via FREECI env (0 = Bob-pure, CI frozen; 1 = + free 12 CI):
  FREECI=0 python3 work/fit_mg1_full_energyonly.py   # Bob's exact recipe
  FREECI=1 python3 work/fit_mg1_full_energyonly.py   # pure-energy + freed CI
Watch: tail -f work/mg1_full/fit_<tag>_progress.log
"""
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import numpy as np
import ing11_params as IP
import make_report as R
from gf_fit import Forward, fit_lm
from gf_table import matched_table

# RCG does its I/O (ING11/OUTG11) inside RUN. To run two variants CONCURRENTLY
# they MUST use separate dirs (else their RCG calls clobber each other's OUTG11
# -> FileNotFoundError). RUNDIR env overrides; seeds are still read by abs path.
RUN = os.environ.get("RUNDIR", os.path.join(ROOT, "work", "mg1_full"))
CANON = os.path.join(ROOT, "work", "mg1_full")
NIST = os.path.join(ROOT, "data", "nist", "MgI_levels.tsv")
NL = os.path.join(ROOT, "data", "nist", "MgI_lines.tsv")

IE = 61671.05                      # Mg I ionization energy (cm^-1)
FREECI = os.environ.get("FREECI", "0") == "1"
TAG = "energyonly_freeci" if FREECI else "energyonly"
# STAGE=strong (phase 1 of a staged fit) writes to a distinct tag so it doesn't
# clobber the final phase-2 output; the orchestration seeds phase 2 from it.
if os.environ.get("STAGE") == "strong":
    TAG += "_phase1"
LOG = os.path.join(RUN, f"fit_{TAG}_progress.log")
nistl = R.load_nist_lines(NL)

# SEED env is the STARTING POINT of the LM. ABINIT env is the ridge-prior CENTRE
# (the physical HF background the ridge pulls toward); it defaults to SEED for the
# single-phase fit. STAGED FIT phase 2 MUST decouple them: seed from the phase-1
# result (good basin) while keeping the ridge centre at scaled-HF (the prior is a
# physical statement about HF, unchanged by phase 1). Set ABINIT=ING11.scaled there.
_seedname = os.environ.get("SEED", "ING11.scaled")
SEED = os.path.join(CANON, _seedname)
ABINIT = os.path.join(CANON, os.environ.get("ABINIT", _seedname))
OUTGINE = os.path.join(CANON, "OUTGINE.abinitio")   # for physical param names
# ZETA polish: loosen the ridge on spin-orbit so observed fine structure can pull
# zeta to truth (our HF zeta is below Bob's, a scale can't reach it; data can).
ZETA_SIGMA = float(os.environ["ZETA_SIGMA"]) if os.environ.get("ZETA_SIGMA") else None
# Release zeta ONLY for these configs (comma list of cfgkeys), the ones with
# observed resolved fine structure. Broad release ill-conditions the fit.
ZETA_CONFIGS = (os.environ["ZETA_CONFIGS"].split(",")
                if os.environ.get("ZETA_CONFIGS") else None)
# EAV polish: loosen the ridge on the CENTROID of named configs so the data can
# pull it far from the scaled-HF centre. Needed for the doubly-excited valence
# configs (3p3d, 3p4s, 3p4d) the parser fix exposed, whose HF centroid is ~7-20 kK
# off. Surgical (EAV_CONFIGS) -- a broad release de-anchors the well-placed bulk.
EAV_SIGMA = float(os.environ["EAV_SIGMA"]) if os.environ.get("EAV_SIGMA") else None
EAV_CONFIGS = (os.environ["EAV_CONFIGS"].split(",")
               if os.environ.get("EAV_CONFIGS") else None)

ND = ["3s3d", "3s4d", "3s5d", "3s6d"]
FREE_CI_PAIRS = ({
    frozenset(("3s2", "3p2")), frozenset(("3s2", "3d2")),
    frozenset(("3p2", "3d2")),
} | {frozenset(("3p2", nd)) for nd in ND}
  | {frozenset(("3d2", nd)) for nd in ND}) if FREECI else None

# FREE_CI env: free SPECIFIC named CI pairs, e.g. "3s2-3p2" or "3s2-3p2,3s2-3d2".
# Unlike the broad FREECI experiment, this targets one or a few. free_ci_sigma via
# FREE_CI_SIGMA (default 1.0). Used to let the data set the 3s2-3p2 CI (which sets
# the ground 3s2 energy) rather than transplanting Bob's value into a pinned zero.
if os.environ.get("FREE_CI"):
    FREE_CI_PAIRS = {frozenset(p.split("-"))
                     for p in os.environ["FREE_CI"].split(",")}
FREE_CI_SIGMA = float(os.environ.get("FREE_CI_SIGMA", "1.0"))


def log(msg):
    with open(LOG, "a") as f:
        f.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")


def gfrms(outg11):
    rows, _ = matched_table(outg11, nistl)
    d = [r["d"] for r in rows
         if r["nist"] >= -1 and r["acc"].strip()[:1] in ("A", "B")]
    return (np.sqrt(np.mean(np.square(d))), len(d)) if d else (float("nan"), 0)


def levelstats(fwd):
    res = np.array(fwd.energy_resid())
    return (float(np.sqrt(np.mean(res ** 2))), float(np.median(np.abs(res))),
            len(res))


def main():
    open(LOG, "w").close()
    t0 = time.time()
    cap = float(os.environ.get("MAXE", "60000"))   # basis-reliable cap (drops 3s.10d/11d)
    # STAGED FIT, phase 1 (STAGE=strong): free the EAV only of configs that have a
    # DIRECTLY-OBSERVED level at/below the cap (their centroid is pinned by data) --
    # plus 3p2, the one reliable above-cap perturber. This freezes the ~35 targetless
    # ("weak") EAVs (3p3d, 3p4s, the 3p.nd/3p.ns doubly-excited manifolds) at
    # scaled-HF for phase 1, so the well-determined structure converges in a
    # well-conditioned basin before the flat directions are released. Phase 2
    # (STAGE unset, seeded from ING11.phase1) frees the full set to polish + cov.
    _nist = R.load_nist(NIST)
    if os.environ.get("STAGE") == "strong":
        obs = {R._cfgkey(n["config"]) for n in _nist if n["E_obs"] <= cap}
        obs.add(R._cfgkey("3p2"))
    else:
        obs = {R._cfgkey(n["config"]) for n in _nist}
    # FIT_GF_SCALES=1: free EAV ONLY + 2 global G/F screening scales (Bob-style:
    # freeze structure, fit centroids + 2 screening factors). Far fewer params (75
    # vs 168), well-conditioned. Else the old free-all-EAV+P set.
    fit_scales = os.environ.get("FIT_GF_SCALES") == "1"
    free_kinds = {"EAV"} if fit_scales else {"EAV", "P"}
    # HYBRID (FREE_G env): individually free the G^k of named configs (e.g. the
    # 3s.np exchange series Bob frees), while the 2 scales handle the REST. Best of
    # both: per-config freedom where it matters + global screening elsewhere.
    free_p_keys = set()
    if os.environ.get("FREE_G"):
        import param_labels as PL
        want = {R._cfgkey(c) for c in os.environ["FREE_G"].split(",")}
        free_p_keys = {p["key"] for p in
                       PL.physical_params(SEED, OUTGINE, strict=False)
                       if p.get("kind") == "G" and p.get("cfg") in want}
    fwd = Forward(RUN, SEED, NIST, NL, free_kinds=free_kinds,
                  abinitio_ing11=ABINIT, obs_configs=obs,
                  free_ci_pairs=FREE_CI_PAIRS,
                  free_ci_sigma=(2.0 if FREECI else
                                 (FREE_CI_SIGMA if FREE_CI_PAIRS else None)),
                  max_energy=cap, outgine=OUTGINE, zeta_sigma=ZETA_SIGMA,
                  zeta_configs=ZETA_CONFIGS, fit_gf_scales=fit_scales,
                  free_p_keys=free_p_keys, eav_sigma=EAV_SIGMA,
                  eav_configs=EAV_CONFIGS)
    nci = sum(1 for p in fwd.params if p["kind"] == "CI")
    log(f"PURE-ENERGY (lambda=0)  max_energy={cap}  freeCI={FREECI} ({nci} CI)  "
        f"seed={_seedname}  zeta_sigma={ZETA_SIGMA}")
    log(f"  {len(fwd.params)} free params; {len(fwd.level_targets)} level targets "
        f"(<= IE)")
    e0, m0, n0 = levelstats(fwd)
    g0, ng0 = gfrms(fwd.outg11)
    log(f"SEED  levelRMS={e0:.1f} median|d|={m0:.1f} N={n0}  "
        f"gfRMS(A/B)={g0:.3f} (off-objective)")

    # lam=0 -> gf block zeroed; ridge keeps it well-posed
    # Loose convergence tol for EXPLORATORY runs (default 1e-6): stops the LM
    # once chi2 has plateaued instead of crawling the last sub-cm^-1 digits for
    # hundreds of extra evals. Override with TOL env (1e-8 for a production fit).
    tol = float(os.environ.get("TOL", "1e-6"))
    # PARALLEL JACOBIAN: NJAC workers compute the finite-difference columns
    # concurrently (~NJAC x faster, bit-identical). 0/unset = serial (old path).
    njac = int(os.environ.get("NJAC", "0"))
    pjac = None
    wdirs = []
    if njac > 0:
        from parallel_jac import build_worker_dirs, cleanup_worker_dirs, ParallelJac
        from gf_fit import _scale_of
        wdirs = build_worker_dirs(RUN, njac)
        pjac = ParallelJac(fwd, wdirs, fwd.init_kwargs, _scale_of(fwd),
                           lam=0.0, ridge=1.0, diff_step=1e-3, progress=log)
        log(f"(parallel Jacobian: {njac} PROCESS workers)")
    best, res, scale, cov = fit_lm(fwd, lam=0.0, ridge=1.0, maxiter=3000,
                                   presolve=False, progress=log, tol=tol,
                                   parallel_jac=pjac)
    if pjac is not None:
        pjac.close()
        from parallel_jac import cleanup_worker_dirs
        cleanup_worker_dirs(wdirs)
    log(f"(convergence tol={tol})")
    e1, m1, n1 = levelstats(fwd)
    g1, ng1 = gfrms(fwd.outg11)
    log(f"FIT   levelRMS={e1:.1f} median|d|={m1:.1f} N={n1}  "
        f"gfRMS(A/B)={g1:.3f}  nfev={fwd.neval}  "
        f"cov={'ok' if cov is not None else 'singular'}  ({time.time()-t0:.0f}s)")

    # Persist the fitted model. Use fwd.run(best) to (re)write the fitted ING11 AND
    # OUTG11 -- it is the ONLY correct writer when fit_gf_scales is on: the 2 SCALE
    # params are GLOBAL multipliers with no single ING11 lineno/col, so a direct
    # IP.write(fwd.params) raises KeyError('lineno'). run() strips the scales, applies
    # them to the frozen G/F slots, and writes the line-bound params. (A prior direct
    # IP.write here silently crashed every fit_gf_scales run AFTER the FIT log line,
    # leaving a stale on-disk model.) fwd.ing11 is the written ING11; fwd.outg11 the
    # RCG output.
    fwd.run(best)
    import shutil
    shutil.copyfile(fwd.ing11, os.path.join(RUN, f"ING11.{TAG}"))
    # COPY (not move) OUTG11 so fwd.outg11 stays in place for _block_levels() below.
    shutil.copyfile(fwd.outg11, os.path.join(RUN, f"OUTG11.{TAG}"))

    # worst levels, so we can see where a pure-energy fit still struggles. Uses the
    # IDENTITY matching (each target is a (cfgkey, termkey, Jstr); resolve it to the
    # current computed level with that identity) -- the old slot-based _level_at was
    # removed when the engine switched to identity matching.
    by_id = fwd._levels_by_identity()
    rows = []
    for k, eobs in fwd.level_targets:
        L = by_id.get(k)
        if L is not None:
            rows.append((L["E_calc"] - eobs, k[0], k[1], k[2]))
    if rows:
        offv = np.median([r[0] for r in rows])
        rows = [(r[0] - offv,) + r[1:] for r in rows]
        rows.sort(key=lambda r: -abs(r[0]))
        log("worst 10 levels (offset-removed):")
        for r in rows[:10]:
            log(f"   d={r[0]:+8.1f}  {r[1]:10} {r[2]:5} J={r[3]}")
    log(f"wrote ING11.{TAG} + OUTG11.{TAG}  -- DONE")


if __name__ == "__main__":
    main()
