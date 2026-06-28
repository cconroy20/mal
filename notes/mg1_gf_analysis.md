# Mg I gf: tightening the fit — analysis & findings

**Date:** 2026-06-28
**Goal:** improve the Mg I fitted log gf (vs NIST) from the warm-up RCE loop.

## UPDATE: combined energy+gf fit works (Tier-2 prototype)

Putting well-measured NIST gf into the objective alongside the energy levels
FIXES the degradation below. RCE itself can't do this (it fits energies only,
no transition data), so the fit is done by a Python optimizer wrapping RCG as a
black-box forward model (RCG is ~8 ms/call for Mg I -> cost is a non-issue):
trial radial params -> ING11 -> RCG -> parse levels+gf -> combined chi^2
(energy term, cm^-1, offset-removed; + lambda * gf term, NIST-accuracy weighted,
log gf >= -1 only). Tool: `tools/gf_fit.py` (+ `tools/ing11_params.py` for the
ING11 read/write, round-trip byte-exact).

Result (Mg I, seed = RCE-fitted ING11, lambda=1, Nelder-Mead), strong+A/B lines
(log gf >= -1), verified independently via gf_table on the final OUTG11:

    ab initio (RCG)        strong-gf RMS 0.114   level RMS  (unreferenced)
    RCE energy-only fit    strong-gf RMS 0.165   level RMS ~1038 cm^-1
    combined gf-fit        strong-gf RMS 0.071   level RMS    26.5 cm^-1

The combined fit improves BOTH at once: gf RMS 0.071 (below the ab-initio floor
of 0.114, and 57% better than energy-only RCE), AND levels to 26.5 cm^-1. The
weak lines deliberately excluded from the objective (e.g. 3s4s 3S -> 3s5p 3P,
NIST class C/D) are unchanged, as intended. This is the Tier-1.5/Tier-2 thesis
demonstrated end to end on one ion: a gf-aware fit recovers gf accuracy the
energy-only fit threw away.

### lambda sweep (tools/gf_fit_sweep.py, seed = RCE-fitted ING11)

    lambda  levelRMS(cm^-1)  strong-gf RMS  acc-wtd gf RMS
     0.1         9.0            0.160          0.168
     0.3        18.8            0.092          0.106
     1          26.5            0.081          0.086
     3          37.0            0.093          0.069   <- chosen
     10        167.8            0.089          0.071

Trade-off is clear: too small (lambda<=0.1) ~ energy-only (gf degrades); too
large (lambda>=10) blows up levels with no gf gain. lambda=1..3 is the sweet
spot. **Chose lambda=3**: minimizes the accuracy-weighted gf RMS (0.069, the
most defensible metric since it down-weights poorly-measured lines), levels
still excellent (37 cm^-1). Independently via gf_table on the final OUTG11:
strong+A/B gf RMS **0.063** (vs ab-initio 0.114 = 45% better), acc-wtd 0.098.
The report PDF gf page now shows this fit ("fitted (energy+gf)"), with the title
reporting strong-line RMS (log gf >= -1), not the weak-line-dominated all-line
RMS. Reproduce: `GFFIT=1 ./work/run_mg1.sh` (lambda via GFLAMBDA), then
`tools/report_mg1.sh`.

NEXT: add per-line gf MC uncertainties on top of this loop (the UQ deliverable);
try ridge-toward-ab-initio (--ridge) for ill-conditioned ions; carry to Fe II.

---

## TL;DR (original analysis that motivated the above)

1. The report's old gf line-matcher was **wrong** — it bucketed by (term-pair,
   J-pair) and took the nearest wavelength within 5%. That collides Rydberg
   series members of the same term (e.g. `3p^2 1S` vs `3s5s 1S`, `3p^2 1D` vs
   `3s3d 1D`) and manufactured large phantom residuals. Replaced with matching
   by **eigenvector-composition identity of both end levels** (config+term+J) —
   the same robust identity the RCE level fit is built on.
2. With the honest matcher, the surprising result: **the RCE fit DEGRADES gf
   relative to ab initio.** On strong, well-measured lines (log gf >= -1, NIST
   class A/B) the RMS goes from **0.114 (ab initio) -> 0.165 (fitted)**;
   accuracy-weighted RMS **0.122 -> 0.169**; plain RMS 0.19 -> 0.24.
   The energy fit buys better LEVELS at the cost of worse dipole matrix
   elements (gf is built from the eigenvectors, which the energy fit reshapes).
3. The large residuals are concentrated in the **singlet system and high
   Rydberg states**, which are also the worst-fit LEVELS (Δ_fit ~ 200–370 cm^-1
   for 3s4s 1S, 3s5s 1S, 3s4d 1D, 3s5p, 3p^2 1S) while the triplet ladder fits
   to ~0. Same root cause for both: the singlet/Rydberg block is under-fit, with
   `3p^2 1S` (Δ_fit = -333) the smoking gun (its CI drags the whole singlet
   system via 3s^2 and 3s_nd 1D).
4. Freezing the under-constrained Rydberg exchange integrals (3s4d G2, 3s5p G1)
   at their ab-initio values — they otherwise run to unphysical values, e.g.
   3s4d G2 -> ~0, 3s5p G1 -> 23% of ab initio — does **not** help gf (strong-line
   RMS 0.153–0.165; overall worse). So the gf degradation is broader than two
   pathological parameters; it is the eigenvector reshaping itself.

## What this means for the project

- A naive "fit the levels, take the gf" pipeline can make gf WORSE than the
  ab-initio RCG values for some lines. This is exactly the motivation for the
  Tier-2 modernization: the fit needs to be **gf-aware / regularized**, not a
  pure energy least-squares. Options to explore:
  - regularize radial params toward ab initio (ridge / L2), so the fit can't
    reshape eigenvectors arbitrarily to chase a few energies;
  - include line data in the objective (Tier-2 Bayesian: spectrum in the
    likelihood), coupling gf to the fit directly;
  - per-line model: report ab-initio gf where the fit demonstrably degrades it.
- The honest evaluation metric matters: weak (log gf < -1) and NIST class C/D
  lines carry most of the apparent RMS but little spectroscopic weight. Use the
  accuracy-weighted / strong-line RMS as the headline number.

## Tooling added

- `tools/parse_cowan.py::identify_lines()` — attach robust upper/lower level
  identities (config, term) to each E1 line by matching its E_low/E_up to the
  eigenvector-composition level table. The OUTG11 transition-block labels are
  unreliable (block-header config); this uses the dominant eigenvector instead.
- `tools/make_report.py::match_gf_by_identity()` — gf↔NIST matching by both end
  levels' identity (replaces the old `_match_gf` term-pair+wavelength matcher).
  The report's gf page now uses it.
- `tools/gf_table.py` — per-line gf comparison table (fitted & ab initio vs
  NIST) with plain, accuracy-weighted, and strong+A/B RMS. The diagnostic used
  to locate which lines drive the RMS.
- `tools/build_ine20.py --freeze-params 'cfg:PARAM,...'` — hold specified
  single-config physical params (e.g. `3s4d:G2,3s5p:G1`) at ab-initio values
  instead of freeing them. The surgical inverse of `--free-ci-pairs`. (Tried for
  Mg I gf; did not help — kept as a tool for under-constrained ions.)
- `work/run_mg1.sh` — now snapshots a pristine `OUTGINE.abinitio` (build reads
  this, not the clobbered OUTGINE) and honors `FREEZE='cfg:PARAM,...'`.

## Learning from Bob's Mg I deck (b1200e/o.com, c1200ez/oz.log)

Compared his actual production decks to ours to understand his approach.

**Scale.** Bob: 122 configs (61 even + 61 odd), ~185 observed levels, ~1340
radial params in the even-block Hamiltonian. Us: 9 configs, ~28 levels, ~35
params. So "Bob has more free parameters" is misleading -- see below.

**He fits a huge config space but FREEZES almost all of it at Hartree-Fock.**
Of ~1340 even params he frees only ~57; the rest are FIXEDHF (held at ab-initio
HF, often with a fixed scale ~0.6-0.85). His free set is overwhelmingly EAVs:
EVEN free = 55 EAV + 1 F + 1 zeta; ODD free = 57 EAV + 16 G + 13 zeta + 1 F.
He frees a Slater G/zeta only where a multiplet splitting is resolved, and frees
~ZERO CI integrals. EAVs are the best-determined params (our Jacobian SVD agrees).

**We do the opposite:** free 12 CI + 15 Slater on 28 levels -- exactly the
ill-determined directions Bob pins. That drove our rank-deficiency + gf degrade.

**BUT, tested on our small basis (LM, seed=RCE fit, lambda=3):**
  free ALL (35)        rank 8  levelRMS 38.5  gfRMS(A/B) 0.063
  EAV only (8)         rank 6  levelRMS 642   gfRMS 0.143
  EAV+Slater/zeta (23) rank 8  levelRMS 122   gfRMS 0.141
  EAV+CI (20)          rank 8  levelRMS 78.5  gfRMS 0.084
In OUR 9-config model freeing CI HELPS (free-ALL/EAV+CI best) -- opposite of Bob.
Reason: with only 9 configs the CI integrals stand in for missing interacting
configs. Bob doesn't free CI because he INCLUDES the perturbers explicitly
(3d4s, 3d2, 3p3d, 3d4p, ... -- exactly the doubly-excited 3d configs we lack and
that perturb our worst-fit terms 3s_nd 1D / 3p2). All our recipes stay rank-
deficient (global offset + weakly-observed EAVs).

**Conclusion / strategy.** Bob's quality comes from a LARGE, COMPLETE config
basis, not from clever freeing. His freeze-at-scaled-HF discipline is what makes
a large basis tractable. The lever for us is COMPLETENESS: enlarge the basis
(add 3d4s, 3d2, 3p3d, 3d4p, higher Rydberg), then -- to keep it well-posed --
freeze radial structure at HF and free mainly EAVs (+ selective zeta/G), i.e.
ridge-toward-ab-initio (our --ridge), which Bob's hand-selection approximates.
Freeing CI is a small-basis crutch to retire once the basis is complete.
NEXT: grow the Mg I config list toward Bob's, re-fit Bob-style, watch gf.

## Ridge priors toward ab-initio HF (generalizing Bob's FIXEDHF)

Rather than hard-freeze parameters (Bob's FIXEDHF binary), use a per-parameter
GAUSSIAN PRIOR toward the ab-initio HF value -- strictly more general (hard
freeze = sigma->0; free fit = ridge 0) and it makes the rank-deficient fit
well-posed. gf_fit now: prior CENTRES = ab-initio ING11 values (--abinitio),
per-kind widths PRIOR_SIGMA = {EAV 1.0 (loose), P/Slater 0.30, CI 0.15
(tightest)}, scaled by --ridge. resid_vector appends sqrt(ridge)*(theta-theta_HF)
/(sigma_kind*|theta_HF|) per param.

Mg I (LM, seed=RCE fit, lambda=3), free-all, HF-centred prior:
  ridge  rank  levelRMS  gfRMS(A/B)  cov
  0.0     8/35   38.5     0.063      singular
  0.1    35/35  600       0.171      ok   (too weak: levels drift)
  1.0    35/35   46.1     0.077      ok   <- operating point
  10.    35/35   64.7     0.093      ok
ridge~1 keeps ~free-fit quality (gf 0.077 vs 0.063, level 46 vs 38) but is now
FULL-RANK with a finite covariance -- the prerequisite for MC gf uncertainties,
and a principled stand-in for hand-selecting free vs fixed. Knobs: --ridge,
--abinitio (prior centre), --free-kinds (optional hard restriction), PRIOR_SIGMA.
