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

## Full-basis (122-config) fit: completeness ALONE is not enough; CI matters

Grew Mg I to Bob's 122-config basis (tools/gen_in36.py, work/run_mg1_full.sh) and
fit Bob-style: free EAV+Slater of OBSERVED configs (170 params), freeze ALL CI +
unobserved-config params at HF, ridge prior toward HF (ridge=1, lam=3, LM).
Required first fixing parse_compositions for multi-chunk EIGENVECTORS blocks
(commit 5e567ff) -- it mislabelled levels on large bases (called the ground state
'3p4p'); now correct, gf matches to NIST 1->79.

RESULT (1217 evals, 62 min, well-posed cov):
  full-basis seed (ab initio): levelRMS 788   gfRMS(A/B) 0.178
  full-basis fit (EAV+P,ridge1): levelRMS 322  gfRMS 0.220   <- WORSE gf!
  our 9-config fit:            levelRMS 37    gfRMS 0.063
  Bob:                          levelRMS 12

Completeness alone made things WORSE, not better. Diagnosis (not a matching bug
-- matches are correct): most levels fit well (median |resid| ~20 cm^-1) but a
few have huge residuals, dominated by the GROUND STATE: full-basis ab-initio 3s2
1S sits ~1700 cm^-1 BELOW where the bulk of levels want the zero. The large basis
introduces strong CI (3s2-3p2-3d2 mixing) that DEPRESSES 3s2 -- and we FROZE all
CI at HF, so the fit cannot correct it. The worst gf residuals are the singlet
3s.nd 1D -> 3s3p 1P series (all too weak), driven by the same uncorrected
singlet-system CI.

KEY LESSON: freezing all CI works for Bob because his HF Slater/CI are SCALED
(his FIXEDHF carries scale factors ~0.6-0.85) so the frozen values already give
the right ground state; our raw HF CI does not. So either (a) free the key CI
integrals (3s2-3p2, 3s2-3d2, ... -- the ones our 9-config experiment found
HELP), or (b) adopt Bob's HF SCALE FACTORS for the frozen CI/Slater (extract from
his b*.log / hf*.dat). Completeness + HF-prior is necessary but needs the CI
treated, not frozen raw.

Also: the energy zero-point/offset (currently median residual) is fragile on the
large basis -- should anchor on the ground state, but that only matters once the
ground-state CI depression is fixed.

NEXT: extract Bob's HF scale factors, or free the dominant CI integrals on the
full basis and re-fit; compare gf to Bob.

## Tested: a single global scale factor (Bob's insight, simplified) -- doesn't work

Idea: replicate Bob's per-parameter FIXEDHF scale factors with ONE (or two)
global scale(s) on the frozen integrals, fit alongside EAVs. Tested cheaply by
scaling ab-initio integrals directly + re-running RCG (no optimizer needed):

  CI (R^k) global scale s:   3s2->3s3p 3P gap = 20376 (s=1) -> 19115 (s=0.80)
     [obs 21850] -- scaling CI DOWN moves the gap the WRONG way.
  Slater (F/G) global scale: 3P gap 20376->19613, 1P gap 35895->37742 as s:1->1.15
     [obs 21850 / 35051] -- 3P and 1P want OPPOSITE corrections; one scale can't.

At ab-initio the singlet-triplet SPLITTING (1P-3P, set by 3s3p G^1 exchange) is
already too large; a global scale moving F and G together can't shrink just the
exchange. Confirmed against our SUCCESSFUL 9-config fit: the params that moved
went in BOTH directions by very different amounts -- 3s5s-3p2 CI +274%, but
3s2-3p2 CI -59%; 3s3p G^1-type down 55%. A single (or two) global scale(s)
cannot reproduce that.

CONCLUSION: Bob's insight is genuinely PER-INTEGRAL, not a global screening
factor -- the HF errors point different ways for different integrals. The path
the evidence supports on the full basis = SELECTIVE freeing (Bob-style): EAVs +
a few low-config Slater (esp. 3s3p G^1) + the dominant CI (esp. 3s2-3p2, which
our 9-config fit cut 59%), ridge prior holding the rest at HF. NOT global scaling.

## What Bob actually freed (from c1200ez/oz.log) -- the general rules

Parsed his parameter status (col5 = uncertainty/step: nonzero => free; FIXEDHF =>
held at HF). Of ~200 params per parity he freed:

  EVEN: 57 free = 55 EAV + F2(22) + ZETA 2
  ODD:  87 free = 57 EAV + 30 non-EAV, the non-EAV being:
        G1(12)=3s-3p exchange, G1(13)=3s-np exchange for n=4..15 (one per member,
        values decay 23423->11.8 cm^-1, INDEPENDENT group codes), ZETA 3 series,
        plus G1(27),G1(27),F2(28),G1(28) for the 3p-nd / 3d configs.

GENERAL LESSONS (the craft, generalized):
1. Free EAVs almost universally -- one centroid per config that has an observed
   level. (~95% of his free params.) These are the best-determined; data fixes
   each config's absolute position.
2. Free EXCHANGE integrals G^k (esp. G^1) -- they set the singlet-triplet
   SPLITTING, which observed levels directly pin. He frees G^1 for the WHOLE
   Rydberg series (every observed 3s.np), not just low members; values decay
   smoothly with n. <-- this is exactly the integral our diagnosis found wrong
   (1P-3P splitting too large at HF).
3. Free SPIN-ORBIT ZETA -- sets fine structure (the J-splitting within a term),
   again directly observed. Freed per series member.
4. Free DIRECT F^k only for SAME-SHELL configs (F2(22)=3p^2, F2(28)=3d-nd type)
   where the direct integral is large and the term structure constrains it.
   NOT freed for 3s.nl configs (one electron in a closed-ish core -> F^k small).
5. Free ~NO CI (R^k) -- he holds interaction integrals at (scaled) HF. [Our raw
   HF CI doesn't work frozen; see prior section -- but the lesson stands that CI
   is the LAST thing to free, after exchange/spin-orbit/EAV.]
6. Each integral freed INDEPENDENTLY (distinct group codes), not tied across the
   series -- though the fitted values come out following a smooth Rydberg trend.

=> RECIPE for our full-basis fit: free EAV (observed configs) + G^1 exchange for
every observed 3s.np and 3s.nd + ZETA + F^2 for same-shell (3p^2,3d^2), ridge
prior on the rest. This is far more targeted than "free all EAV+P" (which freed
many irrelevant F^k and missed that exchange is the key lever), and it matches
both Bob and our successful 9-config fit. CI: revisit via scaled-HF, separately.
