# Mg I gf: tightening the fit — analysis & findings

**Date:** 2026-06-28
**Goal:** improve the Mg I fitted log gf (vs NIST) from the warm-up RCE loop.

## ===== 2026-06-30 (v3b): THE 3s.nd 1D GAP — PROBED, MECHANISM FOUND, NOT YET CURED =====
Goal for this pass: "do at least as well as Bob." Localized WHERE the gap is and
tested (forward-model probes, isolated RCG run dir) every easy explanation. Many of
the session's own earlier assertions were FALSIFIED here — trust the tested results
below, not the narrative that preceded them.

WHERE THE GAP IS (RMS decomposition of our staged fit, levels <=60000):
  6 SINGLET 3s.nd 1D levels carry RMS 305; the OTHER 96 levels are RMS 43 (median
  2.8) — already ~Bob. So "beat Bob" == "fix the 3s.nd 1D series", nothing else.
  On the SAME 6 levels: OURS RMS 305 (3s3d 1D +520, 3s4d +429, 3s5d +215, ...),
  BOB RMS 59 (+77, -13, +52, ...). NB Bob does NOT nail them either — it's his
  worst series too; we're just 5x worse.

MECHANISM — CONFIRMED by probes:
  - Scaling CI(3s.nd-3p2/3d2) moves the 1D levels but NEVER the 3D (3D delta=0 in
    every probe). Clean proof: it's 1D-perturber (3p^2/3d^2 1D) MIXING, exactly the
    adversarial-analysis hypothesis, now tested not asserted.
  - The per-member CI response is DIAGONALLY DOMINANT (offdiag/diag ~0.25 at +10%):
    each 3s.nd 1D is moved mainly by its OWN CI -> a joint CI fit WOULD be well-posed.

FALSIFIED (do not re-try these — tested and wrong):
  - "Bob curates the 3s.nd-3p2 CIs" — NO. His are FIXEDHF 0.800, ~= ours. 3s3d-3p2:
    Bob 10513.8 vs our 10656.7 (1.4% apart). Same CI, 5x different residual.
  - "Match Bob's 3p2 EAV (59468)" — makes the 1D series WORSE (3s3d 1D +520 -> +740).
  - "Free F2(3p2) to place the 1D perturber independently" — NO. Scaling F2(3p2)
    moves the OBSERVED 3p2 3P/1S strongly but the 3p2 1D perturber only ~700 cm^-1
    over a huge F2 swing, and the 3s.nd 1D barely moves. No F2 escape hatch.

THE ACTUAL TENSION — tested:
  Our fit RAISES 3p2 EAV by +1930 cm^-1 (55.57 -> 57.50 kK) to match the OBSERVED
  3p2 3P@57813 & 1S@68275 (fitted residuals -34/-29, good). That raise pushes the
  UNOBSERVED 3p2 1D perturber up (our fit puts it ~65700), dragging the 3s.nd 1D
  series up with it. Forward-model proof of the tension: forcing 3p2 EAV back to
  ab-initio 55.57 improves the 1D series (3s3d 1D 520->253) but blows up 3p2 3P/1S
  to -1945/-1753. ONE 3p2 EAV parameter cannot fit both the observed 3p2 levels AND
  keep the 1D perturber low, and F2 gives no second handle. (Honest nuance: vs
  AB-INITIO the fit actually IMPROVES the 1D-series RMS 500->305 by fixing the one
  bad member 3s3d [-1202->+520]; it DEGRADES the other 5, e.g. 3s4d 1D -6 ->+429.)

THE OPEN CRUX (NOT resolved): Bob has the SAME basis, SAME CIs, and ALSO raises his
3p2 EAV to 59468 — yet his 1D series is 5x better. Something in his fuller deck
relaxes the tension; we don't know what. THE decisive test is a full RCE
REPRODUCTION: transcribe Bob's COMPLETE fitted deck (all ~100 integrals from
c1200{e,o}z.log) into our ING11, run OUR RCG, and check whether we reproduce his
1D residuals (isolates his-parameters vs our-forward-model). Deferred — it's a
large, error-prone transcription; do it deliberately in a fresh pass.
PRAGMATIC ALTERNATIVE (untested): re-weight the objective to stop the LM sacrificing
the 1D series (up-weight 3s.nd 1D or down-weight 3p2 3P/1S) — may recover much of
the gap without explaining Bob.
PROBE HARNESS: the scratchpad probe.py pattern (copy full RCG deck to an isolated
dir, perturb ONE ING11 param by IP.parse key, re-run rcg, read 1D/3D residuals) is
the right tool — cheap (~1 RCG call/probe), and it caught two false-negatives from
wrong param-key formats (param_labels keys != IP.parse keys: IP.parse CI keys have
NO 'Mg I' prefix, e.g. '3s3d    - 3p2|CI0'). Recreate under work/mg1_probe.

## ===== 2026-06-30 (v3): THE PARSER WAS DROPPING 83% OF THE LEVELS =====
**Headline: the "data starvation → enlarge the basis" diagnosis was wrong. The
basis was Bob's full 122-config deck all along; `parse_compositions` was silently
parsing only 192 of 1105 computed levels (64 of 122 configs).** Fixing one line
recovered the rest. This supersedes the "ENLARGE THE BASIS is the master lever"
conclusion in earlier sessions — that was chasing a shadow cast by a parser bug.

**The bug** (`tools/parse_cowan.py`, `parse_compositions`): RCG prints each
J-block's eigenvalues 11-per-line **separated by blank lines**, ended by
`CONFIG. NO.`. The eigenvalue-read loop had `elif s.strip()=="" and evs: break`
— it stopped on the FIRST blank line, keeping only the first ~11 eigenvalues of
every J-block and discarding the rest. Small bases (9-config: one row per block,
no wrap) never triggered it, so it hid through all the bulletproofing. On the
full basis it dropped every high-lying doubly-excited level — 3p3d, 3p4s, 3d.nl —
i.e. exactly the above-IE valence configs Bob fits and the perturbers the
singlet-1D mixing needs. **Fix:** terminate the eigenvalue read only on the real
section headers (`CONFIG. NO.` / `EIGENVECTORS` / `G-VALUES`) or the next
`EIGENVALUES(J=`), never on a blank line.

**Effect of the fix (OUTG11.abinitio, full 122-config basis):**
  - levels parsed:   192  → 1105   (per-J counts now equal the RCG matrix dims)
  - distinct configs: 64  →  122   (the whole deck, incl. 3p3d @ 75–86 kK)
  - NIST matches < IE: 98 →  131
  - above-IE valence now matchable: +27 (3p2 1S@68275, 3p4s, the four 3p3d
    @80693–85925 Bob fits). These are cap-exempt (n≤4) so the fit picks them up
    automatically — level targets 87 → 111 at cap 60000, no driver change.

**Regression guard:** `tools/test_parse_cowan.py` — per-J-block invariant
`n_eigenvalues == n_CONFIG.NO` (both = the matrix dimension; diverge iff the read
truncates), plus full-basis (≥122 configs / ≥1000 levels, 3p3d present) and
9-config-unchanged (exactly 30 levels) checks. Verified it FAILS under the old
parser (first block 11 vs 47). Runs standalone: `python3 tools/test_parse_cowan.py`.

**On the earlier "-14000 cm⁻¹ basis-limit artifact" and "Bob uses quantum-defect
extrapolation to n=59":** both were artifacts of the truncated view / secondhand
notes. Checked against Bob's own RCE logs (`c1200{e,o}z.log`): Bob fits **257
levels, max E_obs 83537**; ABOVE the IE he fits exactly **4 levels, all 3p3d**
(valence, in the shared deck), NOT high Rydberg. His highest Rydberg member is
**3s.15p** (deck stops at 3s.12p; 13–15p leak from the truncated orbital set, same
as ours). The "n=59 / 98899" figure was the NIST catalog range, not Bob's fit. So
Bob's page-7 reach comes from (a) not capping the naturally-produced 3s.13–15p and
(b) fitting the handful of above-IE 3p3d — both now available to us post-fix.

**The energy cap, corrected:** raising the fit cap from ~58000 to 60000 admits the
well-placed n≤8 high-Rydberg levels (they sit ≤190 cm⁻¹ at ab-initio; the cap had
been needlessly discarding ~38 good levels). Stop at 60000, NOT the IE: the topmost
3s.9d/10d Rydberg members are the last members of the truncated series and come out
2000–12000 cm⁻¹ low (basis-edge, same wall Bob has) — including them explodes RMS.

Head-to-head vs Bob, E_obs ≤ 60000, offset-removed (pre-parser-fit numbers; re-fit
with the +24 exposed targets pending): median |E-O| ours 2.4 vs Bob 0.0; RMS 53 vs
18. Bulk is near-Bob; the ~3× RMS gap is the singlet 3s.nd 1D / near-IE tail
(unobserved-perturber wall — curated CI, not more data). Now that 3p3d/3d² are
VISIBLE, the perturber levers the adversarial analysis wanted are finally testable.

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

## Cross-species check: Fe II (open 3d^7) -- the rules generalize but EMPHASIS shifts

Parsed Bob's Fe II logs (c2601e/o.log; NPAR 2645 even / 2393 odd, ground = 3d^7).
Free counts (corrected for per-config repeats):

  Fe II EVEN: ~462 fixed, free = ALPHA x47, BETA x46, ZETA x9, F2/F4(11) [3d-3d]
              several, plus scattered G^k/F^k for 3d-nl; relatively FEW EAVs (~13).
  Fe II ODD:  similar (ALPHA/BETA dominant, F/G of the d-shell, few EAV).

CONTRAST WITH Mg I:
  Mg I (closed-shell + 1-2 valence e-): free set ~95% EAV; key non-EAV = G^1
     EXCHANGE (singlet-triplet splitting) + ZETA. No ALPHA/BETA. No F^k(dd).
  Fe II (open 3d^7): free set dominated by ALPHA/BETA (Trees effective-operator
     CI corrections for the d^n shell) + F^2/F^4(3d,3d) (which BUILD the d-shell
     term structure) + G^k + ZETA. EAVs a minority.

GENERAL PRINCIPLE (holds for both): free the radial integrals that the OBSERVED
LEVEL STRUCTURE most directly constrains, freeze everything else at (scaled) HF.
WHICH integrals those are is dictated by the open shell:
  - closed/few-electron  -> EAV centroids + exchange G^k + spin-orbit ZETA.
  - open d^n (iron group) -> intra-shell F^2/F^4(dd) (term structure) + ALPHA/BETA
    (far-config CI as effective operators) + G^k + ZETA; EAVs secondary.
The 'free exchange/spin-orbit, free same-shell direct F, freeze CI' rules from Mg
I are a SPECIAL CASE; the iron-group case adds the dd-direct integrals and the
Trees alpha/beta operators because the open shell makes them the dominant,
observable degrees of freedom.

IMPLICATION for our engine: the free-set SELECTOR must be config/shell-aware, not
a fixed kind list. A good general heuristic: free EAV(observed) + ZETA + G^k +
F^k of the OPEN/same shell + (for d^n/f^n) ALPHA/BETA; ridge-prior the rest. This
covers Mg-I-like and Fe-II-like ions with one rule. (Cowan RCE supports
ALPHA/BETA; our gf_fit/ing11_params must learn to read/free them.)

## Ruleset selector implemented; full-basis OUTGINE group-code walk still blocking

Implemented the physics ruleset (notes/bob_fit_ruleset.md) as build_ine20
--ruleset: _ruleset_free(name,cfg,observed,open_ls,is_ion,orbitals) decides
free/freeze by family + shell. VERIFIED CORRECT in isolation: on full-basis Mg I
block 0 it would free 31 EAV + 49 G^k + 31 ZETA + 3 F^k (and never CI/ALPHA/BETA
since Mg I is closed-shell) -- exactly the ruleset intent. CLI infers ground
config/open-shells/observed-configs from NIST automatically (--ion flag for BETA).

Fixed two real bugs found along the way:
- _cfgkey multi-digit n: '3d10d'->'3d.10d' (was '3d'), '3s11d'->'3s.11d'. Same
  bug in BOTH build_ine20._ORB and make_report._ORB; fixed both to agree (they
  must match for config-key matching). 9-config report unchanged.
- _is_groupcode_line: was matching FLOAT eigenvalue lines ('129.5661...') because
  they contain integers >=100; now excludes any line with a '.'.

REMAINING BLOCKER: on the 122-config OUTGINE the group-code SECTION walk doesn't
align name<->code across the (apparently multi-section / interleaved) full-basis
layout -- the ruleset freeings compute correctly but don't all land in the output
INE20 (audit sees only a tiny slice freed). The 9-config OUTGINE is a single
clean block and works; the full-basis layout differs and needs the block/group-
code walker generalized (likely multiple group-code runs per parity block, or
the param-name header vs group-code ordering diverges at scale). 9-config
pipeline regression-tested OK (AVDEV->0).
NEXT: generalize _build_focused's group-code section handling for large OUTGINE;
then re-run full-basis Mg I with --ruleset and the HF ridge prior.

## Full-basis RULESET fit (first run) -- confirms Rule 2 (scaled HF) is load-bearing

Ran the bulletproofed build_ine20 --ruleset on full Mg I: freed 243 = 66 EAV +
100 G^k + 74 ZETA + 3 F^k (no CI/ALPHA-BETA; correct for closed shell), then RCE.
RESULT: RCE did NOT converge cleanly (AVDEV oscillated ~4.5-56 kK), and gf got
WORSE: seed strong+A/B 0.178 -> ruleset-fit 0.439. Ground 3s2 1S = -8714 cm^-1,
3s2->3s3p 3P gap = 24217 (obs 21850, now 2367 TOO LARGE).

Diagnosis: SAME root cause as the earlier "free all EAV+P" full-basis failure --
the frozen CI is at RAW HF, which structurally distorts the gap; RCE freed the
right things (EAV/G/ZETA per ruleset) but can't repair a gap error that lives in
the frozen CI. The fit thrashes trying to compensate through the wrong knobs.
(Also: subst_fitted_params only handled 61 per-config lines, not the freed
within-config G/ZETA/F -- needs generalizing -- but the RCE instability is the
real blocker, not the substitution.)

CONCLUSION: this is exactly Rule 2 of our own ruleset -- freeze at SCALED HF
(~0.85), not raw HF. The ruleset's FREE set is validated (correct families), but
the FROZEN background must be scaled first. The two are inseparable: completeness
+ correct free set + scaled-HF frozen background. NEXT: apply a per-integral (or
uniform ~0.85 CI/Slater) HF scale to the OUTGINE before --ruleset + RCE; extract
Bob's actual scales from b*.log/hf*.dat for the per-integral version.

## Full-basis RCE instability ROOT CAUSE: build_ine20 multi-line J-block bug

The full-basis ruleset RCE oscillated (AVDEV ~5 kK) NOT because of the HF scale --
it's a build_ine20 level-SUBSTITUTION bug on large J-blocks. A single (parity,J)
block in the 122-config OUTGINE has MANY value lines (e.g. J=0 even: ~90 computed
T-values across lines 479-491) followed by MULTIPLE flag lines (492-494...). The
value/flag substitution loop assumes ONE value line + ONE flag line per J, does
j+=2, Jval+=1 -- so after the first line it treats each continuation value line
as a new J and mis-substitutes. Result: unmatched high-energy computed values
(~130 kK) get POSITIVE include-flags -> RCE told to fit levels to 130000 cm^-1 ->
AVDEV ~5 kK, oscillating. (Confirmed: INE20 value/flag pairs past the first show
129-131 kK values all flagged positive.)

This is the SAME class as the eigenvector-chunk and group-code bugs: the 9-config
layout is one-line-per-J, the 122-config layout wraps, and the parser assumed the
small layout. FIX NEEDED: rewrite the value/flag substitution to gather ALL value
lines for a J then ALL flag lines, substitute across the full set, before moving
to the next J. (HF-scale question is downstream of this -- can't judge the fit
until the targets are correct.)

Good piece kept: build_ine20 --max-energy <cm^-1> caps the fit to the bound
spectrum (our NIST cache carries 90+ autoionizing/high-Rydberg levels above the
IE that shouldn't be fit over an incomplete basis). Verified: 284 -> 204 levels
below the Mg I IE. Necessary but not sufficient -- the J-block bug dominates.

## DECISIVE: full-basis fit, correct deck + ridge, STILL loses to 9-config -> scaled HF is THE remaining lever

Ran the full-basis fit via gf_fit (RCG forward model + HF ridge prior, NOT RCE),
on the now-correct deck, free EAV+Slater of observed configs (167 params),
CI frozen at raw HF, ridge=1, lam=3.

  seed (ab initio)   gfRMS(A/B) 0.178
  full-basis fit     gfRMS(A/B) 0.220   levelRMS 322   <- still worse
  our 9-config fit   gfRMS(A/B) 0.063
  Bob                levelRMS 12

RIDGE WORKED: chi2 26183->9066, cov=ok, NO divergence (RCE diverged to 229 kK on
the same deck -- regularization is exactly why gf_fit succeeds where RCE fails).
DECK IS CORRECT (clean optimizer gives the same 0.220 as before the parser fixes,
confirming those fixes were faithful). Free-set is the ruleset's.

So with completeness + correct deck + correct free-set + working regularization
all in place, the full basis STILL loses to 9-config. By ELIMINATION the remaining
cause is the one thing not yet done: CI frozen at RAW HF (ruleset Rule 2 = SCALED
HF). Signature confirmed: full-basis-fit 3s2->3s3p gap = 23435 (obs 21850, +1585),
ground 3s2 = -4884 -- the raw-HF-CI distortion the frozen CI can't let the fit
correct. This is now ISOLATED beyond doubt.

NEXT (unambiguous): apply the HF scale to frozen CI. Global scale already shown to
fail (errors are per-integral), so extract Bob's PER-INTEGRAL HF scale factors
from b1200*.log / hf1200z.* and apply them, then re-fit. That is the last missing
ingredient of the ruleset (Rule 2).

## Scaled-HF fit: fixes ENERGIES, not gf -> the limiter is per-integral scale granularity

Applied Bob's HF screening (CI x0.75, Slater x0.8; tools/scale_hf.py) to the
frozen background, re-fit (ridge, ruleset free-set). Results:
  scaled-HF fit: levelRMS 324, gfRMS(A/B) 0.212  (raw-HF was 0.220)
  Rule 2 WORKED for energies: 3s2->3s3p gap 23435 -> 21476 (obs 21850, now great);
    ground depression gone (-2926). chi2 13951->5843 (vs raw 26183->9066).
  But gf barely moved (0.225->0.212), still >> 9-config 0.063.

Apples-to-apples: on the 7 lines BOTH fits cover, 9-config RMS 0.063 vs full-basis
0.225. So full-basis is genuinely worse ON THE SAME LINES, not just covering more.
Worst residuals = singlet 3s.nd 1D -> 3s3p 1P and 3s2 -> 3s4p 1P (all ~0.5 dex
too weak).

KEY COMPARISON vs Bob's OWN gf on those exact lines:
  3s2->3s4p 1P  (nist -0.95):  Bob d=-0.13,  us d=-0.57
  3s4d 1D->3s3p 1P (nist -0.50): Bob d=-0.14, us d=-0.51
Bob gets them RIGHT (~0.13) with the SAME complete basis and SAME freeze-CI
philosophy. So the gap is NOT the basis, NOT the CI-freezing rule, NOT ridge --
it is that Bob's frozen parameter VALUES are better than our uniform-per-family
scale. Bob uses PER-INTEGRAL scales (G^k at 0.6 vs 0.8, CI at 0.7 vs 0.8); our
flat CI 0.75 / Slater 0.8 is too coarse for these sensitive singlet 1D/1P/1S
interactions. (Also: our 9-config win was partly because it FREED CI to absorb
missing-config error -- a small-basis crutch that happened to nail 7 lines.)

CONCLUSION: the engine is complete & correct; the remaining gap to Bob is
PER-INTEGRAL HF SCALE FIDELITY. NEXT: read Bob's exact per-integral scale (and
his fitted G^k) from c1200*.log and transcribe them onto our integrals (he gives
fitted_value, HF_value, and scale per line), rather than a flat per-family scale.
That is the difference between "Bob-style" and "Bob's actual numbers".

## Transcription test: it's the frozen CI STRUCTURE, not single scales or G^1

Compared our HF integrals to Bob's directly (from c1200*.log). 3s3p G^1(12):
our HF 27885 cm^-1 vs Bob HF 29818 (6% -- different HF codes, minor); Bob FREED
and fit it to 23423 (0.785xHF). Tested setting our 3s3p G^1 to Bob's exact 23423:
gf moved only 0.212->0.205. The worst lines (3s2->3s4p 1P d=-0.53, 3s4d 1D->3s3p
1P d=-0.50) DID NOT improve. So it is NOT the 3s3p exchange value.

Also found: our scaled-fit barely moves G^1 from the scaled-HF prior center
(22308->22285), even though it's in the free set with sigma=0.30 that would
ALLOW the move -- because the gf targets that constrain those singlet lines
aren't pulling it (the data signal is weak/absent in our objective for the 1D/1P
Rydberg mixing).

FINAL SCORECARD (strong+A/B gf RMS):
  ab-initio raw HF        0.178
  full-basis fit raw      0.220
  full-basis fit scaled   0.212
  our 9-config fit        0.063   <- still best, on its 7-line subset

DIAGNOSIS (settled): the persistent ~0.5 dex on singlet 3s.nd 1D->3s3p 1P and
3s2->3s4p 1P comes from the 1D/1P EIGENVECTOR MIXING (3s.nd / 3s.np Rydberg with
3p^2/3d^2), which is set by the FROZEN CI (R^k) STRUCTURE. Bob gets these right
(~0.13 dex) with his per-integral frozen CI; our uniform-scaled CI doesn't
reproduce that mixing, and freezing CI means the fit can't fix it. Our 9-config
win was precisely because it FREED the 3s2-3p2 / 3s.ns-3p2 CI (small-basis crutch)
-- which reshaped exactly this mixing.

STRATEGIC CONCLUSION: to match Bob we need his EXACT per-integral frozen CI
values (5567 of them), transcribed from c1200*.log -- a large, fragile parse
(config inheritance, no names on CI lines) with uncertain payoff -- OR we accept
that for OUR engine, selectively FREEING the key low CI integrals (3s2-3p2,
3s2-3s.nd, 3p2-3s.nd) is the pragmatic lever (it worked at 9 configs and directly
targets the 1D/1P mixing), even though Bob freezes them. The ruleset's
"never free CI" is right for BOB (his CI values are curated); for us, freeing the
FEW dominant low-l CI integrals is the equivalent of his curation. This is the
real, generalizable lesson: free what you can't otherwise get right.

## Path 2 TESTED: selectively freeing the dominant low-l CI -- the prior isn't the limiter, the DATA SIGNAL is

Implemented free_ci_pairs in gf_fit.Forward (config-NAME matching: CI keys carry
the readable pair, e.g. "3s2 - 3p2|CI0", cleaner than build_ine20's numeric index
which only handled 1-digit configs). Driver work/fit_mg1_full_freeci.py: scaled-HF
background + EAV+P ruleset free set + 12 freed CI = the dominant singlet-system
integrals {3s2-3p2, 3s2-3d2, 3p2-3d2, 3p2/3d2 - 3s.nd(n=3..6)}, the exact
integrals that set the 1D/1P/1S mixing. Two runs, both ridge=1 lam=3 LM, well-posed cov:

  scaled-HF, CI frozen           gfRMS(A/B) 0.212  levelRMS 324
  + free 12 CI, ridge sigma=0.15 gfRMS(A/B) 0.208  levelRMS 313   (chi2 13951->5427)
  + free 12 CI, ridge sigma=2.0  gfRMS(A/B) 0.215  levelRMS 317   (chi2 13951->5627)

RESULT (tight ridge 0.15): freeing the CI moved the WORST lines in the RIGHT
direction but only slightly -- 3s2->3s.4p 1P d=-0.70->-0.58, 3s.4d 1D->3s.3p 1P
-0.52->-0.46 -- while 3s.5d 1D->3s.3p 1P got slightly WORSE (-0.35->-0.41); net
RMS 0.212->0.208. The freed CI barely moved: 3s2-3p2 +4%, 3s3d-3p2 -4%,
3d2-3p2 -13%, rest <2%. CONTRAST: our 9-config fit cut 3s2-3p2 by 59%.

WHY so little movement? Two candidates: (a) the CI ridge prior (sigma=0.15,
tightest) is PINNING them near scaled-HF; (b) the gf DATA SIGNAL is too weak to
pull them (only ~20 A/B lines; the level-energy chi2 term, which dominates,
doesn't care about the 1D/1P mixing -- same symptom flagged for G^1 above: "the
data signal is weak/absent in our objective for the 1D/1P Rydberg mixing").
DISCRIMINATING TEST (DONE): re-fit with the freed-CI ridge LOOSENED to sigma=2.0
(effectively free, via new free_ci_sigma override that relaxes ONLY the
hand-picked CI, leaving the other ~5500 pinned).

VERDICT -- IT IS THE DATA SIGNAL, NOT THE PRIOR. With the prior essentially
removed the CI moved EVEN LESS (3s2-3p2 +2% vs +4%; 3d2-3p2 -7% vs -13%) and gf
got slightly WORSE (0.208->0.215): the loose ridge just let the fit wander a bit
in unhelpful directions. The worst singlet lines are unchanged (3s2->3s.4p 1P
d=-0.59, 3s.4d 1D->3s.3p 1P d=-0.49). So when the prior gets out of the way the
fit STILL won't move the CI -- because nothing in the objective pulls them. The
energy-level term (which dominates chi2) is satisfied by the complete basis
WITHOUT any CI shift; the handful of singlet gf lines carry negligible weight.

WHY the 9-config fit could move CI 59% but the full basis can't: at 9 configs the
freed CI stood in for ~113 MISSING configs, so the ENERGY LEVELS themselves had a
strong residual that only a large CI move could remove -- a strong signal forced
it. On the complete basis those configs are physically present, the energy levels
fit without a CI shift, and only the (weak, under-weighted) gf term wants the
move. Freeing CI is therefore a SMALL-BASIS CRUTCH, not a transferable lever: its
power came from energy-level signal that completeness removes.

STRATEGIC REDIRECT: the gf gap is NOT fixable through the parameter PRIOR (ridge,
freeze/free, scale) -- all of those have now been exhausted (raw HF, scaled HF,
freeze-CI, free-CI tight, free-CI loose: gf stays 0.21+/-0.01). The lever that
remains is the OBJECTIVE itself. The energy-only (+ token gf) objective doesn't
constrain the singlet 1D/1P eigenvector mixing, so no reparameterization the
energy fit chooses will fix those lines. Options, in order:
  1. UP-WEIGHT / ADD the singlet 1D/1P gf lines as explicit targets (raise lambda
     on just those lines, or add 3s.nd 1D->3s3p 1P + 3s2->3s4p 1P NIST gf to the
     objective). This directly injects the missing signal -- the Tier-2 thesis
     (gf-aware fit) applied surgically to the lines the energy fit can't see.
     With CI now freeable (this session's tooling), the fit HAS the knob to
     respond once the objective asks. THIS is the natural next experiment.
  2. Accept ab-initio gf for these lines (per-line model: report HF gf where the
     energy fit demonstrably can't improve it). Pragmatic fallback.
  3. Transcribe Bob's exact per-integral CI (5567 values) -- still available but
     now clearly LOWER priority: even perfect frozen CI is a fixed structure; the
     general engine needs the objective fix (1), which generalizes to any ion.

## CORRECTION (looked at the ENERGY LEVELS, not just gfRMS) -- the two misfits are the SAME physics

Prompted to check the level fits (I'd been quoting gfRMS in isolation and carrying
levelRMS as a benign scalar). KEY ERRORS FOUND in my own prior analysis:

(a) I'd been analyzing ING11.scaled (the PRE-FIT seed) for the low-block
    residuals, not ING11.scaledfit (the actual fitted model). The fit DOES move
    most low levels well: 3s3p -1844->+1, 3s3d -645->+48, 3s4p -158->-76 cm^-1
    (offset-removed). So "low levels are 1800 too deep" was a SEED artifact; the
    claim in the REDIRECT above that "the energy levels fit without a CI shift" is
    therefore HALF RIGHT -- true for 3s3p/3s3d, FALSE for 3p2.

(b) The real fitted-model (ING11.scaledfit) level residuals, honestly:
    median|d|=59, RMS=504, but the RMS is inflated by TWO distinct things:
    - HIGH-N SLOT-TRACKING ARTIFACT: top "outliers" are 3s10d 1D/3D at |d|~1900,
      Ecalc~55133 vs Eobs~60435 -- a rank-based slot grabbed the wrong near-
      degenerate high-Rydberg level (the 3D triplet all read identical -1901 =
      degenerate collision). These are ABOVE where the basis is reliable and
      should be cut with build_ine20 --max-energy (the cap exists, not applied to
      this fit). NOT a physics error; it pollutes both levelRMS and any line
      matching at high n.
    - REAL PHYSICS OFFENDERS (after cutting high-n): 3p2 1S/3P at -1360..-1650
      (frozen-CI over-repulsion, 3p2 EAV is free but can't climb out because the
      3s2-3p2 / 3d2-3p2 CI shove it down and ground 3s2 EAV is the pinned zero --
      classic surgical-CI signature), AND the SINGLET 3s.nd 1D / 3s.np 1P levels:
      3s3d 1D +573, 3s4d 1D +434, 3s4p 1P -432 cm^-1.

(c) THE PUNCHLINE: those singlet 1D/1P levels (3s3d 1D, 3s4d 1D, 3s4p 1P) that
    are 400-600 cm^-1 off in ENERGY are the EXACT SAME levels that are the
    endpoints of the gf lines 0.5 dex too weak. The energy misfit and the gf
    misfit are ONE problem -- the 1D/1P/3p2 mixing -- seen in two observables. The
    level residuals show the MECHANISM more clearly than the gf does: 3p2 sits
    -1400 cm^-1 deep, and 3p2 is the perturber that sets the singlet mixing. So
    the CI *was* the right lever; the reason free-CI didn't move it is that the
    energy term, though it HAS a 3p2 signal, is dominated by the ~190 well-fit
    high-n levels (and polluted by the high-n artifacts), so the 3p2 + singlet
    signal is a small fraction of chi2.

REVISED NEXT STEPS (supersedes the REDIRECT list above):
  0. FIRST clean the measurement: apply --max-energy (~IE, ~61671 cm^-1 for Mg I)
     so the high-n slot artifacts stop inflating levelRMS and corrupting matches.
     HONEST levelRMS (scaledfit, offset-removed, after cuts):
        all 99 levels      RMS 504  median|d| 59  max 1966 (3s10d slot artifact)
        Eobs<60000 (90)    RMS 294  median|d| 57  max 1382 (3p2 -- REAL)
        Eobs<55000 (29)    RMS 193  median|d| 73  max  616 (3s3d 1D -- REAL)
     So even cleaned, full-basis energy fit ~290 RMS / 60 median, dominated by 3p2
     and singlet 1D/1P -- vs 9-config 37 and Bob 12. The high-n cut helps (504->294)
     but the core misfit (3p2, singlet 1D/1P) is genuine, not an artifact.
  1. The 3p2 energy residual (-1400) is a CONCRETE, energy-side handle on the same
     CI that controls the gf. Two ways to use it: (i) FREE the 3p2-coupling CI AND
     up-weight 3p2's levels (its 4 levels are swamped by 190 Rydberg) so the
     energy term itself forces the CI move; or (ii) per-line gf up-weight as
     before. (i) is more principled -- it uses real energy data, not just gf.
  2-3 unchanged (per-line ab-initio fallback; Bob CI transcription last).

## TWO BIG FINDINGS: (1) gf-OFF beats gf-aware; (2) our LM harness < Cowan RCE (matching+Jacobian, not the LM core)

Prompted by two questions: "is the gf term helping?" and "is our LM worse than
Cowan's shipped RCE optimizer?". Ran a Bob-style PURE-ENERGY fit (lambda=0, gf
term OFF) + a head-to-head against RCE on a matched deck.

### Finding 1: turning the gf term OFF improves BOTH levels and gf.
Bob fits energy levels ONLY; the lambda*gf term is OUR Tier-2 addition. Compare
(full basis, scaled-HF bg, ridge=1, EAV+P free, IE-capped):
  combined  (lambda=3, gf-aware)  levelRMS 324  median|d| 59   gfRMS(A/B) 0.212
  PURE ENERGY (lambda=0)          levelRMS 242  median|d| 3.4  gfRMS(A/B) 0.174
median |d| 59 -> 3.4 cm^-1 (!!) and gf 0.212 -> 0.174 -- the BEST gf of any
full-basis fit, achieved with gf OUT of the objective. So the gf term was
actively HURTING: chasing ~20 gf lines pulled the eigenvectors away from the
energy optimum, degrading the bulk levels AND (net) the gf. Bob's instinct (fit
energies only, let good levels -> good eigenvectors -> good gf) is vindicated on
our engine. [This RETIRES the "up-weight the singlet gf lines" plan from the
prior section -- the gf term is the problem, not the cure. The remaining gf gap
(0.174 vs Bob ~0.13) is now a LEVELS/eigenvector problem, to be won on energies.]

### Finding 2: Cowan RCE converges cleanly and fits more levels; our LM core is
fine but the HARNESS around it (level matching + Jacobian) is the weak link.
Ran Cowan's actual RCE on a matched deck (build_ine20 --ruleset --max-energy
61671 from OUTGINE.scaled; 112 free = 40 EAV+39 G+32 ZETA+1 F; CI frozen). RCE
CONVERGED cleanly -> AVDEV (mean abs dev) = 117 cm^-1 over its 204 levels. No
divergence: the earlier RCE blow-ups were the (now-fixed) J-block bug + no
regularization, NOT the RCE algorithm.
  RCE:     204 levels, AVDEV(MAD) 117 cm^-1, converged in 5 iters (~0.02 min)
  our LM:   99 levels, MAD 185, RMS 242-923, but median|d| 3.4 cm^-1
KEY READ: our LM's OPTIMIZATION is NOT broken -- median 3.4 cm^-1 means it nails
the bulk. What's worse than RCE is the BLACK-BOX HARNESS:
  (a) SLOT-TRACKING BY ENERGY-RANK is fragile: it loses/swaps near-degenerate
      high-n levels -> manufactured ~2000 cm^-1 residuals (the 3s10d artifacts)
      that wreck MAD/RMS. Re-evaluating even gives RMS 242->923 as re-slotting
      catches more mismatches. RCE matches by internal EIGENVECTOR IDENTITY every
      iteration and never has this. <- DOMINANT deficit.
  (b) FINITE-DIFFERENCE Jacobian (diff_step 1e-3, ~9 evals/param/step): noisy
      (chi2 jitters +-5 around 2334), slow (2000 evals/50 min) vs RCE's analytic
      dE/dparam = <psi|dH/dparam|psi> (Hellmann-Feynman, ~free from eigenvectors),
      5 iters/1 sec.
  (c) the matcher SILENTLY DROPS levels it can't track -> we fit 99 of 204, a
      smaller/different problem than RCE.

ANSWER to "is our LM worse than RCE": YES, but the deficit is the matching +
Jacobian HARNESS, not the LM core. Fix = (1) replace rank-slot matching with
RCE-style eigenvector-identity matching re-done each iteration (kills the
artifacts, recovers the dropped levels); (2) optionally an analytic Hellmann-
Feynman Jacobian (faster, denoises). NOT "abandon LM for RCE" -- we still need
the Python optimizer for regularization (ridge) + UQ (covariance), which RCE
lacks; we need to bring its matching/Jacobian rigor INTO our optimizer.

ARTIFACTS: ING11.energyonly / OUTG11.energyonly (Bob-pure LM fit);
INE20.rce_compare / RCEOUT / PARVALS / rce_compare.log (RCE head-to-head);
work/fit_mg1_full_energyonly.py (lambda=0 driver, RUNDIR + max_energy + FREECI).
NB: concurrent fits in one run dir CLOBBER each other's OUTG11 -- run sequentially
or in separate dirs WITH the full RCG input deck (rwfn.dat alone is insufficient).

## IDENTITY MATCHING: levelRMS 242 -> 27 cm^-1 (the rank-slot harness WAS the deficit)

Replaced rank-slot level matching with RCE-style EIGENVECTOR-IDENTITY matching,
re-resolved every evaluation (gf_fit: setup_targets stores (cfgkey,termkey,Jstr)
-> E_obs; new _levels_by_identity() maps current OUTG11 levels by identity;
energy_resid/resid_vector look up by identity, not frozen rank). Re-ran the
pure-energy fit (lambda=0, scaled bg, EAV+P, IE-capped):

  rank-slot,   lambda=0   levelRMS 242   median|d| 3.4   gfRMS 0.174   N=99
  IDENTITY,    lambda=0   levelRMS  27.2 median|d| 1.8   gfRMS 0.163   N=98

levelRMS 242 -> 27 (9x). CONFIRMED: the old "242" was almost ALL rank-slot
artifacts (near-degenerate high-n levels swapping energy-rank -> bogus ~2000
cm^-1 residuals, e.g. 3s10d). With identity matching the worst residual is +175
(3s4d 1D), only 5 levels exceed 50 cm^-1, and 3p2 dropped from -1400 to -31. Our
LM now matches the 9-config fit (37) and approaches Bob (12) -- on energies.

ANSWER to the optimizer question, FINAL: our LM core was NEVER the problem. The
deficit vs Cowan RCE was entirely the rank-slot matching harness. Fixed, our
regularized LM gives 27 cm^-1 AND ridge + covariance (which RCE lacks). RCE's
AVDEV 117 is MAD over 204 levels incl. high-n; not directly comparable, but our
fit is now clearly excellent, not a laggard.

REMAINING gf-relevant residual is now CLEAN and POINTS AT THE SINGLET MIXING:
top-10 worst levels are the 3s.nd 1D series (3s4d/5d/3d/6d/7d 1D J=2, +60..+175
cm^-1) -- the SAME configs whose gf is off. With artifacts gone, the energy fit
and the gf gap are visibly ONE problem (the 3s.nd 1D / 3p2-3d2 mixing), now
diagnosable without artifact contamination. The 1D series sitting too HIGH in
energy is the lever to chase next (on energies, gf-off).

TWO LOOSE ENDS:
- COVERAGE: identity matching fits 98 of the ~204 NIST levels below IE; RCE fits
  204. The other ~106 NIST levels have no SEED computed level whose dominant
  eigenvector identity matches (cfgkey/termkey normalization mismatch, or the
  dominant component differs from the NIST label). Separate harness target:
  raise seed coverage (looser identity, or match by 2nd component / energy
  fallback) so we fit the same population RCE does.
- SPEED: loosening LM tol (1e-8 -> 1e-6, now TOL env) did NOT shorten the run
  much (still ~2868 evals / 60 min). Cause: cost is dominated by the FINITE-DIFF
  JACOBIAN (167 params x ~9 evals = ~1500 evals PER step), not the number of
  steps, so fewer steps barely helps. Real speed levers: (a) analytic
  Hellmann-Feynman Jacobian dE/dp = <psi|dH/dp|psi> (RCG has the pieces); (b)
  fewer free params (Bob-faithful ruleset = 112 not 167); (c) cheaper finite
  diff. tol still worth keeping loose for exploration but it's a minor lever.

## ADVERSARIAL ANALYSIS of the next step: the 1D residual is an UNOBSERVED PERTURBER, not a fittable CI

Ran a 4-proposal / 4-red-team / synthesize workflow (8 agents) on "how to fix the
remaining 1D-high/3D-low residual" instead of reflexively freeing CI. The red team
(verified against the actual ING11.energyonly artifact, not prose) overturned the
obvious answers:

RESIDUAL (current best, identity-matched energy-only): non-D fits to 8.7 cm^-1;
the misfit is almost entirely 3s.nd 1D HIGH (+175 4d, +90 5d, +77 3d, +69 6d, +62
7d; mean +94) and 3s.nd 3D LOW (mean -25). ~85% of the error is on the SINGLET.

WHY THE OBVIOUS LEVERS FAIL:
1. G^2(3s,nd) EXCHANGE is NOT the lever. It is ALREADY a free P-param and the fit
   already moved it to ~0.67-0.80*HF (NOT pinned by the sigma=0.30 ridge). The
   1D-3D splitting is too NARROW (3d calc 1450 vs NIST 1554; 4d 834 vs 1057 ->
   needs LARGER G^2) yet the free energy fit chose SMALLER G^2. Loosening its
   ridge lets it shrink further -> 1D WORSE. And a single exchange integral shifts
   1D/3D SYMMETRICALLY about the centroid; it cannot zero an ~85%-singlet-only,
   asymmetric residual. So G^2 cannot be the fix (refutes the "loosen/smoothness-
   prior the 3s.nd G^2" proposals).
2. FREEING the 1D-channel CI (3s.nd-3p^2/3d^2) is a near-verbatim Path-2 re-run.
   Path-2 already showed these CI move <4% (and LESS when the prior is loosened),
   because nothing in the energy objective pulls them. That failure is UPSTREAM of
   the gf-off + identity-matching fixes, so it would recur.

THE DECISIVE INSIGHT (verified against NIST): the 3p^2 1D perturber that drives the
3s.nd 1D series IS UNOBSERVED. NIST has 3p^2 3P (57813) and 3p^2 1S (68275) but NO
valence 3p^2 1D. Its position is pure MODEL OUTPUT, set by the very frozen CI we'd
be blaming. So freeing 3s.nd-3p^2 CI adds a near-FLAT direction with NO DATA
GRADIENT (the optimizer can place the unobserved 1D anywhere) -- which is mechan-
ically WHY Path-2's CI never moved, independent of objective. The "1D peaks at 4d
= near-resonance with 3p^2 1D" argument is CIRCULAR (the resonance position is a
free consequence of the parameters, not a datum).

REFRAME: the misfit is not a wrong INTEGRAL, it is a MISPLACED UNOBSERVED PERTURBER
(3p^2 1D, and likely 3d^2 1D). The physical lever is the PERTURBER's ENERGY, which
no Mg I level constrains.

RECOMMENDED NEXT STEP (survivor of the adversarial pass): PIN the unobserved
perturber position from external structure -- Bob's converged 3p^2 1D (from his
c1200*.log fitted model) or an MCHF value -- rather than fitting a CI the data
can't pull. Add 3p^2 1D (and 3d^2 1D) as a weakly-weighted PSEUDO-LEVEL target at
the external value, or pin its EAV. GENERALIZABLE RULE: "when a Rydberg term's
residual implicates an unobserved same-LS-symmetry perturber, the lever is the
perturber's ENERGY (fix it from a higher-fidelity model / external data), NOT the
series-perturber CI -- and the spin-partner term with no perturber (here 3s.nd 3D)
is the built-in control that calibrates the radial G^k/EAV."

CHEAP DE-RISK FIRST (pre-registered, refutable): a forward-model probe perturbing
G^2(3s,nd) AND the 3s.nd-3p^2/3d^2 CI, reading BOTH the 1D-3D split and the 1D-only
offset. CONFIRM perturber-route if the 1D excess is a perturber signature G^2 alone
can't remove (and moving CI/perturber-energy does move it); KILL it if a radial-
only (G^2+EAV) setting flattens the 1D excess after all. WHAT WOULD KILL THE WHOLE
DIRECTION: if the 1D excess turns out to be an identity/eigenvector-labeling
fragility at the 1D/3D J=2 boundary (a parse artifact), not a physical perturbation
-- worth ruling out since 3s4d 1D (+175) is exactly where 1D/3D mixing is strongest.

SECONDARY (cheap, frozen, no-objective-needed): transcribe Bob's actual per-integral
3s.nd-3p^2/3d^2 CI VALUES and 3p^2 1D position from c1200*.log and FREEZE them in --
tests "is it the CI VALUE / perturber position" without needing the objective to
pull anything. This is the one demonstrated route to Bob's quality and sidesteps
the no-data-gradient problem entirely.

HONEST CAVEAT: even pinning 3p^2 1D may not fully close it if 3d^2 1D also matters
or the CI MAGNITUDE (not just perturber position) is off. But it attacks the actual
degeneracy (an unobserved latent) instead of a knob the objective provably won't move.

### KILL-CHECK PASSED + mechanism CONFIRMED in the eigenvectors
Parsed raw LS eigenvector purities of the even-J=2 block (OUTG11.energyonly):
  3s.nd 1D levels: purity 0.71-0.99 (NOT 50/50) -> labeling is robust, the
    residual is NOT an identity-mislabel artifact (skeptic hypothesis RULED OUT).
  3s3d 1D: 0.78 self + 0.19 3p^2 1D;  3s4d 1D: 0.74 + 0.12 3p^2 1D  <- the
    low-n 1D members carry 12-19% of the UNOBSERVED 3p^2 1D perturber.
  EVERY 3s.nd 3D: purity 0.99-1.00, ZERO 3p^2/3d^2 admixture (3p^2/3d^2 have no
    3D term) -> the clean spin-partner control, exactly as predicted.
So the 1D-high residual IS the 3p^2 1D perturber mixing (12-19%), and that
perturber is unobserved -> its energy/CI is a data-free latent. Diagnosis airtight;
the lever is the perturber position, not a freeable CI. Proceed to pin 3p^2 1D
(and 3d^2 1D) from external structure (Bob c1200*.log / MCHF).

### CORRECTION (user caught it): the lever IS a free knob we STARVED -- F^2(3p^2), via the energy cap
The adversarial synthesis concluded "3p^2 1D is unobserved -> pin it externally."
That was OVER-complicated. The user asked: isn't there a knob we're not fitting
that directly affects these levels, with priors? YES -- F^2(3p^2) (the 3p2|P2
direct Slater). It is ALREADY in the free set (3p^2 has observed levels) and it
sets the 1S/1D/3P splitting WITHIN 3p^2 -- so fitting it to the OBSERVED 3p^2 3P
and 1S DETERMINES the (unobserved) 3p^2 1D position. There IS a data gradient; it
routes through the observed 1S-3P splitting of the same config.

WHY F^2 WAS STUCK (root cause, verified): the model's 3p^2 1S sits at 54358, but
NIST 3p^2 1S = 68275 -- ABOVE the ionization limit (61671) and thus ABOVE our
max_energy cap. The cap EXCLUDED the one observed level that pins the 1S-3P
splitting. With only 3p^2 3P (a near-centroid) in the fit, F^2 had NO splitting
signal and sat at HF (0.025). Model 1S-3P = -219 (1S BELOW 3P!) vs NIST +10462 --
qualitatively wrong, F^2 ~80x too small. Forward-model probe confirms F^2 is the
lever with the right sign: scanning F^2 0.025->2.0 moves 1S-3P -219 -> +5633
(marching toward +10462). So we had a free knob with strong observed-data signal
and were throwing away its data with the energy cap.

FIX (principled + generalizable): exempt from the max_energy cap any observed
level whose config ALREADY has a sub-cap observed level -- i.e. a valence config
we're already fitting; its above-cap term completes the term structure that pins
its Slater integrals. High-Rydberg configs (no sub-cap anchor) stay capped.
Surgical: admits EXACTLY 3p^2 1S (+1 level, 98->100 w/ identity), no Rydberg flood.
This is the RIGHT general rule: keep term-splitting data that constrains radial
integrals, drop only the unanchored Rydberg tail. The cap was a blunt instrument.

LESSON for the engine + the adversarial method: the workflow correctly killed the
G^2 and free-CI routes, but its "unobserved perturber -> external pin" conclusion
missed that the perturber's position is FIXED BY F^2 which IS constrained by other
OBSERVED levels of the same config -- which we had CAPPED AWAY. The signal wasn't
absent; we discarded it. (Re-fit running: tests whether restoring F^2's signal
fixes 3p^2 1S/1D AND propagates through the 12-19% eigenvector mixing to pull down
the 3s.nd 1D residual.)

### BOB'S VALUE + a SLOT-SLIP CORRECTION to the F^2 probe
What Bob set F^2(3p^2) to (from hf1200z.dat line 176 + c1200ez.log param 141):
  Bob HF F^2(3p^2) = 23179.06 cm^-1;  Bob FITTED = 13430.2 (ratio 0.579);
  status = FREELY FITTED (step 40, NOT FIXEDHF). So Bob DOES free this exact knob
  and pulls it to 0.58*HF -- direct validation of the instinct to fit it.

CORRECTION to my earlier probe: F^2(3p^2) is our 3p2|P0 slot (raw-HF ~21.1 kK,
matching Bob's HF 23179 within the ~9% HF-code diff), NOT P2 (=0.03, a junk/
residual slot I mistakenly scanned). The "1S-3P = -219, wrong sign, F^2 80x too
small" claim was the WRONG SLOT -- DISREGARD it. Re-probing the correct slot (P0):
  F^2 x1.00 (HF):  1S-3P = +11688   (NIST +10462 -> only ~12% too wide)
  F^2 x0.85:       1S-3P = +10021   <- ~nails NIST
  F^2 x0.579(Bob): 1S-3P = +7004    (too small in OUR basis; Bob's ratio isn't
                                     transferable -- his HF/CI differ from ours)
So the 3p^2 block is NOT qualitatively broken; the 1S-3P splitting is ~right at HF,
~12% too wide, wants F^2 ~0.85*HF. The earlier "10000 cm^-1 / wrong-side" panic was
the slot slip.

NET (honest): the instinct holds -- F^2(3p^2) IS a free, data-constrained knob and
the cap WAS starving it of the 3p^2 1S signal (fix admits exactly that level). But
the correction is MODEST (~12% on the splitting), so fixing F^2 will tidy the 3p^2
block and nudge the 3p^2 1D perturber DOWN a bit -- it is unlikely to fully erase
the +94 mean 3s.nd 1D residual by itself. The running re-fit (with 3p^2 1S admitted,
F^2=P0 free) tests how much of the 1D residual it actually removes.

## ROBUST DIAGNOSTICS: tools/param_labels.py + tools/diagnose.py (kill the comparison-error class)
Several level/parameter comparison mistakes this session (wrong eigenvector slot
for F^2(3p2): scanned P2=0.03 not P0=21kK; ground-offset misread; silent empty
_cfgkey matches) all shared a root cause: hand-rolled comparisons in throwaway
snippets that re-derived slot indices / unit factors / identity keys and silently
returned plausible-but-wrong numbers. Fixed at the source with TWO validated tools:

tools/param_labels.py -- AUTHORITATIVE physical labels for ING11 Slater slots.
RCG echoes the physical integral NAMES in OUTGINE block headers ('EAV 3p2 F2(22)
ALPHA ZETA 2'); we read that order (reuse build_ine20._parse_param_names) and join
positionally to the ING11 raw slots, with HARD asserts: >names-than-slots, an
adjustable slot beyond the names, or a unit/value out of range all RAISE instead
of mislabeling. param_value_cm1(ing11,outgine,cfg,'F2') returns the value in cm^-1
for the CORRECT slot or raises. Verified: F2(3p2)=P0=21095 cm^-1 (Bob HF 23179,
~9% code diff), scaled matches RCG's own PARVALS echo (16876 vs 16885.6) to <10.
KNOWN LIMITS (declined, not guessed): open-shell configs whose params overflow the
4 ING11 slots (3d2: F2 F4 ALPHA BETA T(D2) ZETA) and odd-block configs the header
walk hasn't reached (3s3p) are SKIPPED with a recorded reason -- needs the
extended-slot + even/odd-block reader for Fe II later.

tools/diagnose.py -- ONE validated path for the comparisons (reuses the SAME
identity matching the fit uses + param_labels). Subcommands:
  levels    -- per-level residuals, offset-removed (median|ground anchor, the fit's
               convention), worst offenders + per-term means. RAISES on zero matches.
  level     -- one level + its EIGENVECTOR PURITY and top components (the 1D/3D
               kill-check in one command: 3s4d 1D purity 0.88, 0.067 of 3p2 1D).
  splitting -- computed term splitting vs observed (3p2 1S-3P = +10441 vs +10462).
  param     -- physical integral value vs a reference (F2(3p2)=21095, 0.91x Bob).
Use these (import or CLI) instead of ad-hoc snippets for any level/param comparison.

ALSO FOUND (via diagnose levels) -- the cap-exemption rule is TOO BROAD: admitting
'any config with a sub-cap level' re-admitted high-n 3s.10d/11d Rydberg levels
(they have sub-cap members 3s3d etc.) at -14000 cm^-1 slot artifacts, wrecking the
3p2-1S re-fit (RMS 27->68, gf 0.163->0.313). FIX NEEDED: exempt only NON-RYDBERG
valence configs (e.g. 3p2/3d2), not Rydberg series. The F^2(3p2) fix itself WORKED
(3p2 1S-3P splitting now +10441 vs obs +10462, off 21 -- was +11688 at HF).

## RESULT: F^2(3p^2) FIXES the 3p^2 block but does NOT fix the 3s.nd 1D (it makes it worse)
Re-fit energy-only with 3p^2 1S admitted (cap exemption for low-n valence configs,
max_n<=4) and F^2(3p^2) free. Validated via diagnose.py / param_labels.py:

WHAT WORKED (the user's instinct, confirmed):
- F^2(3p^2) MOVED: 21095 (HF) -> 17220 (0.743*HF), pulled toward Bob's 13430 (0.579).
- 3p^2 1S-3P splitting NAILED: +10441 vs observed +10462 (off -21 cm^-1); was
  +11688 at HF. So admitting the 3p^2 1S level gave F^2 its signal and fixed the
  3p^2 internal term structure exactly as predicted.

WHAT DIDN'T (decisive negative on the 1D chain):
- 3s.4d 1D residual got WORSE: +175 -> +434 cm^-1. Its 3p^2 1D admixture is
  unchanged (~6.7%). So relocating the 3p^2 1D perturber (by fixing F^2) did NOT
  pull the 3s.nd 1D series down -- it pushed it UP. The mixing FRACTION is set by
  the CI (frozen), not by the perturber position; moving 3p^2 1D changed the
  energy denominator the wrong way for the 3s.nd 1D members.
=> Fixing the perturber's internal structure is necessary for the 3p^2 block to be
   right, but it does NOT resolve the 3s.nd 1D residual. The 1D problem really is
   the frozen CI MIXING (Path-2 territory), which remains data-gradient-starved.

THE 27->68 RMS REGRESSION was NOT the cap rule (my earlier theory was wrong). It is
the 3s.10d/11d levels: E_obs ~60435-60735 (BELOW the 61671 IE cap) but BEYOND the
122-config basis's reliable reach, so E_calc is ~14000 cm^-1 low -- genuine basis-
limit misfit, not a slot artifact or a cap-exemption bug. RMS vs cap:
  cap 60000 -> 87 levels, RMS 70.5, median|d| 1.8 (clean, 1D-dominated)
  cap 61671 -> 99 levels, RMS 4045  (the 3s.10d/11d basis-limit tail explodes)
FIX: set max_energy ~60000 (drop the basis-unreachable top Rydberg members), OR
grow the basis. With cap 60000 the honest fit is median|d| 1.8, RMS 70.5, the RMS
dominated by the singlet 3s.4d 1D (+433) -- the CI-mixing residual, unchanged.

NET STANDING: best honest model is still ~median 1.8 cm^-1; the headline RMS depends
on the cap (basis reach). The 1D/gf gap is confirmed to be frozen-CI mixing that
neither G^2, nor free-CI (no data gradient), nor F^2(3p^2)/perturber-position can
fix. The remaining principled levers: (a) transcribe Bob's curated 3s.nd-3p^2/3d^2
CI VALUES and FREEZE them (no objective pull needed); (b) accept ab-initio gf for
the few singlet 1D->1P lines (per-line fallback). Tooling this session
(param_labels.py, diagnose.py, the max_n cap rule) is the durable win.

## PARAMETER COMPARISON vs BOB (tools/compare_bob.py)
Joined our fitted params (param_labels, validated slots) to Bob's curated c1200*.log
values by (config, physical-name), cm^-1. EAVs de-offset by the global energy-zero
constant (-3133 cm^-1, itself a striking confirmation: nearly EVERY EAV is ours =
Bob - 3133, ratio 0.94-0.95).

EAVs (config centroids) -- EXCELLENT where the basis is reliable:
  23 configs match to median |diff| 5.4 cm^-1, mean 36; the high-n Rydberg series
  (3s.ng, 3s.ni) agree to <1 cm^-1. Our energy fit recovers Bob's config positions
  almost exactly. Outliers (all understood): 3s.10d/11d (-14000, basis-limit tail);
  3p2 (+840, shifted by the big frozen 3s2-3p2 CI=22).

STRUCTURE integrals (F/G/ZETA) -- where we diverge, but only 2 matched so far:
  3p2 F^2: ours 17220 vs Bob 13430 (1.28x) -- we are 28% HIGH. Bob fit it to
    0.58*HF; we reached only 0.74*HF. (Connects to our 3s.nd 1D residual being
    WORSE than Bob's: his lower F^2 places the 3p^2 1D perturber differently.)
  3p2 ZETA: ours 25 vs Bob 37 (0.67x) -- we under-fit the spin-orbit by a third.

CAVEAT (limits the structure comparison): only 2 F/G/ZETA matched because most of
Bob's freely-fit structure integrals (his 3s.np G^1 exchange series, the odd-block
F/G) live in ODD-PARITY-block configs that param_labels can't yet read (the
even/odd-block + open-shell limitation noted in the tooling section). To compare
the structure integrals comprehensively (the real test of how close we get to
Bob), param_labels needs the even/odd block walk -- worth doing before Fe II anyway.

HEADLINE: on the parameters BOTH tools can compare, our energy fit matches Bob's
EAVs to a few cm^-1 (excellent); the structure integrals we CAN see (3p2 F^2/ZETA)
are 20-30% off, in the direction consistent with our worse singlet 1D fit.

## FULL parameter comparison vs Bob (param_labels now reads BOTH parity blocks)
Extended param_labels._cfg_param_names to parse ALL OUTGINE header regions (the
122-config OUTGINE has TWO: even block near the top, odd block ~line 1993,
separated by matrix data; previously only the even block was read). Coverage:
~60 -> all 122 configs; F2(3p2) unchanged (no regression). Structure-integral
matches vs Bob: 2 -> 18.

RESULT (Bob's freely-fit params, ours from validated slots, cm^-1):
  EAV: 51 configs, median |diff| 5.1 cm^-1 (excellent, comprehensive).
  STRUCTURE INTEGRALS -- a clear SYSTEMATIC BIAS, not random scatter:
    G^1 EXCHANGE (3s.np series, n=3..12): mean ours/Bob = 1.17 (range 0.96-1.27).
      n=3 (3s3p G1) agrees (0.96) but every higher member is ~15-27% HIGH.
    ZETA SPIN-ORBIT (3s.np): mean ours/Bob = 0.63 -- consistently ~37% LOW.
    F^2(3p2): 1.28 (28% high).
  => Our energy fit systematically OVER-estimates exchange G^k (~+17%) and UNDER-
     estimates spin-orbit ZETA (~-37%) across the Rydberg series, while nailing the
     EAVs. This is a BIASED-LEVER signature: with EAVs free and excellent, the fit
     compensates structure errors that a per-family screening/prior would fix
     uniformly. G^1 exchange sets the singlet-triplet splittings that shape the
     eigenvectors -> directly connected to the residual gf gap (our too-large G^1
     => wrong 1S-3P-type splittings => degraded singlet gf).

WHY (likely): (a) our HF G^k start ~6% off Bob's HF (different HF codes, noted) but
that's small; the ~17% is a FIT bias. (b) Bob FREES G^1 per member AND his frozen
CI background is scaled (0.8); our frozen raw/loosely-scaled CI leaves the exchange
to absorb mixing error -> biased high. (c) ZETA low likely because our ridge centre
(raw HF) + sigma pulls it down, or the fine-structure signal is weak per member.

ACTIONABLE: the systematic, same-sign, similar-magnitude bias across a whole family
is exactly what a PER-FAMILY prior centre (Bob's FIXEDHF scales: G^1~0.66, ZETA~1.0,
CI~0.8) would correct -- the generalizable "scaled-HF prior centre per family" rule
(ruleset Rule 2), now with Bob's actual numbers as targets. tools/compare_bob.py is
the instrument to measure convergence toward them. NEXT candidate: set the ridge
PRIOR CENTRE (not just seed) to Bob's per-family scaled-HF and re-fit; watch the
G^1/ZETA ratios move toward 1.0 and the singlet gf improve.

## WHY THE G1/ZETA BIAS PERSISTS: wrong prior CENTRE + no data signal (NOT a too-tight prior)
Checked how far the fit moved each P-param from its ridge prior centre (sigma_P=0.30
fractional). The params sit AT the centre, not straining against it:
  3s.np ZETA: moved 0.00 sigma (every member) -- did not move at all.
  3s.np G^1:  moved 0.04 / 0.53 / 0.12 / -0.01 sigma -- tiny, far inside the band.
A TOO-TIGHT prior would park a param at the sigma BOUNDARY fighting a data pull.
These are at the CENTRE => no data pull; they rest wherever the centre is. So the
bias is NOT prior-width; it is a WRONG PRIOR CENTRE that the (weak per-member) data
cannot correct.

The centre is scaled-HF, and it is itself biased:
  ZETA(3s3p): centre 22.5 vs Bob HF 39.0, Bob FITTED 41.0 -> centre = 0.55*Bob_fit.
    (Our scaled-HF halves zeta; fine-structure signal per member too weak to fix it
    -> zeta comes out ~37% low, stuck at the bad centre.)
  G^1(3s3p): Bob HF 29818, our centre 22308 = 0.75*ourHF ~ 0.95*Bob_fit (good at n=3).
  G^1(3s4p): our centre 2483 = 1.09*Bob_fit -- ALREADY above Bob for n>=4. The FLAT
    0.75 screening that fits n=3 over-shoots higher members: the right screening is
    n-DEPENDENT, ours is constant.

CONCLUSION: loosening sigma will NOT help (same null as the CI loosen test -- no data
gradient to exploit; a looser prior just adds variance around the same wrong centre).
The fix is the prior CENTRE, per family, matching the n-dependence:
  - set ZETA centre to ~Bob's (unscaled HF is closer: Bob fit 41 vs HF 39, scale ~1.0,
    NOT 0.55) -- our zeta HF-screening is simply wrong for spin-orbit (should be ~1.0).
  - set G^1 centre via Bob's actual per-member FIXEDHF scale (decays from ~0.785 at
    n=3), not a flat 0.75.
i.e. scale_hf.py must use PER-FAMILY (and ideally per-member) screening for the
PRIOR CENTRE: ZETA~1.0, G^1~0.66-0.785(n), F^2~0.8, CI~0.8. compare_bob.py measures
convergence. This is the concrete, generalizable version of ruleset Rule 2.

## CORRECTION: the ZETA bias IS the ridge blocking a real signal (not "no signal") -- quantified
Pressed on "why won't the fitter move zeta to the right value?" -- and my earlier
"weak/no data signal" claim was WRONG for zeta. Forward-model probe:
  3s3p 3P J2-J0 splitting: model@zeta=22.5 -> 34.0 cm^-1; obs 60.8. Setting zeta=41
  (Bob's) -> 61.0 = obs. The 3P J=0,1,2 levels ARE in the fit. So zeta has a STRONG,
  clean, OBSERVED gradient (~1.3 cm^-1 splitting per cm^-1 zeta) and moving it 22.5->41
  fixes the fine structure almost exactly. The fitter SHOULD move it.

WHY IT DOESN'T -- the ridge penalty outweighs the energy gain (ridge=1):
  energy chi2 GAIN from fixing the 27 cm^-1 fine-structure error = (27/SIGMA_E)^2
    = (27/50)^2 = 0.29  (the residual is only 0.5*SIGMA_E -- below the noise floor).
  ridge COST of the move (zeta 22.5->41 = 0.82 fractional = 2.7 sigma_P at sigma=0.30)
    = (0.82/0.30)^2 = 7.5.
  Ridge penalty 7.5 >> energy gain 0.29 (26x). The regularized objective is MINIMIZED
  by leaving zeta WRONG. The optimizer is correct; the PRIOR is mis-specified.

So this IS a too-tight/mis-centred prior (contra the CI case, where signal was truly
absent and loosening did nothing). Here the signal is strong; the prior just out-votes
it. TWO independent fixes, EITHER works for zeta:
  (a) fix the CENTRE: zeta prior centre should be ~unscaled HF (~39-41), not the
      0.55-scaled 22.5. Our scale_hf zeta screening (~0.55) is simply wrong; spin-
      orbit should NOT be screened down (Bob's scale ~1.0). With the right centre, no
      move needed.
  (b) loosen sigma on zeta (and exchange): SIGMA_E=50 makes a 27 cm^-1 fine-structure
      error "free", so the ridge dominates. Either lower SIGMA_E (care about fine
      structure) or raise sigma_P for zeta/G so a 2-3 sigma move toward truth is cheap.
GENERAL LESSON for the engine: with ridge=1, sigma_P=0.30, SIGMA_E=50, any parameter
whose data signal moves the levels by < ~SIGMA_E and which sits > ~1 sigma from its
prior centre is FROZEN by the prior, even with a clean gradient. The whole G^1/ZETA
bias family is this: real signal, out-voted by a prior centred on a biased scaled-HF.
The earlier "G^1/ZETA stuck because no signal" note is SUPERSEDED -- it is the prior
(centre + width vs SIGMA_E balance), and it IS fixable.

## RESULT: per-family prior CENTRE fix (zeta un-screened, G 0.66, F 0.8, CI 0.8)
Made scale_hf.py NAME-AWARE (via param_labels) so each Slater slot gets its family
scale: ZETA x1.0 (was wrongly x0.8 -> the bias), F x0.8, G x0.66, CI x0.8. Regenerated
ING11.scaled (the prior CENTRE) and re-fit energy-only.

ENERGY: chi2 178 -> 41 (4.3x); levelRMS 30.8, median|d| 2.0 cm^-1. Much better basin.

PARAMETERS snapped toward Bob (compare_bob, ours/Bob mean):
  G^1 exchange: 1.17 -> 0.99 (NAILED)
  ZETA spin-orb: 0.63 -> 0.79 (big improvement, still 21% low)
  F^2 direct:   1.28 -> 1.16
The G^1/ZETA systematic bias is largely CORRECTED -- the prior-centre diagnosis was
right and the fix works.

FINE STRUCTURE: 3s3p 3P J2-J0 = 34 -> 42 cm^-1 (obs 60.8). Improved but NOT fixed,
because the fit moved zeta only to 28.1 (0.69*Bob's 41), not 41. The ridge STILL
partially holds zeta: un-screening moved the centre to ~31.5, but Bob's 41 is still
~1 sigma above even that, so the energy-gain-vs-ridge-cost tradeoff still under-shoots.
Releasing zeta fully would need EITHER a looser zeta sigma OR a higher zeta scale
(our HF zeta is itself below Bob's HF -- un-screening alone can't reach 41).

gf: strong A/B gfRMS 0.163 -> 0.190 (slightly WORSE). The worst lines are the SAME
singlet ones (3s2->3s4p 1P d=-0.57, 3s.5d 1D->3s3p 1P d=-0.50) -- the frozen-CI
mixing gap, untouched by radial params (as established). F^2 still 1.16x high may
have nudged the singlet eigenvectors unfavorably. So: the per-family prior is a clear
ENGINE win (chi2 4x, params converge to Bob, biases corrected) but does NOT crack the
singlet gf -- confirming again that gap is frozen-CI mixing, not radial structure.

NET: tools/scale_hf.py (name-aware per-family) + the prior-CENTRE insight is the
generalizable result. The remaining levers for the gf gap are unchanged: freeze Bob's
curated 3s.nd-3p^2/3d^2 CI values, or per-line ab-initio gf for the few singlet lines.
NEXT TUNING (cheap): raise zeta scale / loosen zeta sigma to finish the fine structure;
nudge F^2 toward Bob (0.8 centre). But these polish energies, not the gf gap.

## STALL ISOLATED: Bob's 3s2-3p2 CI can't be transplanted in isolation (couples to his pinned zero)
The Bob-CI + zeta fit STALLED (chi2 oscillated 250-360, spiked to 359208, never
reached the prior run's 41). User flagged it wasn't converging. Isolation:
  TEST B (Bob CI only, NO zeta release): STALLS IDENTICALLY -- same 358->284->...
    oscillation, same 359208 spike at eval 8.
  => BOB'S CI is the culprit, NOT the zeta release (which was a red herring;
     the surgical zeta fix was correct but irrelevant to this stall).

MECHANISM (confirmed, one RCG each): Bob's 3s2-3p2 CI = 25972 vs our scaled 23524
(+10%) drags the ground 3s2 from E_calc -3313 -> -3868, a -555 cm^-1 shift. But the
ground 3s2 EAV is the PINNED ENERGY ZERO (not free). So the fit cannot absorb the
shift by moving the ground; it must compensate through the global offset + every
other EAV shifting ~555 together -- a near-degenerate collective direction ("shift
all vs shift ground"). That stiff valley is what makes LM overshoot to 359208 and
oscillate. Ill-conditioned, not wrong-physics.

KEY LESSON: you CANNOT transplant a single one of Bob's CI values in isolation. His
CI values are self-consistent with HIS EAVs and HIS energy zero; we copied ONE CI
but kept OUR pinned ground, so the stronger CI just mis-places the ground 555 cm^-1
and the energy fit thrashes absorbing it. To use Bob's CI we'd need his CONSISTENT
set (CI + EAVs together), or to free the ground reference, or to re-anchor.

FIX OPTIONS:
  1. REVERT to scaled-HF CI (the per-family value 23524) -- back to the clean
     chi2=41 fit. Bob's exact CI isn't usable piecemeal. Simplest; loses the
     ground-residual test (which was the point of trying Bob's CI).
  2. Co-adjust: when setting Bob's CI, also shift the global energy reference /
     re-anchor so the ground lands right -- removes the degenerate direction.
  3. FREE the ground-config EAV (currently pinned) so the fit can place the ground
     against Bob's CI -- but then the energy zero floats (needs re-anchoring on a
     different reference). More invasive.
  4. Damp the offset direction explicitly (the real Tier-2 fix): the "shift-all"
     mode is a known null-ish direction; a small prior/constraint on it would
     stabilize LM regardless of CI. Generalizable.
Recommended: (1) to recover a clean fit now, then (4) as the principled engine fix
if we want Bob's CI later. The ground-residual question (does Bob's CI close the
-70?) is unanswerable in isolation -- it's entangled with the pinned zero.

## TWO-SCALE (s_G, s_F) screening: our deck wants s_G~0.8, NOT Bob's 0.6
User asked: can we have two fitted scale params, one for G (exchange) one for F
(direct)? Bob's 101 scaled FIXEDHF structure integrals use exactly two factors:
G=0.6, F=0.8. Tested whether those transfer to our deck.

Seed-RMS sweep was MISLEADING (monotonic in s_G -- it's offset-dominated, which the
EAV fit removes). The TRUE diagnostic is a STRUCTURE splitting the EAV fit can't
touch: 3s3p 1P-3P (set by G^1 exchange; obs 13201 cm^-1):
  s_G=0.6 -> 9349 (too small);  0.7 -> 11135;  0.8 -> 12916 (~obs);  0.9 -> 14693
  => optimal s_G ~ 0.81 for OUR deck (interp to 13201), NOT Bob's 0.6.
REASON: our HF G^1(3s3p)=27885 vs Bob's HF 29818 (different HF codes). Bob scales
HIS HF by 0.6; we need ~0.8 on OURS to reach the SAME physical value. The SCALE
FACTOR is not transferable -- only the target physical value is. F similarly ~0.8.

IMPLICATION: two global scales (s_G~0.8, s_F~0.8) is the right STRUCTURE, but the
values must be FIT to our deck, not copied from Bob. The splitting data gives a
strong, clean gradient on s_G (unlike the offset-confounded seed RMS), so fitting
s_G/s_F as 2 free params should converge well. This generalizes: every ion's HF
differs, so the screening scales must be fit per-ion, not transplanted.
NEXT: wire s_G, s_F as 2 fitted global multipliers on the frozen G/F integrals.

## LEVEL COUNT vs BOB: we fit ~the same (~85-106), but a different top-end
Bob (c1200{e,o}z.log): 85 observed levels (62 even + 23 odd), up to ~82000 cm^-1
  -- ABOVE the IE (61671), including autoionizing valence (3p2 1S @68275 etc).
Us (cap 60000): 89-106 levels. NIST has 229 <= IE, 106 <= 60000.

So count is comparable; the difference is WHICH top-end levels:
  (a) HIGH-N RYDBERG (3s.10d/11d, 60000-61671): we EXCLUDE (basis-limit -- our
      122-config basis places them ~14000 cm^-1 off). Bob's basis reaches them.
      Fitting these needs a BIGGER BASIS, not a higher cap (else -14000 artifacts).
  (b) ABOVE-IE VALENCE (3p2 1S @68275, low-l, reliable, real structure signal):
      we SHOULD fit these -- the cap-exemption (max_n<=4) admits them, and 3p2 1S
      pins F^2(3p2). Bob fits them.
ACTION: keep the max_n<=4 cap exemption (admits above-IE valence like 3p2 1S);
the high-n Rydberg gap is a basis-completeness issue (enlarge basis to match Bob's
top end), separate from the cap. We are NOT under-fitting in count -- we're capping
the basis-unreachable tail, which is correct, plus admitting the valence top-end.

## TWO FITTED G/F SCALES: nails the exchange splitting, but 2 globals are too rigid for full RMS
Implemented s_G (exchange) + s_F (direct) as 2 FITTED global multipliers on the
FROZEN G/F integrals (Forward fit_gf_scales=True; free EAV only). Combined best-model
fit: 76 params (73 EAV + 1 freed 3s2-3p2 CI + 2 scales), cap IE+exemption, parallel.

RESULTS (cleaner+faster: 291s, converged):
  s_G -> 0.669, s_F -> 0.741  (FIT to our deck; ~Bob's 0.6/0.8! data chose Bob-like
    screening when balancing all levels -- vs the splitting-only probe's 0.8.)
  WINS:
    3s3p 1P-3P splitting = +13191 vs obs +13201 (off 10!) -- the fitted G-scale +
      zeta NAILED the exchange splitting (was 9349 at s_G=0.6). Your 2-scale idea works.
    ground 3s2 residual = -12 (freed CI holds it).
  REGRESSION:
    levelRMS ~62-70 (basis-reliable caps) vs the 168-param free-all-P model's ~31;
    median|d| 7-10 vs 2-3.

CONCLUSION: 2 global scales fix the SYSTEMATIC screening (exchange splitting -- the
thing 2 scales CAN capture) but are too RIGID for the full level set; they can't
capture the per-config variation that 94 individual free G/F integrals did. Bob gets
away with 2 scales ONLY because he ALSO hand-curates 101 integrals individually; we
have just the 2 globals. So:
  free-all-P:  168 params, RMS 31, but ill-conditioned/thrashy (the 94 free G/F/zeta).
  2-scales:     76 params, RMS 62, clean+fast but underfit (too rigid).
RIGHT ANSWER is BETWEEN: Bob's recipe = free EAV + 2 scales + free a CURATED handful
of the most-varying individual G/F (not all 94, not zero). i.e. the ruleset: free
EAV + the key low-l exchange (3s.np G^1) + zeta, scale the REST with s_G/s_F, freeze
CI (+ free the few that matter like 3s2-3p2). That hybrid should get both clean
conditioning AND low RMS. NEXT: free EAV + 2 scales + ~the 3s.np G^1 series
individually (Bob frees exactly these), rest scaled.
NB: cap at IE (61671) re-admits 3s.10d/11d basis artifacts (-11665); use cap ~58000
for honest RMS, OR enlarge basis. The max_n exemption correctly adds low-l valence.

## HYBRID (EAV + 2 scales + individual 3s.np G^1 + freed CI): same wall -- it's the singlet CI
Ran the Bob-faithful hybrid: 80 params (73 EAV + 2 G/F scales + 4 individual 3s.np
G^1 + 1 freed 3s2-3p2 CI), surgical zeta, cap 58000, parallel. The wiring works:
  s_G fit -> 0.66 (Bob-like); the 4 G^1 moved individually (3s3p 18404->22549 etc);
  3s3p 1P-3P splitting = +13191 vs obs +13201 (off 10, NAILED); ground 3s2 ~-12.
  (The FIT-line gfRMS 0.663 was a transient bad-trial-point artifact; the real
   converged gfRMS is ~0.23.)
BUT: chi2 ~206, levelRMS ~62-70 (cap 58-59k) -- ESSENTIALLY THE SAME as the pure-
2-scale model, NOT better. Adding individual G^1 barely helped.

WHY: the worst levels AND gf are the SINGLET 3s.nd 1D series (3s.4d 1D->3s.3p 1P
d=-0.59, 3s2->3s.4p 1P d=-0.57) -- the SAME frozen-CI-mixing problem from the whole
investigation (unobserved 3p2 1D perturber). Exchange G^1 / scales / zeta fix the
TRIPLET + exchange structure (splitting nailed) but DO NOT touch the singlet 1D,
which is set by the frozen 3s.nd-3p2/3d2 CI. The free-all-P model reached RMS 31
NOT by fixing the singlet but by using its 94 free params as a CRUTCH to absorb the
error elsewhere (worse-conditioned, less physical).

CONVERGENT CONCLUSION (reached now from the clean Bob-faithful direction): the
remaining Mg I error -- both levels and gf -- is the SINGLET 3s.nd 1D / frozen-CI
mixing, and NO amount of EAV/scale/G^1/zeta freedom touches it. The ONLY levers
left are the ones long identified: (a) free the 3s.nd-3p2/3d2 CI (but data-gradient-
starved -- unobserved perturber), (b) pin the 3p2 1D / 3d2 1D perturber energy from
external data, (c) accept ab-initio gf for those few singlet lines, (d) bigger basis.
The Bob-faithful parameter recipe is now COMPLETE and clean (80 params, splitting
nailed, ground fixed, s_G/s_F~Bob) -- it just can't cross the singlet-CI wall, which
is a STRUCTURE/DATA problem, not a free-set problem.

## WHY BOB'S RMS BEATS OURS (12 vs 62): he fits 213 levels, we fit 50 -- DATA STARVATION
Compared Bob's per-level residuals (E-O column, c1200*.log) to ours directly.
  Bob: 213 distinct levels fit, E 0-98899 cm^-1, RMS(E-O) 11.7, MEDIAN 0.3 cm^-1.
       His worst are the SAME singlet 3s.nd 1D (+50 to +77) -- he has our problem too,
       just milder.
  Us:  50 levels (cap 58000), bulk(non-1D) median 9.5 / RMS 37, 1D RMS 245.

TWO gaps, the BULK one is bigger and was missed before:
  (1) BULK (non-1D): Bob median 0.3 vs us 9.5 -- ~30x worse. Our worst non-1D are
      3s.4d 3D -98, 3s.5d 3D -62, 3s.5p 1P -80 -- EASY levels Bob nails to ~0.
  (2) SINGLET 1D: Bob ~50-77 vs us 245 -- ~5x (the frozen-CI wall, real but smaller).

ROOT CAUSE of the bulk gap = DATA STARVATION. Bob fits 213 levels (up to 99000,
well above IE); we fit 50 (cap 58000). With 4x more levels constraining the same
radial params, his fit is FAR better determined; ours is under-constrained -> sloppy
bulk. We cap at 58000 because higher levels (3s.10d/11d+) come out as basis-limit
artifacts (-14000) our 122-config basis can't place. So:
  - we THROW AWAY ~75% of available levels (229 exist below IE alone) to dodge artifacts,
  - and we CAN'T just raise the cap -- the basis can't support those high-n levels.

CORRECTED CONCLUSION: the dominant gap to Bob is NOT the singlet CI -- it's BASIS
COMPLETENESS / DATA STARVATION. Our parameters are under-constrained because the
basis can't reach the high-n levels Bob fits, forcing an aggressive cap. The lever
that most closes 62->12 is a BIGGER BASIS (more configs: high-n Rydberg + doubly-
excited 3d configs), which would (a) admit the 150+ levels we cap away (better-
determined fit) AND (b) fix the 3s.10d/11d artifacts AND (c) supply the perturbers
for the singlet 1D mixing. Completeness is the master lever; the parameter recipe
(2 scales + curated G^1 + freed CI + zeta) is now clean and correct but data-starved.
