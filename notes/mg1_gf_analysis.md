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
energy-only fit threw away. NEXT: tune lambda; add per-line gf MC uncertainties
on top of this loop; try ridge-toward-ab-initio (--ridge) for ill-conditioned
ions; carry to Fe II.

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
