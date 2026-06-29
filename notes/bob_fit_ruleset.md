# RULESET: Which Radial Parameters to FREE vs FREEZE in a Cowan Semi-Empirical Fit

A physics-driven, implementable policy distilled from Bob Kurucz's production fits
of 10 ions (Na I, Mg I, Al I, Si I, Ca I, Ca II, Ti II, Cr I, Fe I, Fe II) plus
our own Mg I energy+gf experiments. Derived by an adversarial workflow: 24
candidate rules proposed across 4 lenses, each attacked against all 10 species;
5 survived; synthesized below. (Source data: notes/bob_fit_survey.json,
notes/bob_fit_detail.json, kurucz_ref/<xxyy>/c*.log.)

Parameter families: **EAV** (config centroid), **F^k** (direct Slater), **G^k**
(exchange Slater), **ZETA** (spin-orbit), **ALPHA/BETA** (Trees effective-operator
CI corrections for open d^n/f^n), **R^k / CI** (off-diagonal configuration-
interaction Slater integrals).

**Governing principle:** free only the radial integrals that the observed level
structure directly and locally constrains; freeze everything else at SCALED
Hartree-Fock; and capture configuration interaction through an explicit, complete
basis (plus diagonal Trees operators for open shells) — never by floating
off-diagonal CI integrals. What "directly constrains" means is dictated by the
open shell, so the free set is shell-aware, not a fixed kind-list.

---

## Rule 1 — Free one EAV per observed configuration (the universal backbone)

- **Rule:** Free the EAV (centroid) of every configuration with >=1 observed
  level; freeze EAVs of configs with no observed level.
- **Physics:** The centroid is the best-determined quantity from energies — data
  fixes each config's absolute placement directly and orthogonally to term
  structure. (Our Jacobian SVD independently shows EAVs are best-conditioned.)
- **Scope:** All ions.
- **Confidence:** 10/10, no exceptions. (EAV-share of the free set shrinks in
  open-d ions only because term-structure integrals grow, not because EAV drops.)

## Rule 2 — Freeze the structural background at SCALED HF (not raw HF)

- **Rule:** Every F^k/G^k/ZETA not explicitly freed is held FIXEDHF at a *scaled*
  HF value (Bob's scales ~0.6-0.85; verified ~0.79-0.82 on 3d-3d F2/F4),
  per-integral not global. Default = freeze; freeing must be earned by Rules 3-6.
- **Physics:** HF gives the right radial-integral *shape* but overestimates
  *magnitude* (correlation screens it). A per-integral screening factor corrects
  the bulk error cheaply and keeps a large basis tractable. A single GLOBAL F+G
  scale provably fails (Mg I test): at HF the singlet-triplet splitting is already
  too large while the direct term wants the opposite correction.
- **Scope:** All ions.
- **Confidence:** 10/10. Raw-HF freeze is a known failure: in the Mg I 122-config
  basis, freezing all CI/Slater at raw HF left 3s2 ~1700 cm^-1 too low and
  degraded gf. The fix is the scale factor, not freeing.

## Rule 3 — Free EXCHANGE G^k (esp. G^1) where a singlet-triplet/multiplet splitting is observed

- **Rule:** Free G^k for every low config whose exchange splitting is observed;
  free independently per observed Rydberg member (distinct group codes; fitted
  values decay smoothly with n). Freeze unobserved/high-Rydberg members at scaled
  HF. G^1 (3s-3p, 3s-np exchange) is the primary lever.
- **Physics:** G^k sets the multiplet (singlet-triplet) separation that observed
  levels pin directly and HF over-localizes — and it moves OPPOSITE the direct
  term, so no F+G global scale can fix it. Freeing G reshapes eigenvectors, so
  free only where data demands.
- **Scope:** All ions. PRIMARY non-EAV lever for closed-shell + 1-2 valence e-.
- **Confidence:** 10/10 mechanism. **Principled exception:** closed core + a
  SINGLE s/d valence electron with no resolved exchange multiplet frees ZERO G^k
  — Ca II (all 126 G fixed at scale 0.77-0.84); there ZETA is the primary lever.

## Rule 4 — Free SPIN-ORBIT ZETA where fine structure is observed AND splitting exceeds the precision floor

- **Rule:** Free ZETA on a level only if (a) its J-structure is resolved AND
  (b) the HF spin-orbit magnitude is above the list's precision floor (and above
  the residual scatter of the levels it splits). Always free the valence/open-
  shell ZETA (3d ZETA1 is typically the single most-freed parameter in iron-group
  ions). For a Rydberg series, free member-by-member only while HF ZETA stays
  above the floor; freeze the small-splitting tail (scale ~1.0).
- **Physics:** ZETA sets J-splitting within a term — near-diagonal, cleanly
  observed; HF spin-orbit is good to ~5-20%, so low-risk to free where resolved.
- **Scope:** All ions.
- **Confidence:** 10/10 broad claim (ZETA freed in every species). **Key
  qualification:** do NOT trigger on "observed" alone — gate on magnitude vs the
  ion's OWN floor. Ca II (ionized) freed down to ~1.2 cm^-1; Mg I to ~0.1; but
  Al I FROZE its fully-observed 3s2np 2P series for n>=5 (~<=4 cm^-1), freeing
  only 3p,4p. A fixed cm^-1 threshold or neutral-vs-ion switch is insufficient.

## Rule 5 — Free DIRECT F^k only for the open/same shell that builds the observed term structure

- **Rule:** Free F^k only where the direct integral is large and term structure
  constrains it — the open/same-l shell: F2/F4(dd) in open-d, F2(pp) in p2,
  F2(nl,nl) for same-shell pairs. Do NOT free F^k for a single electron over a
  closed-ish core (3s.nl), where the direct integral is small.
- **Physics:** Intra-shell direct integrals set the spread of terms within the
  open shell, tightly constrained by the observed multiplet pattern. Cross-shell/
  Rydberg direct integrals are small and weakly constrained -> freeze.
- **Scope:** Open-shell ions (p2, d^n, f^n); minimal in closed-shell ions.
- **Confidence:** 10/10. Closed-shell free almost no F (Na 1, Al 3, Mg 2, Ca I 8);
  open shells free many (Si I 18; Ti II 14, Cr I 7, Fe I 15, Fe II 33, F2/F4(dd)
  dominant). Note: by raw count G^k still >= F(dd) in Ti II/Cr I/Fe I — F(dd) is
  co-equal, not dominant over G.

## Rule 6 — For OPEN p2/d^n/f^n shells, add Trees ALPHA (always) and BETA (selectively)

- **Rule:** Only with an open p2/d^n/f^n shell, and only on its lowest richly-
  observed config(s):
  - **ALPHA** (L(L+1) Trees correction): free for every open p2/d^n/f^n config
    with observed structure — the universal open-shell effective-CI lever (freed
    in Si I, Ti II, Cr I, Fe I, Fe II). Its HF value is 0, so it repairs what HF
    structurally cannot.
  - **BETA** (seniority Trees correction; d^n/f^n only, none for p2): free
    SPARINGLY. In Bob's set, freed only for the IONS (Ti II, Fe II); the open-d
    NEUTRALS Cr I and Fe I keep BETA fully FIXED even on the ground config.
  - Rydberg/poorly-observed configs: keep ALPHA/BETA FIXED at a transferred
    (nonzero, hand-guess) constant — not 0.
- **Physics:** Trees alpha/beta are effective two-body operators absorbing
  far-configuration CI within the d^n/f^n shell as diagonal corrections — the
  tractable substitute for floating thousands of off-diagonal R^k.
- **Scope:** Open p2/d^n/f^n ions only.
- **Confidence:** ALPHA 5/5 open-shell species. BETA ions-yes / neutrals-no in
  this sample (mechanism for the ion/neutral split not yet established — flagged).

## Corollary — Do NOT free off-diagonal CI (R^k)

Across all 10 species Bob freed ~no R^k. CI is captured by an explicit COMPLETE
basis + (for open shells) the diagonal Trees operators, never by floating the
off-diagonal interaction integrals. Our own small-basis fits that freed CI were a
crutch for a too-small basis; with a complete basis + Rules 1-6, CI stays frozen
at scaled HF. (Caveat: requires the basis actually be complete near the levels of
interest — the lesson behind the Mg I full-basis ground-state depression.)

---

## DECISION PROCEDURE (for a new ion)

Given the ground configuration and the list of observed levels:

1. **Build a complete-enough basis** (single + key double excitations spanning the
   observed energy range), à la Bob's per-ion config list.
2. **Free one EAV per config that has an observed level** (Rule 1).
3. **Scale all HF F^k/G^k by ~0.85, ZETA by ~1.0; freeze everything by default**
   (Rule 2).
4. **Free G^k** for low configs / observed Rydberg members with a resolved
   singlet-triplet or multiplet splitting (Rule 3). If the ion is closed-core +
   single valence e- with no exchange multiplet (Ca II-like), skip G entirely.
5. **Free ZETA** where fine structure is resolved AND the HF spin-orbit magnitude
   exceeds the ion's precision floor; always free the valence/open-shell ZETA
   (Rule 4).
6. **If an open p2/d^n/f^n shell:** free F2/F4 of that shell (Rule 5) and ALPHA
   on its observed configs; free BETA only for ions, cautiously (Rule 6).
7. **Never free off-diagonal R^k** (Corollary). Apply the ridge prior toward
   scaled HF on everything left frozen so the fit stays well-posed.

## Open questions flagged by the adversarial pass

- The ZETA precision floor is ion-dependent (effective Z / ionization stage); a
  fixed cm^-1 threshold or neutral/ion switch is insufficient (Al I counterexample).
- BETA freed for ions but not neutrals in this sample; mechanism unknown.
- "Scaled HF" scale factors (~0.6-0.85) are per-integral; extracting Bob's actual
  per-integral scales (from b*.log / hf*.dat) remains to be done for our engine.
