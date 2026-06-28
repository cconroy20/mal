# mal — Modern Atomic Linelists

A project to take over, modernize, and extend the production of atomic line
lists in the lineage of Robert L. Kurucz's work — the energy levels and
oscillator strengths (gf-values) that underpin stellar-atmosphere opacity and
spectral synthesis.

## Why

Kurucz's line lists are critical infrastructure for stellar spectroscopy, but
the gf-*generating* engine (his private branch of R. D. Cowan's atomic-structure
codes) is not publicly available. This project rebuilds the generating capability
around the **public Cowan code** (the exact lineage Kurucz forked from ~1970) and
modernizes the one stage Kurucz never updated in ~50 years: the semi-empirical
least-squares fit (RCE).

## Strategy (two tiers)

- **Tier 1 — industrialize + UQ (the spine):** apply validated modern fitting
  methods (regularization; orthogonal operators where they scale) plus per-line
  uncertainty quantification across all ions → a complete, better,
  uncertainty-bearing replacement for `gfall`. The hard/novel part is automating
  the per-ion "black magic" of the RCE fit.
- **Tier 2 — Bayesian + spectral coupling (the research bet):** a fully Bayesian
  fit that folds the solar/stellar spectrum *into* the likelihood, coupling lines
  through shared radial parameters → spectrum-informed gf with posteriors. Forward
  model: the in-house ATLAS12/SYNTHE port (`~/kurucz/atlas12`). Confirmed novel.

See `docs/` and the planning notes for the full rationale and research basis.

## Status

- [x] **Cowan kernel builds & runs on macOS/arm64 (gfortran 15).** All four
      programs (RCN, RCN2, RCG, RCE) compile; RCG generates the CFP decks and
      runs to normal exit. See `notes/build_notes.md`.
- [ ] Run a known example through the full RCN→RCN2→RCG→RCE chain.
- [ ] Reproduce a known Kurucz Fe II gf (first deliverable).
- [ ] Prototype a modern regularized fit + per-line gf UQ on Fe II.

## Layout

    build/          editable Cowan source + build scripts
      src/          Fortran (from NIST/Kramida branch; see build_notes.md for edits)
      build.sh      compile all four programs with gfortran -> build/bin/
      make_cfp.sh   generate the binary CFP decks RCG needs (run once)
      cfp/          generated decks + run scratch (gitignored)
    cowan_nist/     vendored NIST/Kramida package (zip & Windows .exe/.dll removed)
      extracted/for/    pristine Fortran source (RCN36K.F, RCN2K.F, rcg11k.f, RCE20K.F)
      extracted/code/   perl utilities, .bat reference scripts, ING11.CFP, SENIOR
      extracted/work/   example input decks (IN2, IN36, ...)
      extracted/*_DOC.txt, readme.cowan.htm   the Cowan documentation
    tools/          pipeline/orchestration code (to come)
    work/           per-ion working calculations (to come)
    notes/, docs/   build notes, methodology, references

## Building

    bash build/build.sh       # -> build/bin/{rcn,rcn2,rcg,rce}
    bash build/make_cfp.sh    # -> build/cfp/{FOR072,FOR073,FOR074}

Requires `gfortran`. Source provenance and the (physics-neutral) modifications
vs. the pristine NIST source are documented in `notes/build_notes.md`.

## Provenance & license

Cowan code: A. Kramida (NIST) branch, "A suite of atomic structure codes
originally developed by R. D. Cowan adapted for Windows-based personal
computers," doi:10.18434/T4/1502500, NIST open license. Original author
R. D. Cowan (LANL); spec is Cowan, *The Theory of Atomic Structure and
Spectra* (1981).
