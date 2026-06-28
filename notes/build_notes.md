# Cowan code — macOS/gfortran build notes

**Date:** 2026-06-28
**Machine:** Apple Silicon (arm64) macOS, gfortran 15.2.0 (Homebrew GCC)
**Source:** Kramida/NIST branch, `Cowan_PC_2021.zip`, doi:10.18434/T4/1502500
(FOR/ subdir). Pristine source preserved in `cowan_nist/extracted/`; editable
build copy in `build/src/`.

## Result: all four programs build (rcn, rcn2, rcg, rce). ✅

## Compiler flags
Kramida's `gfortran_O3.bat` flags MINUS `-malign-double` (x86-only; not valid on
arm64), PLUS `-std=legacy -fallow-argument-mismatch`:

    gfortran -fshort-enums -ftracer -fno-backslash -O3 \
             -std=legacy -fallow-argument-mismatch -I.

## Source modifications vs pristine NIST source
1. **rcg11k.f, SUBROUTINE SORT2 (~line 1525):** replaced local
   `IMPLICIT REAL*8 (A-H,O-Z)` with `INCLUDE 'RCGPAR.F'`. SORT2 was the only
   routine dimensioning arrays with `KLAM` that lacked the include supplying
   `PARAMETER (KLAM=500000)` (sibling routines at ~1593, ~1614 have it).
   **Physics-neutral** — supplies the same constant + IMPLICIT typing.
2. **Symlink** `RCGPAR.F -> rcgpar.f` so the uppercase `INCLUDE 'RCGPAR.F'`
   resolves on case-sensitive filesystems (Linux); on macOS HFS it's moot.

## Warnings (benign, not errors)
Many `-fallow-argument-mismatch` warnings in rcg11k.f and RCE20K.F: legacy F77
patterns — REAL(8)/INTEGER(4) aliasing in sort routines, the `SORT2(...DUM,DUM)`
placeholder-arg convention (too-few-elements), one scalar/rank-1 `ERAS`,
INTEGER(4)→INTEGER(2) `NMAX`. These are intentional in the original code and
permitted by the flag. Flagged here in case any numerical regression appears.

## Post-build requirement (per README)
A freshly compiled RCG needs the binary CFP "decks" generated once via the
`make_cfp` step (the Windows `MAKE_CFP.BAT` runs RCG on a special input to write
FOR072/073/074). Must replicate on Unix before RCG produces spectra.

## TODO
- Generate CFP decks on Unix; confirm RCG runs end-to-end.
- Run a known example (WORK/IN2 etc.) and validate against expected output.
- Then: reproduce a known Kurucz Fe II gf (first deliverable, step 2).
