#!/bin/bash
# Build R. D. Cowan's atomic-structure suite (Kramida/NIST branch source)
# on modern Unix (macOS/Linux) with gfortran.
#
# Source: NIST PDR "A suite of atomic structure codes originally developed by
#   R. D. Cowan adapted for Windows-based personal computers" (Kramida, v.2021),
#   doi:10.18434/T4/1502500. Fortran sources from the for/ subdir.
#
# Flags: Kramida's recommended gfortran flags (gfortran_O3.bat) MINUS -malign-double
#   (x86-only; absent on Apple Silicon arm64), PLUS -std=legacy -fallow-argument-mismatch
#   (modern gfortran rejects F77 argument rank/type mismatches by default).
#
# Local source modifications vs pristine NIST source (see notes/build_notes.md):
#   - rcg11k.f SORT2: replaced local 'IMPLICIT REAL*8 (A-H,O-Z)' with
#     INCLUDE 'RCGPAR.F' (the routine was missing the include that supplies
#     PARAMETER KLAM=500000; matches sibling routines). Physics-neutral.
set -e
cd "$(dirname "$0")/src"
FLAGS=(-fshort-enums -ftracer -fno-backslash -O3 -std=legacy -fallow-argument-mismatch -I.)
mkdir -p ../bin
[ -f RCGPAR.F ] || ln -sf rcgpar.f RCGPAR.F   # case-insensitive INCLUDE portability

echo "Building Cowan suite with: gfortran ${FLAGS[*]}"
gfortran "${FLAGS[@]}" -o ../bin/rcn  RCN36K.F  && echo "  rcn  OK"
gfortran "${FLAGS[@]}" -o ../bin/rcn2 RCN2K.F   && echo "  rcn2 OK"
gfortran "${FLAGS[@]}" -o ../bin/rcg  rcg11k.f  && echo "  rcg  OK"
gfortran "${FLAGS[@]}" -o ../bin/rce  RCE20K.F  && echo "  rce  OK"
echo "Done. Binaries in build/bin/"
