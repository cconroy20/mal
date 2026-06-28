#!/bin/bash
# Worked example: drive the Cowan chain RCN -> RCN2 -> RCG on Sn7+ (4d^7, 4d^6 5p),
# the vendor's shipped example (IN36 + IN2). Produces energy levels, log gf, and
# (via our RCN dump patch) the converged radial wavefunctions in rwfn.dat.
#
# Prereqs: build/build.sh and build/make_cfp.sh have been run.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/.."
BIN="$ROOT/build/bin"
CFP="$ROOT/build/cfp"
EX="$ROOT/cowan_nist/extracted"
RUN="$HERE/sn7plus"

mkdir -p "$RUN"; cd "$RUN"
cp -f "$EX/work/IN36" in36          # RCN configuration input
cp -f "$EX/work/IN2"  IN2           # RCN2 control deck
cp -f "$CFP/FOR072" "$CFP/FOR073" "$CFP/FOR074" .   # CFP decks for RCG
cp -f "$EX/code/SENIOR" .
echo "$RUN/" > cowan.cfg            # RCG locates its working dir via this
rm -f rwfn.dat                      # RCN appends; start clean

echo "[1/3] RCN  (in36 -> out36, tape2n, rwfn.dat)"
"$BIN/rcn"  < /dev/null > rcn_run.log  2>&1
[ -f TAPE2N ] || ln -sf tape2n TAPE2N
echo "[2/3] RCN2 (IN2 + TAPE2N -> ING11)"
"$BIN/rcn2" < /dev/null > rcn2_run.log 2>&1
echo "[3/3] RCG  (ING11 + CFP decks -> OUTG11 with levels + log gf)"
"$BIN/rcg"  < /dev/null > rcg_run.log  2>&1

echo "Done. Outputs in $RUN:"
echo "  out36   : $(wc -l < out36) lines (RCN radial integrals)"
echo "  rwfn.dat: $(grep -c '^#' rwfn.dat) orbitals (converged P_nl(r))"
echo "  OUTG11  : $(awk '/ELEC DIP SPECTRUM/{s=1} /SUMS2,SUMGF/{s=0} s' OUTG11 | grep -cE '^ *[0-9]+ +-?[0-9]+\.[0-9]+ +[0-9]') E1 transitions (levels + log gf)"
