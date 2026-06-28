#!/bin/bash
# Warm-up example: drive RCN -> RCN2 -> RCG on Mg I.
# Configurations: 3s^2 (even) and 3s3p (odd) over a neon-like 1s2 2s2 2p6 core.
# Produces the Mg I resonance line 3s2 1S0 - 3s3p 1P1 (and the 3P term that, once
# we add RCE + spin-orbit mixing, yields the 4571 A intercombination line).
#
# The committed deck is work/mg1/in36; this script supplies the rest of the
# environment and runs the chain. Prereqs: build/build.sh, build/make_cfp.sh.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/.."
BIN="$ROOT/build/bin"; CFP="$ROOT/build/cfp"; EX="$ROOT/cowan_nist/extracted"
RUN="$HERE/mg1"
cd "$RUN"

# in36 is committed; bring in the rest of the run environment
cp -f "$EX/work/IN2" IN2
cp -f "$CFP/FOR072" "$CFP/FOR073" "$CFP/FOR074" .
cp -f "$EX/code/SENIOR" .
echo "$RUN/" > cowan.cfg
rm -f rwfn.dat

echo "[1/3] RCN"
"$BIN/rcn"  < /dev/null > rcn_run.log  2>&1
[ -f TAPE2N ] || ln -sf tape2n TAPE2N
echo "[2/3] RCN2"
"$BIN/rcn2" < /dev/null > rcn2_run.log 2>&1
echo "[3/3] RCG"
"$BIN/rcg"  < /dev/null > rcg_run.log  2>&1

echo "Done. Computed E1 lines (3s2 - 3s3p):"
awk '/ELEC DIP SPECTRUM     \(ENERGIES/{s=1} /SUMS2,SUMGF/{s=0}
     s && /LAMBDA\(A\)/{hdr=1; next}
     s && hdr && /^ *[0-9]+ +-?[0-9]/{printf "  %s\n",$0}' "$RUN/OUTG11" | head
