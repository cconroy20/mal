#!/bin/bash
# Warm-up example: drive the full Cowan chain on Mg I, including the RCE
# semi-empirical fit to the observed (NIST) levels.
#   RCN -> RCN2 -> RCG  (ab initio levels, gf; writes OUTGINE template)
#   build INE20 from OUTGINE + NIST levels (substitute observed T-values,
#       free the physical parameters)
#   RCE  (least-squares fit of the radial parameters to the observed levels)
#
# Even configs: 3s2, 3s4s, 3s3d, 3s5s, 3s4d.  Odd: 3s3p, 3s4p, 3s5p.
# The committed deck is work/mg1/in36. Prereqs: build/build.sh, build/make_cfp.sh.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/.."
BIN="$ROOT/build/bin"; CFP="$ROOT/build/cfp"; EX="$ROOT/cowan_nist/extracted"
RUN="$HERE/mg1"
NIST="$ROOT/data/nist/MgI_levels.tsv"
cd "$RUN"

cp -f "$EX/work/IN2" IN2
cp -f "$CFP/FOR072" "$CFP/FOR073" "$CFP/FOR074" .
cp -f "$EX/code/SENIOR" .
echo "$RUN/" > cowan.cfg
rm -f rwfn.dat

echo "[1/4] RCN"
"$BIN/rcn"  < /dev/null > rcn_run.log  2>&1
[ -f TAPE2N ] || ln -sf tape2n TAPE2N
echo "[2/4] RCN2"
"$BIN/rcn2" < /dev/null > rcn2_run.log 2>&1
echo "[3/6] RCG  (writes OUTGINE template, OUTG11 with ab initio levels+gf)"
"$BIN/rcg"  < /dev/null > rcg_run.log  2>&1
cp -f OUTG11 OUTG11.abinitio          # keep the pre-fit spectrum for comparison
cp -f ING11  ING11.abinitio           # keep the ab initio parameter deck

echo "[4/6] build INE20 (substitute NIST observed levels, free parameters)"
python3 "$ROOT/tools/build_ine20.py" --outgine OUTGINE --nist "$NIST" \
    --outg11 OUTG11.abinitio --free-ci-pairs 1-6 --out INE20
cp -f INE20 OUTGINE                    # RCE reads the file named OUTGINE

echo "[5/6] RCE  (semi-empirical least-squares fit)"
"$BIN/rce"  < /dev/null > rce_run.log  2>&1

echo "[6/6] fitted gf: substitute RCE params into ING11, re-run RCG"
python3 "$ROOT/tools/subst_fitted_params.py" \
    --parvals PARVALS --ing11 ING11.abinitio --out ING11.fit
cp -f ING11.fit ING11
"$BIN/rcg"  < /dev/null > rcg_fit_run.log 2>&1
cp -f OUTG11 OUTG11.fitted            # the fitted-parameter spectrum (fitted gf)

echo "Done. RCE convergence:"
grep -iE "Iteration No.*AVDEV" rce_run.log | tail -6 | sed 's/^/  /'
