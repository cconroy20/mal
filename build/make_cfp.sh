#!/bin/bash
# Generate the binary CFP (coefficients of fractional parentage) decks
# FOR072/FOR073/FOR074 that a freshly-compiled RCG requires.
# Unix replication of Kramida's MAKE_CFP.BAT. Run once after building RCG.
#
# Output: build/cfp/{FOR072,FOR073,FOR074}
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/../cowan_nist/extracted/CODE"
CFP="$HERE/cfp"
mkdir -p "$CFP"; cd "$CFP"

cp -f "$SRC/ING11.CFP" "$SRC/SENIOR" .
# MAKE_CFP.BAT: drop 8 header lines, strip trailing ';'-comments -> file 'ing11'
tail -n +9 ING11.CFP | sed 's/;.*$//' > ing11
echo "$CFP/" > cowan.cfg          # how the Cowan PC port locates its work dir

"$HERE/bin/rcg" < /dev/null > rcg_cfp.log 2>&1
if [ -s FOR072 ] && [ -s FOR073 ] && [ -s FOR074 ]; then
  echo "CFP decks created in $CFP:"
  ls -la FOR072 FOR073 FOR074 | awk '{print "  "$NF": "$5" bytes"}'
else
  echo "ERROR: CFP decks not created. See $CFP/rcg_cfp.log"; tail "$CFP/rcg_cfp.log"; exit 1
fi
