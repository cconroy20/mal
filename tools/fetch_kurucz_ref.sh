#!/bin/bash
# Fetch Bob Kurucz's per-ion production files from http://kurucz.harvard.edu/atoms
# into kurucz_ref/<xxyy>/ (gitignored; large). These are his ACTUAL Cowan-lineage
# semi-empirical fit decks and outputs -- the ground truth this project stewards.
#
# Directory naming: xxyy = element_number*100 + charge. Mg I = 1200, Fe II = 2601.
#
# Key files per ion (see http://kurucz.harvard.edu/atoms.html for the full legend):
#   b<xxyy>{e,o}.com  RCE least-squares fit INPUT decks (even/odd) -- Bob's
#                     configuration list, observed levels, CI structure, and which
#                     parameters were freed. The authoritative version of what
#                     tools/build_ine20.py reconstructs.
#   b<xxyy>{e,o}.log  RCE fit OUTPUTS: fitted radial (Slater/CI) parameters + levels.
#   c<xxyy>{e,o}.log  all eigenvalues + 3 strongest LS eigenvector components
#                     (Bob's level identifications -- compare to parse_compositions).
#   hf<xxyy>*.dat     Hartree-Fock starting integrals (ab-initio radial params).
#   gf<xxyy>.gam      energy levels (J, Lande g, A-sums, eigenvector components).
#   gf<xxyy>.lines/.pos/.all   the line lists (gfall content for this ion);
#                     .all = computed with LAB gf substituted where better.
#   gf<xxyy>*.lab     laboratory gf values Bob trusted (e.g. Fuhr&Wiese for Fe).
#
# Usage:  tools/fetch_kurucz_ref.sh 1200 2601 [...]
#         tools/fetch_kurucz_ref.sh 2601 --all     # grab the ion's whole dir
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BASE="http://kurucz.harvard.edu/atoms"
CORE=(b%se.com b%so.com b%se.log b%so.log c%se.log c%so.log \
      gf%s.gam gf%s.lines gf%s.pos gf%s.all)

want_all=0
ions=()
for a in "$@"; do [ "$a" = "--all" ] && want_all=1 || ions+=("$a"); done
[ ${#ions[@]} -eq 0 ] && { echo "usage: $0 <xxyy> [xxyy ...] [--all]"; exit 1; }

for ion in "${ions[@]}"; do
  dst="$ROOT/kurucz_ref/$ion"; mkdir -p "$dst"
  echo "[$ion] -> $dst"
  if [ "$want_all" = 1 ]; then
    # scrape the directory index and pull every linked data file
    curl -sk -L "$BASE/$ion/" \
      | grep -oiE 'href="[^"]+"' | sed -E 's/href="//I;s/"//' \
      | grep -viE '\.(html?|gif|css|ico)$|^[?/]' | sort -u \
      | while read -r f; do curl -sk -L "$BASE/$ion/$f" -o "$dst/$f"; done
  else
    for tpl in "${CORE[@]}"; do
      f=$(printf "$tpl" "$ion")
      curl -sk -L "$BASE/$ion/$f" -o "$dst/$f" 2>/dev/null || true
    done
  fi
  # drop 404 stubs (tiny HTML error bodies)
  find "$dst" -type f -size -50c -delete 2>/dev/null || true
  echo "  $(ls "$dst" | wc -l | tr -d ' ') files"
done
