#!/usr/bin/env python3
"""
Fetch observed atomic energy levels from the NIST Atomic Spectra Database (ASD)
and save a clean tab-separated table.  Reusable for any species.

These observed levels are what RCE fits the theoretical Hamiltonian parameters to
(the "semi-empirical" anchor).  We fetch once and cache in the repo so runs are
reproducible and don't depend on live network access.

Usage:
    python3 fetch_nist_levels.py "Mg I"  data/nist/MgI_levels.tsv
    python3 fetch_nist_levels.py "Fe II" data/nist/FeII_levels.tsv

Output columns (tab-separated): config  term  J  level_cm1  parity  is_predicted
  - parity: 'e' (even) or 'o' (odd), derived from the term '*' marker
  - is_predicted: 1 if NIST flagged the value (brackets/?) as not directly
    observed, else 0

Data source: NIST ASD, https://physics.nist.gov/asd  (please cite Kramida,
Ralchenko, Reader & NIST ASD Team).
"""
import sys
import urllib.parse
import urllib.request

ASD = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"


def fetch(spectrum):
    params = {
        "de": 0, "spectrum": spectrum, "units": 0,   # 0 = cm^-1
        "format": 3,                                  # 3 = tab-delimited
        "output": 0, "page_size": 15, "multiplet_ordered": 1,
        "conf_out": "on", "term_out": "on", "level_out": "on", "j_out": "on",
        "temp": "", "submit": "Retrieve Data",
    }
    url = ASD + "?" + urllib.parse.urlencode(params)
    # NIST rejects the default urllib UA with HTTP 403; send a browser-like UA.
    req = urllib.request.Request(url, headers={
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                      "AppleWebKit/537.36 (KHTML, like Gecko) "
                      "Chrome/120.0 Safari/537.36"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def unq(s):
    return s.strip().strip('"').strip()


def parse(raw):
    rows = []
    lines = raw.splitlines()
    if not lines:
        return rows
    # first line is the header
    for ln in lines[1:]:
        if not ln.strip():
            continue
        f = ln.split("\t")
        if len(f) < 5:
            continue
        conf, term, J, prefix, level = (unq(f[0]), unq(f[1]), unq(f[2]),
                                        unq(f[3]), unq(f[4]))
        suffix = unq(f[5]) if len(f) > 5 else ""
        if not conf or not level:
            continue
        # parity from the term '*' (odd) marker
        parity = "o" if "*" in term else "e"
        term_clean = term.replace("*", "")
        # predicted/uncertain if brackets or '?' present in level or markers
        flagged = any(c in (prefix + level + suffix) for c in "[]()?")
        try:
            lev = float(level.replace("[", "").replace("]", "")
                        .replace("(", "").replace(")", ""))
        except ValueError:
            continue
        rows.append((conf, term_clean, J, f"{lev:.3f}", parity,
                     "1" if flagged else "0"))
    return rows


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    spectrum, outpath = sys.argv[1], sys.argv[2]
    raw = fetch(spectrum)
    rows = parse(raw)
    if not rows:
        sys.stderr.write("No levels parsed; raw response head:\n" + raw[:500])
        sys.exit(2)
    with open(outpath, "w") as o:
        o.write("# NIST ASD energy levels for %s  (cm^-1)\n" % spectrum)
        o.write("# config\tterm\tJ\tlevel_cm1\tparity\tis_predicted\n")
        for r in rows:
            o.write("\t".join(r) + "\n")
    print(f"{spectrum}: wrote {len(rows)} levels -> {outpath}")


if __name__ == "__main__":
    main()
