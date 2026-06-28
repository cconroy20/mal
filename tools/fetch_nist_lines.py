#!/usr/bin/env python3
"""
Fetch atomic transition data (wavelengths, log gf, gA, level IDs) from the NIST
Atomic Spectra Database (ASD) "lines" query, and save a clean tab-separated
table. The reference gf we validate our computed gf against.

The NIST lines CGI is finicky; the working parameter set (found empirically)
needs limits_type=0, I_scale_type=1, bibrefs=1, A_out=1, loggf_out=1, format=3,
line_out=1 (transition-probability lines only), and a wavelength window.

Output columns (tab-separated):
  ritz_wl_A  log_gf  gA  acc  conf_i  term_i  J_i  conf_k  term_k  J_k
(wavelength converted nm->Angstrom; only rows with a log_gf are kept.)

Usage:
    python3 fetch_nist_lines.py "Mg I" 2000 9000 data/nist/MgI_lines.tsv
        (wavelengths in Angstrom)

Data source: NIST ASD, https://physics.nist.gov/asd (cite Kramida, Ralchenko,
Reader & NIST ASD Team).
"""
import sys
import urllib.parse
import urllib.request

CGI = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"


def fetch(spectrum, lo_A, hi_A):
    params = {
        "spectra": spectrum, "limits_type": 0,
        "low_w": lo_A / 10.0, "upp_w": hi_A / 10.0, "unit": 1,  # nm
        "submit": "Retrieve Data", "de": 0, "I_scale_type": 1,
        "format": 3, "line_out": 1, "remove_js": "on", "en_unit": 1,
        "output": 0, "bibrefs": 1, "page_size": 15, "show_obs_wl": 1,
        "show_calc_wl": 1, "unc_out": 1, "order_out": 0, "show_av": 2,
        "tsb_value": 0, "A_out": 1, "allowed_out": 1, "forbid_out": 1,
        "conf_out": "on", "term_out": "on", "enrg_out": "on",
        "J_out": "on", "g_out": "on", "loggf_out": 1,
    }
    url = CGI + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def unq(s):
    return s.strip().strip('"').strip()


def parse(raw):
    rows = []
    lines = raw.splitlines()
    if not lines:
        return rows
    hdr = [unq(h) for h in lines[0].split("\t")]
    idx = {name: i for i, name in enumerate(hdr)}
    # wavelength column: prefer Ritz, vac then air; fall back to observed.
    wlcol = None
    is_air = False
    for cand in ("ritz_wl_vac(nm)", "obs_wl_vac(nm)",
                 "ritz_wl_air(nm)", "obs_wl_air(nm)"):
        if cand in idx:
            wlcol = cand
            is_air = "air" in cand
            break
    for ln in lines[1:]:
        f = [unq(x) for x in ln.split("\t")]
        if len(f) < len(hdr) - 1:
            continue

        def get(name):
            return f[idx[name]] if name in idx and idx[name] < len(f) else ""
        loggf = get("log_gf")
        wl = get(wlcol) if wlcol else ""
        if not loggf or not wl:
            continue
        try:
            lam_A = float(wl) * 10.0
            lg = float(loggf)
        except ValueError:
            continue
        rows.append((f"{lam_A:.3f}", f"{lg:.3f}", get("gA(s^-1)"), get("Acc"),
                     get("conf_i"), get("term_i"), get("J_i"),
                     get("conf_k"), get("term_k"), get("J_k")))
    return rows


def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    spectrum, lo, hi, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    rows = parse(fetch(spectrum, lo, hi))
    if not rows:
        sys.stderr.write("No lines with log gf parsed.\n")
        sys.exit(2)
    with open(out, "w") as o:
        o.write(f"# NIST ASD lines for {spectrum}  ({lo:.0f}-{hi:.0f} A, vac)\n")
        o.write("# ritz_wl_A\tlog_gf\tgA\tacc\tconf_i\tterm_i\tJ_i\t"
                "conf_k\tterm_k\tJ_k\n")
        for r in rows:
            o.write("\t".join(r) + "\n")
    print(f"{spectrum}: wrote {len(rows)} lines with log gf -> {out}")


if __name__ == "__main__":
    main()
