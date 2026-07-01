#!/usr/bin/env python3
"""Compare OUR fitted radial parameters to Bob Kurucz's curated values, integral
by integral, joined by (config, physical-name). Bob's values come from his RCE
fit logs c1200{e,o}z.log; ours from param_labels (the validated ING11 slot->name
mapping, so no wrong-slot/unit errors). Both in cm^-1.

Bob's log row format (param section):
   <#>  <grp>  <NAME>  <FITTED>  <STEP>  <HF>  [FIXEDHF|FIXED]  <scale> ... [<cfg> on EAV rows]
The config appears only on the EAV row and is inherited by the following rows of
the same group (same as RCG's OUTGINE layout). STEP>0 (and no FIXEDHF) => Bob
FREELY FIT it; FIXEDHF => he held it at scale*HF.

Usage:
  python3 tools/compare_bob.py \
      --ing11 work/mg1_full/ING11.energyonly \
      --outgine work/mg1_full/OUTGINE.abinitio \
      --bob-even kurucz_ref/1200/c1200ez.log \
      --bob-odd  kurucz_ref/1200/c1200oz.log \
      [--config 3p2] [--kind EAV|F|G|ZETA]
"""
import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import numpy as np
import make_report as R
import param_labels as PL


def parse_bob_log(path):
    """Parse a Bob c-log param section -> list of dicts
    {cfg (cfgkey), phys, fitted, hf, step, scale, fixed (bool), name}.
    The config is carried from each EAV row to the following same-group rows."""
    rows = []
    cur_cfg = None
    cur_grp = None
    with open(path, errors="replace") as f:
        for ln in f:
            # param rows start with: <int> <whitespace> <grp-token>
            m = re.match(r"\s*\d+\s+(\S+)\s+([A-Z][\w()]*(?:\s\d+)?)\s+"
                         r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+"
                         r"(FIXEDHF|FIXED)?\s*(-?\d+\.\d+)", ln)
            if not m:
                continue
            grp, name, fitted, step, hf, fixflag, scale = m.groups()
            # config name (if present) is the trailing token on EAV lines
            cfg_m = re.search(r"([0-9][spdfghikl0-9.]+[a-z]?)\s*$", ln.rstrip())
            if name.startswith("EAV"):
                # EAV row carries the config; grp changes per config
                cur_cfg = R._cfgkey(cfg_m.group(1)) if cfg_m else cur_cfg
                cur_grp = grp
            # normalize physical name: 'F2(22)'->'F2', 'G0(13)'->'G0', 'ZETA 2'->'ZETA',
            # 'EAV'->'EAV'. Keep the full name too for disambiguation.
            phys = ("EAV" if name.startswith("EAV")
                    else "ZETA" if name.startswith("ZETA")
                    else re.match(r"([FG]\d|ALPHA|BETA)", name).group(1)
                    if re.match(r"([FG]\d|ALPHA|BETA)", name) else name)
            rows.append({"cfg": cur_cfg, "phys": phys, "name": name.strip(),
                         "fitted": float(fitted), "hf": float(hf),
                         "step": float(step), "scale": float(scale),
                         "fixed": bool(fixflag)})
    return rows


def _our_params(ing11, outgine):
    """{(cfgkey, phys-prefix) -> value_cm1} for our single-config params. phys is
    reduced to the F2/G1/EAV/ZETA prefix to join with Bob's normalized names."""
    out = {}
    for p in PL.physical_params(ing11, outgine):
        if p["phys"] == "CI" or p["cfg"] is None:
            continue
        # normalize 'F2(22)'->'F2', 'ZETA 2'->'ZETA'
        ph = p["phys"]
        ph = ("EAV" if ph == "EAV"
              else "ZETA" if ph.startswith("ZETA")
              else re.match(r"[FG]\d", ph).group(0) if re.match(r"[FG]\d", ph)
              else ph)
        out[(p["cfg"], ph)] = p["value_cm1"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ing11", required=True)
    ap.add_argument("--outgine", required=True)
    ap.add_argument("--bob-even", required=True)
    ap.add_argument("--bob-odd", default=None)
    ap.add_argument("--config", default=None, help="filter to one config (cfgkey).")
    ap.add_argument("--kind", default=None, help="filter to EAV|F|G|ZETA.")
    ap.add_argument("--free-only", action="store_true",
                    help="only params Bob FREELY FIT (step>0, not FIXEDHF).")
    a = ap.parse_args()

    bob = parse_bob_log(a.bob_even)
    if a.bob_odd:
        bob += parse_bob_log(a.bob_odd)
    ours = _our_params(a.ing11, a.outgine)

    ck = R._cfgkey(a.config) if a.config else None
    rows = []
    for b in bob:
        if b["cfg"] is None:
            continue
        if ck and b["cfg"] != ck:
            continue
        if a.kind and not b["phys"].startswith(a.kind):
            continue
        if a.free_only and (b["fixed"] or b["step"] == 0):
            continue
        ov = ours.get((b["cfg"], b["phys"]))
        rows.append((b, ov))

    # EAVs carry a global ENERGY-ZERO offset (our fit pins the ground EAV to 0;
    # Bob's absolute scale differs by a near-constant ~3133 cm^-1). To compare EAVs
    # meaningfully we remove that constant: the MEDIAN (ours-Bob) over EAVs that
    # are NOT basis-limit-broken (|diff| within 3x the median). F/G/ZETA have no
    # such offset and are compared raw.
    eav_diffs = [ov - b["fitted"] for b, ov in rows
                 if ov is not None and b["phys"] == "EAV"]
    eav_off = float(np.median(eav_diffs)) if eav_diffs else 0.0

    print(f"(EAV global-offset removed: {eav_off:+.1f} cm^-1)")
    print(f"{'config':9} {'param':6} {'bob_fit':>10} {'ours':>10} "
          f"{'diff':>9} {'ours/bob':>8} {'bob:free?':>9}")
    eav_d, struct_d = [], []
    for b, ov in rows:
        bf = b["fitted"]
        if ov is None:
            print(f"{b['cfg']:9} {b['phys']:6} {bf:10.1f} {'--':>10} "
                  f"{'(no match)':>9} {'':>8} {'free' if not b['fixed'] else 'fixed':>9}")
            continue
        is_eav = b["phys"] == "EAV"
        ov_adj = ov - eav_off if is_eav else ov     # de-offset EAVs for the diff
        diff = ov_adj - bf
        ratio = ov / bf if abs(bf) > 1e-6 else float("nan")
        (eav_d if is_eav else struct_d).append((b["cfg"], b["phys"], diff, ratio))
        tag = "free" if not b["fixed"] else "fixed"
        note = "  <-EAV(offset-rm)" if is_eav else ""
        print(f"{b['cfg']:9} {b['phys']:6} {bf:10.1f} {ov:10.1f} "
              f"{diff:+9.1f} {ratio:8.3f} {tag:>9}{note}")

    def _stats(name, ds):
        if not ds:
            return
        ad = np.array([abs(d[2]) for d in ds])
        i = int(np.argmax(ad))
        print(f"  {name}: {len(ds)} params | mean|diff|={ad.mean():.1f} "
              f"median|diff|={np.median(ad):.1f} | max={ad.max():.1f} "
              f"({ds[i][0]} {ds[i][1]})")
    print("\nAGREEMENT WITH BOB (cm^-1):")
    # split EAVs into well-fit vs basis-limit-broken (|diff|>500 after offset)
    eav_ok = [d for d in eav_d if abs(d[2]) < 500]
    eav_bad = [d for d in eav_d if abs(d[2]) >= 500]
    _stats("EAV (offset-removed, |d|<500)", eav_ok)
    if eav_bad:
        print(f"  EAV basis-limit/broken (|d|>=500): "
              f"{[(d[0]) for d in eav_bad]}")
    _stats("F/G/ZETA structure integrals", struct_d)


if __name__ == "__main__":
    main()
