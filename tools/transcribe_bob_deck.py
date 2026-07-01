#!/usr/bin/env python3
"""Transcribe Bob Kurucz's FITTED radial parameters (from his RCE logs
c1200{e,o}z.log) into a copy of OUR ING11, for the RCE-reproduction test: run the
result through OUR RCG and see whether we reproduce his level residuals -- isolating
"his parameters" vs "our forward model" as the source of the 3s.nd 1D gap.

WHAT IT TRANSCRIBES: every FREE parameter of Bob's that maps to a slot in our deck
-- EAV (centroids), Gk/Fk (Slater integrals), ZETA (spin-orbit). CI is left at our
values (Bob's CIs == ours, 0.8*HF, already verified). Out-of-basis params (n>=12,
not in our 122-config deck) are skipped and reported.

ZERO-POINT: Bob pins the ground 3s2 EAV to his own zero; our EAVs sit on a
different absolute zero (median-offset convention). EAVs are only meaningful up to
a common additive constant, so we REMOVE a single zero-point = the mean (Bob - ours)
over the low, well-determined EAVs both decks share, applied to every transcribed
EAV. (Gk/Fk/ZETA are absolute -- no zero-point.) The zero-point choice cannot change
the physics: RCG level energies are invariant under a common EAV shift; only the
per-config DIFFERENCES (which we transcribe faithfully) matter.

Usage:
  tools/transcribe_bob_deck.py --base work/mg1_full/ING11.energyonly \
      --out work/mg1_full/ING11.bobdeck [--dry-run]
"""
import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import numpy as np
import make_report as R
import param_labels as PL
import ing11_params as IP
import parse_bob_params as B

EVEN = os.path.join(ROOT, "kurucz_ref", "1200", "c1200ez.log")
ODD = os.path.join(ROOT, "kurucz_ref", "1200", "c1200oz.log")


def bob_free_params():
    """All of Bob's FREE params (even+odd), keyed by (cfgkey, kind, k)."""
    out = {}
    for fn in (EVEN, ODD):
        for p in B.parse(fn):
            if p["fixed"] or not p["config"]:
                continue
            # skip CI (pair configs) -- we keep ours; take single-config params only
            cfg = p["config"].split("-")[0].strip()
            if " " in cfg or "-" in p["config"]:
                continue
            ck = R._cfgkey(cfg)
            out[(ck, p["kind"], p.get("k"))] = p
    return out


def our_slots(base_ing11):
    """Our writable params keyed by (cfgkey, kind, k) -> physical_param dict."""
    slots = {}
    for p in PL.physical_params(base_ing11, os.path.join(
            ROOT, "work", "mg1_full", "OUTGINE.abinitio"), strict=False):
        if p.get("kind") in ("EAV", "G", "F", "ZETA"):
            slots[(p.get("cfg"), p["kind"], p.get("k"))] = p
    return slots


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="our ING11 to copy & overwrite")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    bob = bob_free_params()
    slots = our_slots(a.base)

    # --- zero-point from shared low EAVs (Bob - ours), in cm^-1 ---
    # The Bob-ours EAV diff is a clean ~+3135 for almost every config (a pure
    # zero-point). A few configs are OUTLIERS because OUR value there is unreliable,
    # not because Bob differs: the high-n d basis-edge members (3s.9d/10d/11d, where
    # our identity labeling is fragile) and 3p.4s (a doubly-excited config our model
    # misplaces). We (a) estimate the zero-point robustly (median), and (b) EXCLUDE
    # those configs from EAV transcription entirely -- writing Bob's value onto a
    # slot our own model can't place is meaningless. 3p2's -1170-from-zeropoint diff
    # is REAL physics (the perturber sits lower) and is KEPT.
    EAV_EXCLUDE = {"3s.9d", "3s.10d", "3s.11d", "3p.4s"}
    shared_eav = []
    for (ck, kind, k), bp in bob.items():
        if kind != "EAV":
            continue
        sp = slots.get((ck, "EAV", None))
        if sp is not None:
            shared_eav.append((ck, bp["value"] - sp["value_cm1"]))
    diffs = np.array([d for ck, d in shared_eav if ck not in EAV_EXCLUDE])
    zp = float(np.median(diffs))
    print(f"shared EAVs: {len(shared_eav)}   zero-point (median Bob-ours, "
          f"outliers excluded) = {zp:.1f} cm^-1   spread(std) = {diffs.std():.1f}")

    # --- build the value overrides keyed by ING11 param key ---
    raw, params = IP.parse(a.base)
    by_key = {p["key"]: i for i, p in enumerate(params)}
    values = [p["value"] for p in params]

    applied = {"EAV": 0, "G": 0, "F": 0, "ZETA": 0}
    skipped = []
    for (ck, kind, k), bp in bob.items():
        if kind == "EAV" and ck in EAV_EXCLUDE:
            skipped.append((ck, kind, k, "excluded (our slot unreliable)"))
            continue
        sp = slots.get((ck, kind, k))
        if sp is None:
            skipped.append((ck, kind, k))
            continue
        key = sp["key"]
        if key not in by_key:
            skipped.append((ck, kind, k, "no ING11 key"))
            continue
        # Bob values are in cm^-1; our ING11 stores kK (value*1000 = cm^-1).
        v_cm = bp["value"] - (zp if kind == "EAV" else 0.0)
        values[by_key[key]] = v_cm / 1000.0
        applied[kind] += 1

    print(f"applied: {applied}   (total {sum(applied.values())})")
    # report in-basis skips (n<12) separately from out-of-basis
    import re
    real_skips = [s for s in skipped
                  if max([int(m) for m in re.findall(r'(\d+)[spdfghik]', s[0])]
                          or [0]) < 12]
    print(f"skipped: {len(skipped)}  (of which in-basis n<12: {len(real_skips)})")
    for s in real_skips[:15]:
        print("   in-basis skip:", s)

    if a.dry_run:
        print("(dry run -- not written)")
        return
    IP.write(raw, params, values, a.out)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
