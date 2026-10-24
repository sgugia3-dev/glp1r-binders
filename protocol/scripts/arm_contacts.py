#!/usr/bin/env python
"""Measure final patch engagement on each p9 arm's output backbone.

The three arms share one fixed target in one frame (the target is not diffused),
so coordinates are directly comparable with no superposition. Reports, per arm,
binder contacts to the 86-89 patch and to the whole target, so the rise in the
patch_ncontacts trace can be attributed to the potential rather than to the
binder happening to dock near the patch anyway.
"""
import glob, math, os, sys

PATCH = [86, 87, 88, 89]
BINDER_LEN = 140


def load(p):
    b, t = [], {}
    for l in open(p):
        if not l.startswith("ATOM"): continue
        rn, an = int(l[22:26]), l[12:16].strip()
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if rn <= BINDER_LEN:
            b.append((rn, an, xyz))
        else:
            t.setdefault(rn - BINDER_LEN, []).append((an, xyz))
    return b, t


def switch(d, r_0=8.0, d_0=6.0):
    x = (d - d_0) / r_0
    if abs(x - 1.0) < 1e-9: return 0.5
    return (1 - x ** 6) / (1 - x ** 12)


print(f"{'arm':22s} {'patch<5A':>9s} {'patch<8A':>9s} {'minDist':>8s} "
      f"{'switchSum':>10s} {'tgt<5A':>7s}")
rows = []
for d in sorted(glob.glob(sys.argv[1] + "/arm*")):
    pdbs = sorted(glob.glob(d + "/des_*.pdb"))
    if not pdbs: continue
    lbl = os.path.basename(d)
    for p in pdbs:
        b, t = load(p)
        bca = [x for r, a, x in b if a == "CA"]
        n5 = n8 = 0
        mind = 1e9
        ssum = 0.0
        for pr in PATCH:
            for an, tx in t.get(pr, []):
                for bx in bca:
                    dd = math.dist(bx, tx)
                    mind = min(mind, dd)
            # CA-CA for the switching sum, matching the potential's definition
            tca = [x for an, x in t.get(pr, []) if an == "CA"]
            for tx in tca:
                for bx in bca:
                    dd = math.dist(bx, tx)
                    ssum += switch(dd)
                    if dd < 5: n5 += 1
                    if dd < 8: n8 += 1
        tot5 = 0
        for pr, ats in t.items():
            tca = [x for an, x in ats if an == "CA"]
            for tx in tca:
                for bx in bca:
                    if math.dist(bx, tx) < 5: tot5 += 1
        print(f"{lbl:22s} {n5:9d} {n8:9d} {mind:8.2f} {ssum:10.2f} {tot5:7d}")
        rows.append((lbl, n5, n8, mind, ssum, tot5))
if rows:
    base = next((r for r in rows if "arm0" in r[0]), None)
    print()
    if base:
        print(f"baseline arm0 patch switchSum = {base[4]:.2f}, min dist {base[3]:.2f} A")
        for r in rows:
            if r is base: continue
            print(f"  {r[0]}: switchSum {r[4]:.2f} ({r[4]-base[4]:+.2f} vs baseline), "
                  f"min dist {r[3]:.2f} A ({r[3]-base[3]:+.2f})")
    print("\nNOTE: one design per arm. This shows the wiring works and the direction,"
          "\nnot a significant effect. Attribution needs many designs per arm.")
