#!/usr/bin/env python
"""Agonist consensus for binder positions 1-10, and compliance of the p14 anchors.

Frame from the supplied table: consensus position 1 = GLP-1 His7 = binder
position 1, so consensus positions 1-10 map directly onto our 10-residue anchor.

Non-canonical residues (X = Aib / Acb / a-Me-Leu) cannot be produced by MPNN, so
where X leads at a position the canonical alternatives from the counts are used.
"""
import csv, glob, json, os, re, sys
from collections import Counter

# position: (allowed canonical residues in count order, top, conservation label)
CONSENSUS = {
    1:  (["H", "Y", "F"],           "H", "variable 82.6% H"),
    2:  (["G", "A", "S"],           "G", "variable; X=Aib leads 34.8% but is non-canonical"),
    3:  (["E", "D"],                "E", "acid class 65.2%; Q excluded to keep the acid requirement"),
    4:  (["G"],                     "G", "INVARIANT 100%"),
    5:  (["T"],                     "T", "INVARIANT 100%"),
    6:  (["F"],                     "F", "INVARIANT 100%"),
    7:  (["T", "I"],                "T", "near-invariant 95.7% T"),
    8:  (["S"],                     "S", "INVARIANT 100%"),
    9:  (["D"],                     "D", "INVARIANT 100%"),
    10: (["Y", "V", "L", "K"],      "Y", "variable; hyd class 56.5%"),
}
INVARIANT = [p for p, (a, _, _) in CONSENSUS.items() if len(a) == 1]
NATIVE = "HAEGTFTSDV"

if __name__ == "__main__":
    print("CONSENSUS RESTRAINTS for binder positions 1-10")
    print(f"  {'pos':>3s} {'GLP-1':>5s} {'allowed':14s} {'top':>3s}  note")
    for p in sorted(CONSENSUS):
        allowed, top, note = CONSENSUS[p]
        print(f"  {p:3d} {p+6:5d} {','.join(allowed):14s} {top:>3s}  {note}")
    print(f"\n  invariant positions: {INVARIANT}")
    print(f"  native GLP-1 anchor: {NATIVE}")
    ok = all(NATIVE[p-1] in CONSENSUS[p][0] for p in CONSENSUS)
    print(f"  native anchor satisfies every restraint: {ok}")

    # compliance of what the free redesign produced
    A = sys.argv[1] if len(sys.argv) > 1 else None
    if not A: sys.exit(0)
    anchors = Counter()
    for f in sorted(glob.glob(os.path.join(A, "p14_1092297_sh*_interactions.csv"))):
        for r in csv.DictReader(open(f)):
            a = r.get("anchor_seq", "")
            if len(a) == 10: anchors[a] += 1
    n = sum(anchors.values())
    print(f"\nCOMPLIANCE of the {n} p14 anchors ({len(anchors)} distinct) "
          f"with this consensus")
    viol = Counter(); full_ok = 0; inv_ok = 0
    for a, c in anchors.items():
        bad = [p for p in CONSENSUS if a[p-1] not in CONSENSUS[p][0]]
        if not bad: full_ok += c
        if not [p for p in bad if p in INVARIANT]: inv_ok += c
        for p in bad: viol[p] += c
    print(f"  satisfy ALL restraints        : {full_ok}/{n} ({100*full_ok/n:.1f}%)")
    print(f"  satisfy the invariant positions: {inv_ok}/{n} ({100*inv_ok/n:.1f}%)")
    print(f"  {'pos':>3s} {'expected':10s} {'violating':>10s}  most common wrong residue")
    for p in sorted(CONSENSUS):
        allowed = CONSENSUS[p][0]
        wrong = Counter(a[p-1] for a, c in anchors.items() for _ in range(c)
                        if a[p-1] not in allowed)
        v = viol.get(p, 0)
        mark = " <-- INVARIANT" if p in INVARIANT else ""
        print(f"  {p:3d} {','.join(allowed):10s} {v:7d} ({100*v/n:3.0f}%)  "
              f"{(wrong.most_common(1)[0][0] + ' x' + str(wrong.most_common(1)[0][1])) if wrong else '-':12s}{mark}")
    print("\n  most common p14 anchors vs native:")
    for a, c in anchors.most_common(6):
        diff = "".join("." if a[i] == NATIVE[i] else a[i] for i in range(10))
        print(f"    {a}  x{c:<4d}  vs native: {diff}")
