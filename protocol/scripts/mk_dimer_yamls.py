#!/usr/bin/env python
"""Build two-copy (chain A + chain B) Boltz inputs for every p15 binder, receptor omitted.

Folding two copies of the binder asks a different question from the complex: does the
designed chain hold its fold on its own, or is it only ordered because the
receptor scaffolds it? Names are kept identical to the complex predictions so
the result can be compared residue-for-residue with its bound conformation.

Identical sequences are folded once and the mapping recorded, since several
designs can share a binder sequence.
"""
import glob, json, os, sys
from collections import defaultdict

ROOT, OUT, MAN = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)

by_seq = defaultdict(list)
for m in sorted(glob.glob(os.path.join(ROOT, "sh*", "manifest.json"))):
    sh = os.path.basename(os.path.dirname(m)).replace("sh", "")
    for y in json.load(open(m))["yamls"]:
        by_seq[y["binder_seq"]].append({"shard": sh, "name": y["name"],
                                        "binder_len": y["binder_len"]})
print(f"{sum(len(v) for v in by_seq.values())} complex designs -> "
      f"{len(by_seq)} unique binder sequences")

man = []
for i, (seq, members) in enumerate(sorted(by_seq.items())):
    nm = f"dimer_{i:04d}"
    with open(os.path.join(OUT, nm + ".yaml"), "w") as fh:
        fh.write("version: 1\nsequences:\n")
        for cid in ("A", "B"):
            fh.write(f"  - protein:\n      id: {cid}\n")
            fh.write(f"      sequence: {seq}\n      msa: empty\n")
    man.append({"mono": nm, "seq": seq, "length": len(seq),
                "n_complex_designs": len(members), "members": members})
os.makedirs(os.path.dirname(MAN), exist_ok=True)
json.dump({"n": len(man), "monomers": man}, open(MAN, "w"), indent=1)
L = {m["length"] for m in man}
print(f"wrote {len(man)} dimer yamls to {OUT}  (lengths {sorted(L)})")
print(f"manifest {MAN}")
