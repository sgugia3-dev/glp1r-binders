#!/usr/bin/env python
"""Copy the two-chain prediction for each viewable design, matched by SEQUENCE.

The single-chain and two-chain runs numbered their inputs independently, so
mono_0042 and dimer_0042 are not necessarily the same molecule. Matching on the
binder sequence is exact and does not depend on either run's ordering.
"""
import glob, json, os, shutil, sys

DATA, DIMROOT, DIMMAN, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
d = json.load(open(DATA))
dim = json.load(open(DIMMAN))["monomers"]
by_seq = {}
for e in dim:
    by_seq.setdefault(e["seq"], e["mono"])
print(f"two-chain entries indexed by sequence: {len(by_seq)}")

os.makedirs(OUT, exist_ok=True)
n = 0
for x in d["designs"]:
    if not x.get("file_complex"):      # only designs already viewable
        continue
    did = by_seq.get(x["seq"])
    if not did:
        print(f"  no two-chain match for {x['tag'][:40]}"); continue
    src = glob.glob(os.path.join(DIMROOT, did, "*_model_0.cif"))
    if not src:
        print(f"  no structure for {did}"); continue
    dst = os.path.join(OUT, f"{x['tag']}_pair.cif")
    shutil.copy(src[0], dst)
    x["file_pair"] = f"pdb/{x['tag']}_pair.pdb"
    x["pair_id"] = did
    n += 1
json.dump(d, open(DATA, "w"), indent=1)
print(f"copied {n} two-chain structures")
