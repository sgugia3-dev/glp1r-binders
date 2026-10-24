#!/usr/bin/env python
"""Build the rfd3 input for the combined A27 / A42 / A86 hotspot set.

Clones prod_r3.json (the validated P2 protocol) key for key and swaps only
select_hotspots, so contig, infer_ori_strategy, is_non_loopy and the per-model
inputs stay exactly as P2 and P3 ran them.

Atoms are the ones those two campaigns actually steered on:
  A27:CD1 and A42:OG1 from P2, A86:CD2 from P3.
"""
import json, os, sys

src, dst = sys.argv[1], sys.argv[2]
HOT = {"A27": "CD1", "A42": "OG1", "A86": "CD2"}

d = json.load(open(src))
print(f"source {os.path.basename(src)}: {len(d)} keys")
out = {}
for k, v in d.items():
    nk = k.replace("prod_r3", "w3_r27_42_86")
    nv = dict(v)
    old = nv.get("select_hotspots")
    nv["select_hotspots"] = dict(HOT)
    out[nk] = nv
    print(f"  {k} -> {nk}")
    print(f"      contig   {nv.get('contig')}")
    print(f"      hotspots {old}  ->  {nv['select_hotspots']}")
json.dump(out, open(dst, "w"), indent=2)

# verify everything except the hotspots is untouched
chk = json.load(open(dst))
assert len(chk) == len(d), "key count changed"
for (ko, vo), (kn, vn) in zip(d.items(), chk.items()):
    for f in vo:
        if f == "select_hotspots": continue
        assert vo[f] == vn[f], (f, vo[f], vn[f])
    assert vn["select_hotspots"] == HOT
print(f"\nverified: only select_hotspots differs from the P2 protocol")
print(f"wrote {dst}")
