#!/usr/bin/env python
"""Pick the top distinct backbones from the band and write a PyMOL render script.

One row per sequence in the band CSV, but several sequences can share a backbone,
so selection is by backbone to avoid presenting the same fold repeatedly.

In the Boltz complex the binder is residues 1-140 and the receptor is renumbered
141-549, so receptor residue N appears at 140+N. Hotspot numbers are translated
accordingly, per source campaign, since the four lineages used three different sets.
"""
import csv, glob, os, shutil, sys

BAND, ROOT, OUTDIR, N = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
HOTS = {"p3": [86, 87, 88, 89], "p10": [27, 42, 86],
        "p12": [27, 42, 86], "p13": [27, 42, 86, 95]}
OFF = 140

rows = [r for r in csv.DictReader(open(BAND))
        if r["gated"] == "True" and r["full"] == "True" and r["score"]]
rows.sort(key=lambda r: -float(r["score"]))
seen, picked = set(), []
for r in rows:
    if r["bb_tag"] in seen: continue
    seen.add(r["bb_tag"]); picked.append(r)
    if len(picked) >= N: break
print(f"{len(rows)} gated+full rows -> {len(picked)} distinct backbones picked")

os.makedirs(OUTDIR, exist_ok=True)
lines = ["bg_color white", "set ray_opaque_background, 1",
         "set cartoon_transparency, 0", "set antialias, 2",
         "set ray_shadows, 0", "set specular, 0.2",
         "set cartoon_fancy_helices, 1", "set opaque_background, 1"]
manifest = []
for i, r in enumerate(picked, 1):
    d = os.path.join(ROOT, "sh" + r["shard"], "boltz_results_yaml",
                     "predictions", r["name"])
    cif = glob.glob(os.path.join(d, "*_model_0.cif"))
    if not cif:
        print(f"  missing structure for {r['name']}"); continue
    dst = os.path.join(OUTDIR, f"cand{i:02d}_{r['src']}.cif")
    shutil.copy(cif[0], dst)
    obj = f"c{i:02d}"
    hot = "+".join(str(h + OFF) for h in HOTS.get(r["src"], []))
    png = os.path.join(OUTDIR, f"cand{i:02d}_{r['src']}.png")
    lines += [
        f"load {dst}, {obj}",
        f"hide everything, {obj}",
        f"show cartoon, {obj}",
        # receptor grey, binder by segment
        f"color grey80, {obj} and resi 141-549",
        f"color skyblue, {obj} and resi 11-140",
        f"color orange, {obj} and resi 1-10",
        f"show sticks, {obj} and resi {hot} and not name N+C+O" if hot else "",
        f"color firebrick, {obj} and resi {hot}" if hot else "",
        f"show spheres, {obj} and resi {hot} and name CA" if hot else "",
        f"set sphere_scale, 0.45, {obj} and resi {hot}" if hot else "",
        f"orient {obj} and resi 1-140",
        f"turn x, -15",
        f"zoom {obj} and resi 1-140, 12",
        f"ray 1100, 850",
        f"png {png}, dpi=150",
        f"disable {obj}",
    ]
    manifest.append({"i": i, "png": png, "cif": dst, **r})

# a close interface view of the best one
if manifest:
    b = manifest[0]
    hot = "+".join(str(h + OFF) for h in HOTS.get(b["src"], []))
    lines += [
        "enable c01", "hide everything, c01", "show cartoon, c01",
        "color grey70, c01 and resi 141-549",
        "color skyblue, c01 and resi 11-140",
        "color orange, c01 and resi 1-10",
        f"show surface, c01 and resi 141-549 within 10 of (c01 and resi 11-140)",
        "set transparency, 0.55",
        f"color firebrick, c01 and resi {hot}" if hot else "",
        f"show sticks, c01 and resi {hot}" if hot else "",
        f"orient c01 and resi {hot}" if hot else "orient c01",
        f"zoom c01 and resi {hot}, 14" if hot else "zoom c01, 10",
        "ray 1300, 1000",
        f"png {os.path.join(OUTDIR, 'interface_best.png')}, dpi=150",
    ]
lines.append("quit")
script = os.path.join(OUTDIR, "render.pml")
open(script, "w").write("\n".join(l for l in lines if l) + "\n")
import json
json.dump(manifest, open(os.path.join(OUTDIR, "manifest.json"), "w"), indent=1)
print(f"wrote {script} for {len(manifest)} structures")
for m in manifest:
    print(f"  cand{m['i']:02d} {m['src']:4s} score {float(m['score']):.3f} "
          f"q={m['q_des']:>3s} motif {float(m['motif']):.2f} "
          f"{m['hits']}/{m['n_total']} sb={m['sb']} hc={m['hc']}")
