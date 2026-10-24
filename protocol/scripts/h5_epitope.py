#!/usr/bin/env python
"""Did the DESIGNED extension actually engage the hotspot patch?

Motif RMSD only says the fixed segment stayed put. It says nothing about where the
designed part went, and this project has repeatedly seen high-confidence designs bind a
different face entirely (one scored ipSAE 0.79 while binding 11-20 A from its intended
hotspots). So epitope adherence is measured geometrically, per design.

Measured per refold:
  - min heavy-atom distance from the DESIGNED region (binder 11..L, motif excluded) to
    each requested hotspot residue
  - how many of the 5 hotspots are within 5 A
  - the same for the motif region (1-10) as a cross-check that it still sits in its pocket
Chains are identified by length, never by label.
"""
import argparse, csv, glob, json, math, os, sys

HOT = [24, 25, 27, 42]   # Ser26 dropped: 0/30 designs reached it (median 6.8 A)
CUT = 5.0


def read_atoms(path):
    hdr, rows, inl = [], [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            hdr.append(line.strip().split(".")[1]); inl = True; continue
        if inl and line.startswith(("ATOM", "HETATM")): rows.append(line.split())
        elif inl and line.startswith("#"): break
    i = {h: k for k, h in enumerate(hdr)}
    need = ("auth_asym_id", "auth_seq_id", "label_atom_id", "Cartn_x", "Cartn_y", "Cartn_z")
    if any(n not in i for n in need): return {}
    out = {}
    for r in rows:
        if len(r) <= max(i[n] for n in need): continue
        if r[i["label_atom_id"]].startswith("H"): continue      # heavy atoms only
        try: n = int(r[i["auth_seq_id"]])
        except ValueError: continue
        out.setdefault(r[i["auth_asym_id"]], {}).setdefault(n, []).append(
            (float(r[i["Cartn_x"]]), float(r[i["Cartn_y"]]), float(r[i["Cartn_z"]])))
    return out


def mindist(a, b):
    best = 1e9
    for p in a:
        for q in b:
            if abs(p[0]-q[0]) > best or abs(p[1]-q[1]) > best: continue
            d = math.dist(p, q)
            if d < best: best = d
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", required=True)
    ap.add_argument("--manifest", required=True)
    # the p14 round mixes backbones from campaigns with DIFFERENT hotspot sets in
    # one prediction directory, so each set must be scored over only its own
    # designs. Without this the glob would score every design four times.
    ap.add_argument("--name-prefix", default=None,
                    help="only score designs whose directory name starts with this")
    ap.add_argument("--dev-csv", required=True, help="h3 deviation csv, to join motif RMSD")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--hotspots", default=None,
                    help="comma-separated target residues; default 24,25,27,42")
    a = ap.parse_args()

    global HOT
    if a.hotspots:
        HOT = [int(x) for x in a.hotspots.split(",")]
    print(f"scoring against hotspots {HOT}")
    man = {y["name"]: y for y in json.load(open(a.manifest))["yamls"]}
    dev = {r["name"]: r for r in csv.DictReader(open(a.dev_csv))}

    rows = []
    cifs = sorted(glob.glob(os.path.join(a.pred_dir, "*", "*_model_0.cif")))
    if a.name_prefix:
        before = len(cifs)
        cifs = [c for c in cifs
                if os.path.basename(os.path.dirname(c)).startswith(a.name_prefix)]
        print(f"name-prefix {a.name_prefix!r}: {len(cifs)}/{before} designs")
    for cif in cifs:
        name = os.path.basename(os.path.dirname(cif))
        rec = man.get(name)
        if not rec: continue
        ch = read_atoms(cif)
        if len(ch) < 2: continue
        tgt = max(ch, key=lambda c: len(ch[c]))
        bnd = [c for c in ch if c != tgt][0]
        motif = set(rec["motif"])
        blen = rec["binder_len"]
        des_atoms = [p for k, v in ch[bnd].items() if k not in motif for p in v]
        mot_atoms = [p for k, v in ch[bnd].items() if k in motif for p in v]
        if not des_atoms: continue
        d_des, d_mot = {}, {}
        for h in HOT:
            if h not in ch[tgt]: continue
            d_des[h] = mindist(des_atoms, ch[tgt][h])
            d_mot[h] = mindist(mot_atoms, ch[tgt][h]) if mot_atoms else 99.0
        nhit = sum(1 for h in d_des if d_des[h] <= CUT)
        row = {"name": name, "binder_len": blen,
               "motif_rmsd": dev.get(name, {}).get("motif_rmsd_vs_af3", ""),
               "receptor_rmsd": dev.get(name, {}).get("receptor_rmsd", ""),
               "n_hotspots_hit": nhit,
               "min_hotspot_dist": round(min(d_des.values()), 2) if d_des else "",
               "designed_engages": nhit >= 3}
        for h in HOT:
            row[f"d_des_A{h}"] = round(d_des.get(h, 99.0), 2)
        rows.append(row)

    if not rows: print("FATAL: nothing measured"); return 3
    rows.sort(key=lambda r: (-r["n_hotspots_hit"], r["min_hotspot_dist"]))
    cols = list(rows[0].keys())
    with open(a.out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows: fh.write(",".join(str(r[c]) for c in cols) + "\n")

    import statistics as st
    print(f"measured {len(rows)} refolds, hotspots {HOT}, contact cut {CUT} A\n")
    print(f"{'design':34s} {'len':>4s} {'motifRMSD':>10s} {'hits':>5s} {'minDist':>8s}")
    for r in rows[:12]:
        print(f"{r['name'].replace('glp1r_ecd_glp1r_','')[:34]:34s} {r['binder_len']:4d} "
              f"{str(r['motif_rmsd']):>10s} {r['n_hotspots_hit']:5d} {r['min_hotspot_dist']:8}")
    hits = [r["n_hotspots_hit"] for r in rows]
    print(f"\nhotspots contacted per design: median {st.median(hits)}, max {max(hits)}")
    print(f"designs engaging >=3 of {len(HOT)}     : {sum(1 for r in rows if r['designed_engages'])}/{len(rows)}")
    print(f"designs engaging >=1          : {sum(1 for h in hits if h >= 1)}/{len(rows)}")
    print(f"designs engaging 0            : {sum(1 for h in hits if h == 0)}/{len(rows)}")
    print(f"\nwrote {a.out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
