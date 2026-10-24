import json, glob, os, sys
hdr = ("design", "len", "motif", "maxCAdev", "brk", "clash", "helix", "nSS")
print("%-46s %4s %9s %9s %4s %6s %6s %4s" % hdr)
for j in sorted(glob.glob(os.path.join(sys.argv[1], "*.json"))):
    d = json.load(open(j)); m = d["metrics"]; im = d["diffused_index_map"]
    sc = d["specification"]["extra"]["sampled_contig"]
    binder = [v for k, v in im.items() if k.startswith("B") and v.startswith("A")]
    flanks = sum(int(x[:-1]) for x in sc.split(",/0,")[0].split(",") if x.endswith("P"))
    mp = sorted(int(v[1:]) for v in binder)
    print("%-46s %4d %9s %9.2f %4d %6d %6.2f %4d" % (
        os.path.basename(j).replace("glp1r_pep10_glp1r_", "")[:-5],
        flanks + len(binder), "%d-%d" % (min(mp), max(mp)),
        m["max_ca_deviation"], m["n_chainbreaks"],
        m["n_clashing.interresidue_clashes_w_sidechain"],
        m["helix_fraction"], m["num_ss_elements"]))
