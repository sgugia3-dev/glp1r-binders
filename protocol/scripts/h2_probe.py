import glob, json, os, sys
indir, w = sys.argv[1], sys.argv[2]
cif = sorted(glob.glob(os.path.join(indir, "_cif", "*.cif")))
if not cif:
    print("no staged cif; looking for gz"); sys.exit(2)
c = cif[0]; print("structure:", os.path.basename(c))
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
e = MPNNInferenceEngine(model_type="protein_mpnn", checkpoint_path=w,
                        is_legacy_weights=True, out_directory="/tmp",
                        write_fasta=False, write_structures=False, device="cpu")
inp = {"bias": None, "bias_per_residue": None, "omit": ["CYS"], "omit_per_residue": None,
       "structure_path": c, "name": "probe", "seed": 37,
       "batch_size": 2, "number_of_batches": 1, "temperature": 0.1}
o = e.run(input_dicts=[inp], atom_arrays=None)
print("type(outs):", type(o).__name__, "len:", len(o) if hasattr(o, "__len__") else "n/a")
def walk(x, p="", d=0):
    if d > 3: return
    if isinstance(x, dict):
        for k, v in x.items():
            t = type(v).__name__
            ln = f" len={len(v)}" if hasattr(v, "__len__") and not isinstance(v, str) else ""
            pv = (repr(v)[:90] if isinstance(v, (str, int, float, bool)) else "")
            print(f"  {p}{k}: {t}{ln} {pv}")
            walk(v, p + k + ".", d + 1)
    elif isinstance(x, (list, tuple)) and x:
        print(f"  {p}[0] is {type(x[0]).__name__}")
        walk(x[0], p + "[0].", d + 1)
walk(o)
