import glob, os, sys
indir, w = sys.argv[1], sys.argv[2]
c = sorted(glob.glob(os.path.join(indir, "_cif", "*.cif")))[0]
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
e = MPNNInferenceEngine(model_type="protein_mpnn", checkpoint_path=w,
                        is_legacy_weights=True, out_directory="/tmp",
                        write_fasta=False, write_structures=False, device="cpu")
inp = {"bias": None, "bias_per_residue": None, "omit": ["CYS"], "omit_per_residue": None,
       "structure_path": c, "name": "probe", "seed": 37,
       "batch_size": 2, "number_of_batches": 1, "temperature": 0.1}
o = e.run(input_dicts=[inp], atom_arrays=None)[0]
print("class:", type(o).__module__ + "." + type(o).__name__)
pub = [k for k in dir(o) if not k.startswith("_")]
print("attrs:", pub)
for k in pub:
    try: v = getattr(o, k)
    except Exception as ex: print(f"  {k}: <err {ex}>"); continue
    if callable(v): continue
    t = type(v).__name__
    ln = f" len={len(v)}" if hasattr(v, "__len__") and not isinstance(v, str) else ""
    s = repr(v)
    print(f"  {k}: {t}{ln} = {s[:160]}")
