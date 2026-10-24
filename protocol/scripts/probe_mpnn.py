import inspect
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
print("engine init:", str(inspect.signature(MPNNInferenceEngine.__init__))[:400])
print("run        :", str(inspect.signature(MPNNInferenceEngine.run))[:400])
src = inspect.getsource(MPNNInferenceEngine)
want = ("bias_per_residue", "bias", "omit_per_residue", "designed_residues")
for key in want:
    hits = [l.strip() for l in src.splitlines() if key in l]
    print(f"\n--- {key}  ({len(hits)} lines)")
    for l in hits[:6]:
        print("   ", l[:130])
