"""Add patch_ncontacts to classic RFdiffusion, additively.

interface_ncontacts maximises contacts to the ENTIRE target. With a 409-residue
receptor and the binder anchored in the TMD core, that mostly inflates contacts
against the large surface near the anchor rather than at the hotspot patch, so it
does not do what 'maximise contacts on the 86-89 region' needs. This restricts
the same PLUMED-style switching function to a contiguous target residue range.

The potential-string parser casts every non-'type' argument with float(), and
float() accepts underscores between digits, so a list argument like
patch:86_87_88_89 would silently become 86878889.0. Hence integer range bounds.
"""
import os, re, shutil

P = "/lustre/home/gsgueglia/diff_test/RFdiffusion/rfdiffusion/potentials/potentials.py"
bak = P + ".bak_patchpot"
if not os.path.exists(bak):
    shutil.copy(P, bak); print("backup ->", bak)
s = open(P).read()
if "class patch_ncontacts" in s:
    print("patch_ncontacts already present; leaving the file alone")
    raise SystemExit(0)

NEW = '''

class patch_ncontacts(Potential):
    \'\'\'
        Maximise the number of contacts between the binder and a CONTIGUOUS PATCH
        of the target, rather than the whole target surface as interface_ncontacts
        does. Same switching function as interface_ncontacts.

        patch_lo / patch_hi are 1-based residue numbers within the target segment,
        in the order the target appears in the contig. They are passed as floats
        because the potential-string parser casts every argument with float().
    \'\'\'

    def __init__(self, binderlen, patch_lo=1, patch_hi=1, weight=1, r_0=8, d_0=6):
        self.binderlen = int(binderlen)
        self.patch_lo  = int(patch_lo)
        self.patch_hi  = int(patch_hi)
        self.r_0       = r_0
        self.d_0       = d_0
        self.weight    = weight

    def compute(self, xyz):
        L = xyz.shape[0]
        # target residues sit after the binder; patch is 1-based within the target
        lo = self.binderlen + (self.patch_lo - 1)
        hi = self.binderlen + (self.patch_hi - 1)
        lo = max(lo, self.binderlen)
        hi = min(hi, L - 1)
        if hi < lo:
            raise ValueError(
                f"patch_ncontacts: empty patch (binderlen={self.binderlen}, "
                f"patch={self.patch_lo}-{self.patch_hi}, L={L})")
        Ca_b = xyz[:self.binderlen, 1]
        Ca_p = xyz[lo:hi + 1, 1]
        dgram = torch.cdist(Ca_b[None, ...].contiguous(),
                            Ca_p[None, ...].contiguous(), p=2)
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0, 6)
        denominator = torch.pow(divide_by_r_0, 12)
        ncontacts = (1 - numerator) / (1 - denominator)
        ncontacts = ncontacts.sum()
        print("PATCH CONTACTS:", ncontacts.sum(), "patch idx", lo, hi)
        return self.weight * ncontacts
'''

# append the class just before the registry
anchor = "implemented_potentials = {"
i = s.index(anchor)
s = s[:i] + NEW.lstrip("\n") + "\n" + s[i:]

# register it
s = s.replace("implemented_potentials = { 'monomer_ROG':          monomer_ROG,",
              "implemented_potentials = { 'patch_ncontacts':      patch_ncontacts,\n"
              "                           'monomer_ROG':          monomer_ROG,", 1)
assert "'patch_ncontacts':      patch_ncontacts" in s, "registry insert failed"

# it needs binderlen injected at runtime
m = re.search(r"require_binderlen\s*=\s*\{", s)
s = s[:m.end()] + "'patch_ncontacts',\n                         " + s[m.end():]

open(P, "w").write(s)
import ast
ast.parse(open(P).read())
print("added + registered patch_ncontacts; file parses")
for pat in ("class patch_ncontacts", "'patch_ncontacts':      patch_ncontacts",
            "'patch_ncontacts',"):
    print(f"  present: {pat!r} -> {pat in open(P).read()}")
