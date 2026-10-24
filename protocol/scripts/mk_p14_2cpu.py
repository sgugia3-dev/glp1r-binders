"""2-cpu variant of the redesign pipeline.

Requesting 2 cores instead of 4 makes more GPUs reachable, because the GPU nodes
run out of cores long before they run out of GPUs and Slurm will not oversubscribe.
Boltz's worker count has to come down with it: 2 dataloader workers plus the main
process would oversubscribe a 2-core allocation and slow folding rather than
speed it up.
"""
import os
S = "/lustre/home/gsgueglia/challenge2026/glp1r2026/scripts"
s = open(os.path.join(S, "p14_redesign.sbatch")).read()
subs = [
    ("#SBATCH --cpus-per-task=4", "#SBATCH --cpus-per-task=2"),
    ("#SBATCH --job-name=p14-chg", "#SBATCH --job-name=p14-chg2"),
    ("--num_workers 2 --preprocessing-threads 2",
     "--num_workers 1 --preprocessing-threads 1"),
]
for old, new in subs:
    n = s.count(old)
    assert n == 1, (old, n)
    s = s.replace(old, new)
    print(f"  {old[:52]!r} -> {new[:46]!r}")
hdr = ("#\n# 2 CORES PER TASK: the GPU nodes exhaust cores before GPUs, and Slurm will not\n"
       "# oversubscribe, so a 4-core request leaves reachable GPUs stranded. Boltz drops to\n"
       "# 1 dataloader worker and 1 preprocessing thread to match the smaller allocation.\n")
s = s.replace("#\n# Charge-complementary redesign round", hdr + "#\n# Charge-complementary redesign round", 1)
dst = os.path.join(S, "p14_redesign_2cpu.sbatch")
open(dst, "w").write(s)
os.chmod(dst, 0o755)
print("wrote", dst)
