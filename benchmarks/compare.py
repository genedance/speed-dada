"""Compare Rust vs R dada2 ASV output: Jaccard similarity and abundance correlation."""
import json, math, sys
from pathlib import Path

OUT = Path("/tmp/bench_out")


def load(path):
    data = json.loads(Path(path).read_text())
    return {e["sequence"]: e["abundance"] for e in data}


def jaccard(a, b):
    u = a | b
    return len(a & b) / len(u) if u else 1.0


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return 1.0
    mx, my = sum(xs)/n, sum(ys)/n
    num = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    denom = math.sqrt(sum((x-mx)**2 for x in xs) * sum((y-my)**2 for y in ys))
    return num / denom if denom else 1.0


rust = load(OUT / "rust_output.json")
r    = load(OUT / "r_output.json")

print("=== ASV count ===")
print(f"  Rust : {len(rust)}")
print(f"  R    : {len(r)}")

jac = jaccard(set(rust), set(r))
print(f"\n=== Jaccard similarity of ASV sets: {jac:.4f} ===")

shared = set(rust) & set(r)
if shared:
    rv = [rust[s] for s in shared]
    rv2 = [r[s] for s in shared]
    pr = pearson(rv, rv2)
    print(f"=== Abundance Pearson r (shared ASVs): {pr:.4f} ===")

only_rust = set(rust) - set(r)
only_r    = set(r)    - set(rust)
if only_rust:
    print(f"\nOnly in Rust ({len(only_rust)}):")
    for s in sorted(only_rust, key=lambda x: -rust[x])[:5]:
        print(f"  {s[:40]}...  ab={rust[s]}")
if only_r:
    print(f"\nOnly in R ({len(only_r)}):")
    for s in sorted(only_r, key=lambda x: -r[x])[:5]:
        print(f"  {s[:40]}...  ab={r[s]}")

print("\n=== Top ASVs ===")
print(f"{'Sequence[:35]':<37} {'Rust':>8}  {'R':>8}")
all_seqs = sorted(set(rust) | set(r), key=lambda s: -(rust.get(s,0)+r.get(s,0)))
for s in all_seqs[:10]:
    print(f"  {s[:35]:<35}  {rust.get(s,'-'):>8}  {r.get(s,'-'):>8}")
