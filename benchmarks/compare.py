"""Compare ASV output and timing across R dada2, dada2rs, and Python dada2."""
import json, math
from pathlib import Path

OUT = Path("/tmp/bench_out")

TOOLS = [
    ("R dada2",     OUT / "r_output.json"),
    ("dada2rs",     OUT / "dada2rs_output.json"),
    ("Python dada2", OUT / "rust_output.json"),
]

STAGES = ["filter_ms", "learn_errors_ms", "derep_ms", "dada_ms", "merge_ms", "chimera_ms"]
LABELS = ["filter", "learn_err", "derep", "dada", "merge", "chimera"]


def load(path):
    data = json.loads(Path(path).read_text())
    asvs = {e["sequence"]: e["abundance"] for e in data["asvs"]}
    return data, asvs


def jaccard(a, b):
    u = a | b
    return len(a & b) / len(u) if u else 1.0


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    num   = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denom = math.sqrt(sum((x - mx) ** 2 for x in xs) *
                      sum((y - my) ** 2 for y in ys))
    return num / denom if denom else float("nan")


# ── Load ──────────────────────────────────────────────────────────────────────

results = []
for name, path in TOOLS:
    if not path.exists():
        print(f"  WARNING: {path} not found — skipping {name}")
        continue
    data, asvs = load(path)
    results.append((name, data, asvs))

if not results:
    print("No benchmark outputs found in /tmp/bench_out/")
    raise SystemExit(1)

# ── Timing table ──────────────────────────────────────────────────────────────

print("=" * 72)
print("  TIMING  (ms)")
print(f"  {'Stage':<14}", end="")
for name, _, _ in results:
    print(f"  {name:>14}", end="")
print()
print("-" * 72)

for stage, label in zip(STAGES, LABELS):
    print(f"  {label:<14}", end="")
    for _, data, _ in results:
        ms = data.get("stages", {}).get(stage, "-")
        print(f"  {ms:>14}", end="")
    print()

print("-" * 72)
print(f"  {'TOTAL':<14}", end="")
for _, data, _ in results:
    print(f"  {data.get('total_ms', '-'):>14}", end="")
print()

# ── ASV counts ────────────────────────────────────────────────────────────────

print()
print("=" * 72)
print("  ASV COUNTS")
for name, _, asvs in results:
    print(f"  {name:<16}  {len(asvs)} ASVs   "
          f"total abundance {sum(asvs.values()):,}")

# ── Pairwise Jaccard & abundance correlation ──────────────────────────────────

print()
print("=" * 72)
print("  PAIRWISE SIMILARITY")
for i in range(len(results)):
    for j in range(i + 1, len(results)):
        na, _, aa = results[i]
        nb, _, ab = results[j]
        jac    = jaccard(set(aa), set(ab))
        shared = set(aa) & set(ab)
        if shared:
            r = pearson([aa[s] for s in shared], [ab[s] for s in shared])
            print(f"  {na} vs {nb}:  "
                  f"Jaccard={jac:.4f}  Pearson_r(abund)={r:.4f}  "
                  f"shared={len(shared)}/{max(len(aa),len(ab))}")
        else:
            print(f"  {na} vs {nb}:  Jaccard={jac:.4f}  (no shared ASVs)")

# ── Top ASVs ──────────────────────────────────────────────────────────────────

print()
print("=" * 72)
print("  TOP ASVs")
all_seqs = sorted(
    set().union(*(set(asvs) for _, _, asvs in results)),
    key=lambda s: -sum(asvs.get(s, 0) for _, _, asvs in results)
)
header = f"  {'Seq[:32]':<34}"
for name, _, _ in results:
    header += f"  {name[:12]:>12}"
print(header)
print("-" * 72)
for s in all_seqs[:12]:
    row = f"  {s[:32]:<34}"
    for _, _, asvs in results:
        v = asvs.get(s, "-")
        row += f"  {v:>12}"
    print(row)
