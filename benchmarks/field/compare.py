#!/usr/bin/env python3
"""Three-way comparison of bench_r.R, bench_dada2rs.R, bench_rust.py."""
import json
import math
import sys
from pathlib import Path

BASE = Path("/tmp/bench_field_out")
TOOLS = [
    ("R dada2",      BASE / "r"        / "r_output.json",        BASE / "r"        / "rss.txt"),
    ("dada2rs",      BASE / "dada2rs"  / "dada2rs_output.json",  BASE / "dada2rs"  / "rss.txt"),
    ("Python dada2", BASE / "python"   / "rust_output.json",     BASE / "python"   / "rss.txt"),
]

STAGES = ["filter_ms", "learn_errors_ms", "derep_ms", "dada_ms", "merge_ms", "chimera_ms"]
LABELS = ["filter",    "learn_err",       "derep",    "dada",    "merge",    "chimera"]


def jaccard(a, b):
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    num   = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denom = math.sqrt(sum((x - mx) ** 2 for x in xs) *
                      sum((y - my) ** 2 for y in ys))
    return num / denom if denom else float("nan")


def parse_rss(path):
    """macOS /usr/bin/time -l output: parse 'maximum resident set size' bytes."""
    if not path.exists():
        return None
    for line in path.read_text().splitlines():
        line = line.strip()
        if "maximum resident set size" in line:
            return int(line.split()[0])
    return None


def fmt_ms(ms):
    if ms is None:
        return "    -"
    if ms < 1000:
        return f"{ms:>7.1f} ms"
    if ms < 60_000:
        return f"{ms/1000:>7.2f} s "
    return f"{ms/60000:>7.2f} m "


def fmt_mb(bytes_):
    if bytes_ is None:
        return "    -"
    return f"{bytes_ / 1024 / 1024:>8.1f} MB"


def main():
    results = []
    for name, json_path, rss_path in TOOLS:
        if not json_path.exists():
            print(f"  WARNING: {json_path} not found — skipping {name}")
            continue
        data = json.loads(json_path.read_text())
        asvs = {e["sequence"]: e["abundance"] for e in data["asvs"]}
        rss  = parse_rss(rss_path)
        results.append((name, data, asvs, rss))

    if not results:
        print("No benchmark outputs found.")
        sys.exit(1)

    print("=" * 80)
    print("  STAGE TIMINGS")
    print(f"  {'stage':<12}", end="")
    for name, _, _, _ in results:
        print(f"  {name:>16}", end="")
    print()
    print("-" * 80)
    for stage, label in zip(STAGES, LABELS):
        print(f"  {label:<12}", end="")
        for _, data, _, _ in results:
            ms = data.get("stages", {}).get(stage)
            print(f"  {fmt_ms(ms):>16}", end="")
        print()
    print("-" * 80)
    print(f"  {'TOTAL':<12}", end="")
    for _, data, _, _ in results:
        print(f"  {fmt_ms(data.get('total_ms')):>16}", end="")
    print()

    print()
    print("=" * 80)
    print("  RESOURCE USE  (peak RSS measured by /usr/bin/time -l)")
    for name, _, _, rss in results:
        print(f"  {name:<16}  {fmt_mb(rss)}")

    # Speedup vs R dada2
    if results[0][0] == "R dada2":
        print()
        print("=" * 80)
        print("  SPEEDUP vs R dada2 (total_ms)")
        ref_total = results[0][1].get("total_ms")
        for name, data, _, _ in results:
            total = data.get("total_ms")
            if ref_total and total:
                print(f"  {name:<16}  {ref_total / total:>6.1f}×")

    print()
    print("=" * 80)
    print("  ASV COUNTS")
    for name, data, asvs, _ in results:
        n_before = data.get("n_asvs_before_chimera", "?")
        n_after  = data.get("n_asvs_after_chimera", len(asvs))
        total    = data.get("total_abundance", sum(asvs.values()))
        print(f"  {name:<16}  pre-chimera={n_before:>5}  "
              f"post-chimera={n_after:>5}  reads_in_asvs={total:,}")

    print()
    print("=" * 80)
    print("  PER-SAMPLE READ COUNTS")
    sample_names = [s["sample"] for s in results[0][1]["samples"]]
    print(f"  {'sample':<22}", end="")
    for name, _, _, _ in results:
        print(f"  {name[:14]:>14}", end="")
    print()
    for i, sn in enumerate(sample_names):
        print(f"  {sn:<22}", end="")
        for _, data, _, _ in results:
            s = data["samples"][i] if i < len(data["samples"]) else None
            cell = f"{s['reads_in']:>5}→{s['reads_out']:<6}" if s else "?"
            print(f"  {cell:>14}", end="")
        print()

    print()
    print("=" * 80)
    print("  PAIRWISE ASV SIMILARITY")
    for i in range(len(results)):
        for j in range(i + 1, len(results)):
            na, _, aa, _ = results[i]
            nb, _, ab, _ = results[j]
            jac    = jaccard(set(aa), set(ab))
            shared = set(aa) & set(ab)
            if shared:
                r = pearson([aa[s] for s in shared], [ab[s] for s in shared])
                print(f"  {na} vs {nb}:  "
                      f"Jaccard={jac:.4f}  Pearson_r(abund)={r:.4f}  "
                      f"shared={len(shared)}/{max(len(aa),len(ab))}")
            else:
                print(f"  {na} vs {nb}:  Jaccard={jac:.4f}  (no shared ASVs)")


if __name__ == "__main__":
    main()
