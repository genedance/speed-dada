"""
Parity test: compare Rust pipeline output against saved R dada2 output.

Usage
-----
    python compare_r_output.py <rust_output.json>

The script loads:
  - tests/parity/fixtures/r_output.json  — reference ASV table from R dada2
  - <rust_output.json>                   — Rust pipeline output (same format)

Asserts:
  - Jaccard similarity of ASV sets >= 0.95
  - Per-ASV abundance Pearson r >= 0.99 (for ASVs present in both)
  - No Rust ASV with abundance > 10 absent from R output
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


def load_asvs(path: Path) -> dict[str, int]:
    """Load ASV -> abundance mapping from a JSON file."""
    data = json.loads(path.read_text())
    return {entry["sequence"]: entry["abundance"] for entry in data}


def jaccard(set_a: set, set_b: set) -> float:
    union = set_a | set_b
    if not union:
        return 1.0
    return len(set_a & set_b) / len(union)


def pearson_r(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n < 2:
        return 1.0
    mx, my = sum(x) / n, sum(y) / n
    num = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    denom = math.sqrt(
        sum((xi - mx) ** 2 for xi in x) * sum((yi - my) ** 2 for yi in y)
    )
    return num / denom if denom > 0 else 1.0


def main(rust_json: Path) -> None:
    fixtures_dir = Path(__file__).parent / "fixtures"
    r_asvs = load_asvs(fixtures_dir / "r_output.json")
    rust_asvs = load_asvs(rust_json)

    r_seqs = set(r_asvs)
    rust_seqs = set(rust_asvs)

    jac = jaccard(r_seqs, rust_seqs)
    print(f"Jaccard similarity : {jac:.4f}")
    assert jac >= 0.95, f"Jaccard {jac:.4f} < 0.95"

    shared = r_seqs & rust_seqs
    if len(shared) >= 2:
        r_ab = [r_asvs[s] for s in shared]
        rust_ab = [rust_asvs[s] for s in shared]
        r_val = pearson_r(r_ab, rust_ab)
        print(f"Abundance Pearson r: {r_val:.4f}")
        assert r_val >= 0.99, f"Pearson r {r_val:.4f} < 0.99"

    false_positives = [
        (s, rust_asvs[s]) for s in rust_seqs - r_seqs if rust_asvs[s] > 10
    ]
    if false_positives:
        for seq, ab in false_positives:
            print(f"  FALSE POSITIVE: {seq[:20]}... abundance={ab}")
        raise AssertionError(f"{len(false_positives)} Rust ASV(s) with abundance > 10 absent from R output")

    print("All parity checks passed")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <rust_output.json>", file=sys.stderr)
        sys.exit(1)
    main(Path(sys.argv[1]))
