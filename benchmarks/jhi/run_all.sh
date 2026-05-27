#!/usr/bin/env bash
# Run all four benchmarks (R dada2, speeddada, Python dada2, rust-only CLI)
# on the 3 JHI paired samples, capturing peak RSS via /usr/bin/time -l.
# Usage: ./run_all.sh [THREADS] [IN_DIR]
set -uo pipefail

THREADS="${1:-16}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORKSPACE="$(cd "$HERE/../.." && pwd)"
OUT_BASE=/tmp/bench_jhi_out
IN_DIR="${2:-/Users/alex/Library/CloudStorage/OneDrive-Genedance/Genedance Corp - Documents/customer_data/COUSIN/BZ20260129886761-799F-1193R-Data-20260316/rawData}"

mkdir -p "$OUT_BASE"/{r,speeddada,python,rust_only}

run_with_rss() {
    local label="$1"; shift
    local rss_file="$1"; shift
    local log_file="$1"; shift
    echo "==== $label ===="
    if /usr/bin/time -l -o "$rss_file" "$@" 2>&1 | tee "$log_file"; then
        echo "  exit: 0"
    else
        echo "  exit: non-zero"
    fi
    grep -i "maximum resident" "$rss_file" 2>/dev/null || echo "  (rss missing)"
    echo
}

# 1. R dada2 (Bioconductor reference)
run_with_rss "R dada2" \
    "$OUT_BASE/r/rss.txt" \
    "$OUT_BASE/r/log.txt" \
    Rscript "$HERE/bench_r.R" "$THREADS" "$IN_DIR" "$OUT_BASE/r"

# 2. speeddada (Rust core via R binding)
run_with_rss "speeddada" \
    "$OUT_BASE/speeddada/rss.txt" \
    "$OUT_BASE/speeddada/log.txt" \
    env RAYON_NUM_THREADS="$THREADS" \
    Rscript "$HERE/bench_speeddada.R" "$THREADS" "$IN_DIR" "$OUT_BASE/speeddada"

# 3. Python dada2 (Rust core via Python binding)
run_with_rss "Python dada2" \
    "$OUT_BASE/python/rss.txt" \
    "$OUT_BASE/python/log.txt" \
    env RAYON_NUM_THREADS="$THREADS" \
    "$WORKSPACE/.venv/bin/python" "$HERE/bench_rust.py" \
        --threads "$THREADS" --in-dir "$IN_DIR" --out-dir "$OUT_BASE/python"

# 4. Rust-only CLI (no bindings)
run_with_rss "rust-only CLI" \
    "$OUT_BASE/rust_only/rss.txt" \
    "$OUT_BASE/rust_only/log.txt" \
    env RAYON_NUM_THREADS="$THREADS" \
    "$WORKSPACE/target/release/dada2-bench" \
        --threads "$THREADS" \
        --in-dir "$IN_DIR" \
        --out-dir "$OUT_BASE/rust_only" \
        --samples "JHI-2025-Q1-A-004,JHI-2025-Q1-A-009,JHI-2025-Q1-A-010"

echo
echo "==== Comparison ===="
"$WORKSPACE/.venv/bin/python" "$HERE/compare.py"
