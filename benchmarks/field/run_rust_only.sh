#!/usr/bin/env bash
# Run only speeddada + Python benchmarks (R dada2 skipped).
# Usage: ./run_rust_only.sh [THREADS]
set -uo pipefail

THREADS="${1:-16}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORKSPACE="$(cd "$HERE/../.." && pwd)"
OUT_BASE=/tmp/bench_field_out
IN_DIR=/Users/alex/Downloads/raw_data_FIELD

mkdir -p "$OUT_BASE"/{speeddada,python}

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
    echo "  peak RSS:"
    grep -i "maximum resident" "$rss_file" 2>/dev/null || echo "  (rss file missing)"
    echo
}

# 1. speeddada (Rust R binding) — native pseudo via wrap__dada_pseudo
run_with_rss "speeddada" \
    "$OUT_BASE/speeddada/rss.txt" \
    "$OUT_BASE/speeddada/log.txt" \
    env RAYON_NUM_THREADS="$THREADS" \
    Rscript "$HERE/bench_speeddada.R" "$THREADS" "$IN_DIR" "$OUT_BASE/speeddada"

# 2. Python dada2 — native dada_pseudo
run_with_rss "Python dada2" \
    "$OUT_BASE/python/rss.txt" \
    "$OUT_BASE/python/log.txt" \
    env RAYON_NUM_THREADS="$THREADS" \
    "$WORKSPACE/.venv/bin/python" "$HERE/bench_rust.py" \
        --threads "$THREADS" --in-dir "$IN_DIR" --out-dir "$OUT_BASE/python"

echo
echo "==== Comparison (skipping R dada2) ===="
"$WORKSPACE/.venv/bin/python" "$HERE/compare.py"
