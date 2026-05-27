#!/usr/bin/env bash
# Run all three benchmarks, capturing peak RSS via /usr/bin/time -l.
# Usage: ./run_all.sh [THREADS]
set -uo pipefail
# Note: no -e — one tool's failure should not stop the others.

THREADS="${1:-16}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORKSPACE="$(cd "$HERE/../.." && pwd)"
OUT_BASE=/tmp/bench_field_out
IN_DIR=/Users/alex/Downloads/raw_data_FIELD

mkdir -p "$OUT_BASE"/{r,speeddada,python}

run_with_rss() {
    local label="$1"; shift
    local rss_file="$1"; shift
    local log_file="$1"; shift
    echo "==== $label ===="
    # /usr/bin/time -l writes resource summary to stderr on macOS.
    if /usr/bin/time -l -o "$rss_file" "$@" 2>&1 | tee "$log_file"; then
        echo "  exit: 0"
    else
        echo "  exit: non-zero — continuing with remaining tools"
    fi
    echo "  peak RSS:"
    grep -i "maximum resident" "$rss_file" 2>/dev/null || echo "  (rss file missing)"
    echo
}

# 1. R dada2 (Bioconductor)
run_with_rss "R dada2" \
    "$OUT_BASE/r/rss.txt" \
    "$OUT_BASE/r/log.txt" \
    Rscript "$HERE/bench_r.R" "$THREADS" "$IN_DIR" "$OUT_BASE/r"

# 2. speeddada (Rust R binding)
run_with_rss "speeddada" \
    "$OUT_BASE/speeddada/rss.txt" \
    "$OUT_BASE/speeddada/log.txt" \
    env RAYON_NUM_THREADS="$THREADS" \
    Rscript "$HERE/bench_speeddada.R" "$THREADS" "$IN_DIR" "$OUT_BASE/speeddada"

# 3. Python dada2 (Rust Python binding)
run_with_rss "Python dada2" \
    "$OUT_BASE/python/rss.txt" \
    "$OUT_BASE/python/log.txt" \
    env RAYON_NUM_THREADS="$THREADS" \
    "$WORKSPACE/.venv/bin/python" "$HERE/bench_rust.py" \
        --threads "$THREADS" --in-dir "$IN_DIR" --out-dir "$OUT_BASE/python"

echo
echo "==== Comparison ===="
"$WORKSPACE/.venv/bin/python" "$HERE/compare.py"
