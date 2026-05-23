#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "=== cargo fmt ==="
cargo fmt --all -- --check

echo "=== cargo clippy ==="
cargo clippy --workspace --all-targets -- -D warnings

echo "=== cargo test ==="
cargo test --workspace

echo "=== cargo build --release ==="
cargo build --release --workspace

echo "=== maturin develop + pytest ==="
cd crates/dada2-py
maturin develop
pytest tests/

echo "All checks passed"
