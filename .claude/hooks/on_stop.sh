#!/usr/bin/env bash
# Fires asynchronously when Claude finishes a session.
set -euo pipefail

cd /home/fromage/dada2_rust
LOG=".claude/hooks/session.log"
mkdir -p "$(dirname "$LOG")"

echo "=== $(date -u +%FT%TZ) SESSION END ===" >> "$LOG"
{
    cargo fmt  --all -- --check  && echo "fmt ok"    || echo "fmt FAILED"
    cargo clippy --workspace -- -D warnings \
                              && echo "clippy ok" || echo "clippy FAILED"
    cargo test --workspace   && echo "tests ok"  || echo "tests FAILED"
} >> "$LOG" 2>&1

echo "Check log: $LOG"
