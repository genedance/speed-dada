---
name: test-runner
description: Run the test suite and report failures with diagnosis.
---

Run the following in order and report results:
1. `cargo test --workspace 2>&1`
2. `cargo clippy --workspace --all-targets -- -D warnings 2>&1`
3. `cd crates/dada2-py && maturin develop && pytest tests/ -v 2>&1`

For each failure, report: test name, error message, likely cause, suggested fix.
