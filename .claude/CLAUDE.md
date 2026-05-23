# dada2_rust — CLAUDE.md

## Project
High-performance Rust reimplementation of the DADA2 amplicon sequencing pipeline,
with Python bindings via PyO3/maturin.

## Build
```bash
cargo build --workspace          # debug
cargo build --release --workspace
cd crates/dada2-py && maturin develop  # Python wheel (dev)
```

## Test
```bash
cargo test --workspace
cd crates/dada2-py && pytest tests/
```

## Lint
```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
```

## Key constraints
- No `unsafe{}` without `# SAFETY:` comment
- No `.unwrap()` or `.expect()` in library code
- Every public function/struct must have a `///` doc-comment
- Each `.rs` file must stay under 500 lines
- `clippy::pedantic` enabled in dada2-core
- All CPU-intensive loops must use `rayon::par_iter()`
- Deterministic output: RNG seed configurable via `seed: u64` field (default 42)
