# Skill: rust-cargo

## Invoke with
"use the rust-cargo skill"

## Cargo workspace patterns
- `cargo build --workspace` builds all crates
- `cargo test --workspace` tests all crates
- `cargo clippy --workspace --all-targets -- -D warnings`
- `cargo fmt --all -- --check` (check only) / `cargo fmt --all` (format)
- For PyO3 crates: `cd crates/dada2-py && maturin develop`

## Common Cargo.toml patterns
```toml
[workspace.dependencies]
# Define versions once; use `.workspace = true` in member crates
```

## Dependency resolution
- Run `cargo update` to refresh lock file
- Run `cargo tree` to inspect dependency graph
