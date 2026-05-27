# Contributing to speed-dada

Thanks for your interest in contributing! This document covers the
expectations for code, tests, and review.

## Code rules

The workspace enforces a small set of conventions — most are checked by CI:

| Rule | Where |
|---|---|
| No `.unwrap()` / `.expect()` in library code (tests excepted) | `crates/speeddada-core`, `crates/speeddada-py`, `crates/speeddada-r` |
| Every public function / struct has a `///` doc comment | All crates |
| No `unsafe { }` without a `// SAFETY:` comment explaining the invariant | All crates |
| Each `.rs` file stays under 500 lines | Hard cap — split into sub-modules |
| `clippy::pedantic` enabled (warning, not error) | `speeddada-core` |
| CPU-intensive loops use `rayon` | All crates |
| Deterministic output: RNG seed configurable via `seed: u64` (default 42) | Wherever randomness appears |

## Build & test

```bash
# Rust
cargo fmt --all -- --check
cargo clippy --workspace --exclude speeddada-py --exclude speeddada-r --all-targets -- -D warnings
cargo test  --workspace --exclude speeddada-py --exclude speeddada-r

# Python (extension module — needs maturin)
cd crates/speeddada-py
maturin develop
pytest tests/

# R
cd r-package/SpeedDada
R -e "roxygen2::roxygenise('.')"      # regenerate man/*.Rd from #' comments
R -e 'install.packages(".", repos = NULL, type = "source")'
R CMD check --no-manual --as-cran .
```

`speeddada-py` is excluded from `cargo build` / `cargo test` because its
`cdylib` needs the Python headers + a custom linker config that maturin
provides. The compiled wheel is validated separately in
`.github/workflows/wheels.yml`.

## Pull requests

1. Open a small, focused PR with a clear title.
2. Run `cargo fmt --all` before committing.
3. Add tests for any new behaviour.
4. Update `CHANGELOG.md` under the **Unreleased** section.
5. CI (rust + r-cmd-check + wheels) must pass before merge.

## Reporting bugs

Please open an issue at
<https://github.com/Genedance/speed-dada/issues> with:

* A minimal reproduction (FASTQ snippet or simulated data is ideal).
* Your platform (`uname -srm`, R / Python / cargo version).
* The full error message and a stack trace if available.
