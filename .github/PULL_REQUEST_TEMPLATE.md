## Summary

<!-- 1–3 bullet points describing what this PR does. -->

## Type of change

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing scripts to fail)
- [ ] Documentation / packaging only

## Checklist

- [ ] `cargo fmt --all -- --check` passes
- [ ] `cargo clippy --workspace --exclude speeddada-py --all-targets -- -D warnings` passes
- [ ] `cargo test --workspace --exclude speeddada-py` passes
- [ ] R: `R CMD check r-package/SpeedDada` passes (if R is touched)
- [ ] Python: `pytest crates/speeddada-py/tests` passes (if Python is touched)
- [ ] `CHANGELOG.md` updated
- [ ] Public items have `///` / `#'` doc comments
- [ ] No `.unwrap()` / `.expect()` in library code (tests excepted)
