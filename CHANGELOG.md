# Changelog

All notable changes to this project are documented in this file. The format
is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.99.0 — 2026-05-27

Initial public release under the **speed-dada** brand (R: `SpeedDada`,
Python: `speeddada`).

Distributed via **GitHub Releases** for this version: pre-built Python
wheels (Linux x86_64 + aarch64, macOS x86_64 + arm64, Windows x86_64) and
the R source tarball are attached to each tagged release. PyPI and
Bioconductor submission will follow once the GitHub channel has had a
testing window.

### Added

- Cross-platform build scripts (`configure`, `configure.win`,
  `src/Makevars.win.in`) so the R package installs on Linux, macOS, and
  Windows out of the box, including aarch64/Raspberry Pi.
- `inst/extdata/` FASTQ fixtures + a Bioconductor-style vignette
  (`SpeedDada-pipeline.Rmd`) demonstrating the full paired-end workflow.
- `man/*.Rd` documentation generated from roxygen2; every exported function
  has a runnable `@examples` block.
- `tests/testthat/` smoke tests for every exported function.
- Python mixed Rust/Python layout under `python/speeddada/` with type stubs
  (`.pyi`) and the PEP 561 `py.typed` marker.
- GitHub Actions workflows: `wheels.yml` (cibuildwheel-equivalent matrix
  via PyO3/maturin-action across Linux x86_64+aarch64, macOS x86_64+arm64,
  Windows x86_64), `rust.yml`, `r-cmd-check.yml`.
- Citation files (`CITATION.cff` at repo root, `inst/CITATION` in the R
  package), top-level `LICENSE`, `CONTRIBUTING.md`.

### Changed

- Rebranded all crates and packages under the **speed-dada** name with
  Alexandre Jousset / Genedance GmbH as sole author/copyright holder.
- Split the monolithic `crates/speeddada-core/src/dada.rs` (754 lines) into
  `dada.rs` + `dada_scoring.rs` + `dada_pool.rs`, each well under the 500-
  line cap.
- Removed the hard-coded Unix `/tmp/dada2_out` default from
  `run_pipeline_samples`; the binding now resolves to
  `std::env::temp_dir()/speeddada_out` so it works on Windows.
- Removed `rustflags = ["-C", "target-cpu=native"]` from the default cargo
  config — produced binaries that crashed on older CPUs and broke
  cross-compilation. Native-CPU codegen stays available via the `bench`
  profile or `RUSTFLAGS=`.
- Deduplicated `write_temp_fastq` / `write_fastq` test helpers into a
  shared `crates/speeddada-core/src/test_util.rs` module.

### Removed

- Stale `.claude/` developer-tooling directory.
- `PLAN.md`, build artefacts that were accidentally tracked, and personal
  paths embedded in benchmark scripts.
- Dead `FastqRecord::truncate()` method.
