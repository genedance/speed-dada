# SpeedDada NEWS

## SpeedDada 0.99.1

* Prebuilt binary R packages now attached to each GitHub Release for
  macOS Apple Silicon (`.tgz`), Windows x64 (`.zip`), and Linux x86_64
  (`_R_x86_64-pc-linux-gnu.tar.gz`). Installing them does not invoke
  the Rust compiler on the user's machine — non-technical users no
  longer need rustup on the supported platforms.
* No functional changes to the pipeline.

## SpeedDada 0.99.0

* Initial public release prepared for Bioconductor submission.
* Drop-in compatible wrappers for the seven dada2 pipeline functions:
  `filterAndTrim`, `learnErrors`, `derepFastq`, `dada`, `mergePairs`,
  `makeSequenceTable`, and `removeBimeraDenovo`.
* Rust core (extendr-api bindings) typically 10-20x faster than R
  `dada2` and ~10x lower peak memory.
* Cross-sample pooling via `dada(..., pool = TRUE)` and pseudo-pooling
  via `dada(..., pool = "pseudo")`.
* Vignette: end-to-end paired-end pipeline using the bundled FASTQ
  fixtures (`vignette("SpeedDada-pipeline", "SpeedDada")`).
* Builds out-of-the-box with a Rust toolchain on Linux (x86_64 / aarch64),
  macOS (x86_64 / arm64), and Windows (x86_64), including Raspberry Pi 5.
