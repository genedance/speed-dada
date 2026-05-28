# SpeedDada NEWS

## SpeedDada 0.99.2 (in development)

### dada2 R-API parity

* New R-facing functions matching the dada2 v1.36 surface:
  `assignTaxonomy`, `assignSpecies`, `addSpecies`, `removePrimers`,
  `plotQualityProfile`, `rc`, `getSequences`, `getUniques`,
  `uniquesToFasta`, `mergeSequenceTables`, and the database helpers
  `buildTaxonomyDb` + `assignTaxonomyDb`.
* `learnErrors` now collects real mismatch evidence (selfConsist-lite,
  same-length pairwise) instead of falling back to a logistic prior,
  and smooths the resulting transition rates with weighted local linear
  regression (LOWESS, degree 1) — equivalent to `dada2::loessErrfun`.

### Cross-platform error-model support

* `learnErrors(errFun = "auto")` (the new default) sniffs the input
  quality profile and picks the right smoother:
  - `loess` — full-range Illumina (MiSeq, HiSeq full quality).
  - `binned` — NovaSeq (4-bin), NextSeq (8-bin), MGI DNBSEQ.
    Piecewise-linear interpolation between observed bins, mirroring
    `dada2::makeBinnedQualErrfun`.
  - `pacbio` — PacBio CCS / HiFi (near-Q40). Canned empirical matrix.
* Oxford Nanopore data is detected and warned about. The pipeline still
  runs but ASVs will be inaccurate — proper ONT support needs banded
  indel-aware alignment (separate body of work).

### Supported platforms

| Platform                         | Status         | Notes                                  |
|----------------------------------|----------------|----------------------------------------|
| Illumina MiSeq / HiSeq           | Supported      | LOESS error model                      |
| Illumina NovaSeq 6000 / X        | Supported      | Auto-binned (4-bin) error model        |
| Illumina NextSeq 1000/2000       | Supported      | Auto-binned (8-bin) error model        |
| MGI / Complete Genomics DNBSEQ   | Supported      | Auto-binned error model                |
| PacBio CCS / HiFi                | Supported      | Canned PacBio error matrix             |
| Oxford Nanopore                  | Not supported  | Warning emitted; output unreliable     |

### Testing

* New `SPEEDDADA_PARITY=1`-gated parity test suite compares results
  cell-by-cell against the reference Bioconductor dada2 package across
  every supported platform fixture.

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
