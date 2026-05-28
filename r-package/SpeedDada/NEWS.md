# SpeedDada NEWS

## SpeedDada 0.99.2

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
  every supported platform fixture. dada2 reference outputs are baked
  into `tests/testthat/fixtures-dada2/` (regenerate with
  `inst/scripts/bake-dada2-snapshots.R` on a machine with dada2
  installed) so the suite runs even when dada2 isn't available.

### Filter bug fixes uncovered by parity tests

* `filterAndTrim`: quality-score truncation now uses `<= trunc_q`
  (matching dada2), not `< trunc_q`. Important on platforms with a Q2
  bin (NovaSeq, NextSeq, MGI): without this fix every read kept its
  Q2 bases and the downstream error model diverged from dada2.
* `filterAndTrim`: reads shorter than `truncLen` are now discarded
  (matching dada2 semantics). Previously they fell through to the
  `minLen` check, which let in many error-prone short reads that
  dada2 had filtered out.

### Cross-platform parity vs dada2 1.38.0

After the filter and error-model fixes in this release, ASV-set
agreement against the reference Bioconductor dada2 1.38.0 on the
synthetic platform fixtures is:

| Platform           | ASV-set Jaccard |
|--------------------|----------------:|
| PacBio CCS         | 1.00 (exact)    |
| Illumina MiSeq     | 1.00 (exact)    |
| Illumina HiSeq     | 1.00 (exact)    |
| Illumina NovaSeq   | 1.00 (exact)    |
| Illumina NextSeq   | 1.00 (exact)    |
| MGI DNBSEQ         | 0.50            |

The MGI case is the one outlier: SpeedDada returns 4 ASVs, dada2
returns 2 (out of 5 underlying truth ASVs in the fixture). dada2
collapses more aggressively on the MGI quality profile; SpeedDada is
actually closer to truth there. We hold the bar at exact equality on
the other five platforms.

### Error-model improvements that closed the binned-quality gap

* **Pooled mismatch evidence**: per-quality mismatch rates are now
  estimated by pooling across all 12 non-self transitions instead of
  fitting each direction independently. The uniform-substitution
  assumption costs a little per-direction bias but makes the estimates
  robust on small samples where individual cells previously collapsed
  to the rate floor and made `dada()` over-split.
* **All-pairs mismatch evidence collection**: the selfConsist-lite
  pairwise pass now compares every pair of dereplicated reads within a
  length bucket (~n²/2), not just bucket[0]-vs-others. On the noisy
  binned-quality fixtures every read is a unique singleton; the
  parent-vs-others heuristic only found ~n pairs which left high-Q
  cells empty. With ~n²/2 pairs the smoother gets enough evidence at
  every (transition, q) bin to populate the matrix.
* **dada2-compatible rate floor/cap**: mismatch rates are clamped to
  `[MIN_ERROR_RATE, MAX_ERROR_RATE] = [1e-7, 0.25]` matching
  `dada2::makeBinnedQualErrfun`.

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
