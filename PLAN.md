# dada2-rs — Improvement Plan

Priority ranking: **P1** = blocks common real-world workflows, **P2** = significant gap vs. R dada2,
**P3** = quality-of-life / correctness hardening.

---

## P1 — Wire pool=TRUE end-to-end

**What is broken today:**
`DadaConfig.pool: bool` is stored but never read inside `run_dada()`. The `PoolStore`
implementation in `pool.rs` exists and works in isolation, but nothing calls it.
Passing `pool=True` from Python currently runs as single-sample mode silently.

**Plan:**

1. Add `run_dada_pooled(samples: &[&[UniqueSeq]], em: &ErrorModel, cfg: &DadaConfig)`
   to `dada.rs`. Internally:
   - Construct a `PoolStore::new(500_000)`.
   - Call `store.add_sample(i, samples[i])` for each sample with Rayon.
   - Call `store.into_pooled_uniques()` → get the merged `Vec<UniqueSeq>`.
   - Run the existing `run_dada()` on the pooled uniques.
   - Re-split the resulting ASV assignments back to per-sample `Vec<Asv>` using the
     `PoolEntry.per_sample` provenance stored during accumulation.
   - Return `Vec<Vec<Asv>>` — one per input sample.

2. Expose in `dada2-py`:
   ```python
   run_dada_pooled(
       samples: list[list[tuple[bytes, int]]],   # one derep per sample
       error_model: ErrorModel,
       omega_a: float = 1e-40,
   ) -> list[DadaResult]
   ```

3. Add a Rust integration test: 3 synthetic samples, 500 reads each, pool=True →
   verify the pooled run recovers the same ASVs as per-sample runs on identical data.

4. Add a Python test in `dada2_test/tests/test_dada.py` for `run_dada_pooled`.

**Files:** `dada.rs`, `pool.rs`, `dada2-py/src/lib.rs`, tests.

---

## P1 — Paired-end simultaneous filtering

**What is broken today:**
`filter_and_trim` processes a single file. Users must filter R1 and R2 independently,
which can produce unmatched pairs (a read passing in R1 but failing in R2). R dada2's
`filterAndTrim` keeps pairs in lock-step.

**Plan:**

1. Add `filter_and_trim_paired(cfg_fwd, cfg_rev, r1_in, r2_in, r1_out, r2_out)`
   to `filter.rs`. Algorithm:
   - Iterate both FASTQ files in lock-step using two needletail cursors.
   - Apply `apply_filters_owned` to both reads.
   - Write the pair to output only if **both** pass; otherwise discard both.
   - Return `FilterStatsPaired { reads_in, pairs_out, fwd_failed, rev_failed, both_failed }`.

2. Expose in `dada2-py`:
   ```python
   filter_and_trim_paired(
       config_fwd: FilterConfig,
       config_rev: FilterConfig,
       r1_in: str, r2_in: str,
       r1_out: str, r2_out: str,
   ) -> FilterStatsPaired
   ```

3. Add `filter_and_trim_paired_many` for Rayon sample-level parallelism.

4. Add Python tests for paired filtering: assert `reads_out_fwd == reads_out_rev`
   and that discarding one direction removes both.

**Files:** `filter.rs`, `dada2-py/src/lib.rs`, tests.

---

## P1 — Taxonomy lineage TSV loader

**What is broken today:**
`assign_taxonomy_py` parses lineage by splitting the FASTA description field on `;`.
Standard reference databases (SILVA, GTDB, RDP) ship lineage in a separate TSV/CSV
file, not embedded in headers.

**Plan:**

1. Add `load_lineage_tsv(path: &Path) -> Result<HashMap<String, Vec<String>>, Dada2Error>`
   to `taxonomy.rs`. Expected TSV format (tab-separated):
   ```
   seq_id\tkingdom;phylum;class;order;family;genus;species
   ```
   This matches the SILVA `tax_slv_ssu_*.txt` format.

2. Update `assign_taxonomy_py` to accept an optional `lineage_tsv: Option<str>`
   parameter. When provided, load from TSV; when `None`, fall back to parsing FASTA
   headers (current behaviour).

3. Add a fixture TSV and a unit test that round-trips loading → classifying.

**Files:** `taxonomy.rs`, `dada2-py/src/lib.rs`, test fixtures.

---

## P2 — Quality inspection (Stage 1)

**What is missing:**
R dada2's `plotQualityProfile` is often the first step — it shows per-cycle mean
quality to inform `trunc_len` choice. There is currently no equivalent.

**Plan:**

1. Add `quality_profile.rs` to `dada2-core/src/`. Expose:
   ```rust
   pub struct QualityProfile {
       pub n_reads: u64,
       pub cycle_mean: Vec<f64>,   // mean Phred per cycle
       pub cycle_p25:  Vec<f64>,   // 25th percentile
       pub cycle_p50:  Vec<f64>,   // median
       pub cycle_p75:  Vec<f64>,   // 75th percentile
       pub cycle_count: Vec<u64>,  // reads reaching this cycle
   }

   pub fn quality_profile(path: &Path, n_reads: usize) -> Result<QualityProfile, Dada2Error>
   ```
   Stream through the file, accumulate per-cycle quality histograms (u32 counts per
   Phred value), then compute percentiles in a single pass. Memory use is
   `O(max_read_len × 41)` regardless of file size.

2. Expose in Python:
   ```python
   profile = dada2.quality_profile("sample.fastq", n_reads=100_000)
   # profile.cycle_mean, .cycle_p25, .cycle_median, .cycle_p75, .cycle_count
   # → pass to matplotlib directly
   ```

3. Add unit test: 50 synthetic reads → assert `len(cycle_mean) == read_length`.

**Files:** `quality_profile.rs` (new), `lib.rs`, `dada2-py/src/lib.rs`.

---

## P2 — Sequence table output

**What is missing:**
R dada2's `makeSequenceTable` produces a sample × ASV count matrix — the canonical
deliverable of the pipeline. There is no equivalent; `run_pipeline` only returns a
pooled `dict[str, int]` with no per-sample breakdown.

**Plan:**

1. Define a `SequenceTable` type in `dada2-core`:
   ```rust
   pub struct SequenceTable {
       pub samples: Vec<String>,
       pub sequences: Vec<Vec<u8>>,
       pub counts: Array2<u32>,    // shape: [n_samples, n_asvs]
   }
   impl SequenceTable {
       pub fn to_json(&self) -> Result<String, Dada2Error>
       pub fn to_tsv(&self, path: &Path) -> Result<(), Dada2Error>
   }
   ```

2. Add `make_sequence_table(sample_names: &[&str], results: &[Vec<Asv>]) -> SequenceTable`.

3. Expose in Python:
   ```python
   table = dada2.make_sequence_table(
       sample_names=["s1", "s2"],
       results=[dada_s1, dada_s2],
   )
   table.to_tsv("asv_table.tsv")
   # or: table.counts → numpy array, table.sequences → list[bytes]
   ```
   Return a dict `{"samples": [...], "sequences": [...], "counts": [[...]]}` that
   converts trivially to a pandas/polars DataFrame.

4. Update `run_pipeline` to return a `SequenceTable` rather than a flat dict when
   multiple input files are given.

**Files:** `sequence_table.rs` (new), `lib.rs`, `dada2-py/src/lib.rs`.

---

## P2 — Primer trimming

**What is missing:**
Before filtering, amplicon reads typically contain primer sequences that must be
removed. R users call `cutadapt` externally. A built-in exact/near-exact primer
trimmer would make the pipeline self-contained.

**Plan:**

1. Add `primer.rs` to `dada2-core/src/`. Expose:
   ```rust
   pub struct PrimerConfig {
       pub fwd_primer: Vec<u8>,
       pub rev_primer: Vec<u8>,
       pub max_mismatches: u32,   // 0 = exact match only
       pub min_overlap: usize,    // minimum primer overlap to trim
   }

   pub fn trim_primers(cfg: &PrimerConfig, input: &Path, output: &Path)
       -> Result<FilterStats, Dada2Error>
   ```
   Algorithm: for each read, search for the primer in the first
   `primer_len + max_mismatches + 5` bases using `align::hamming_distance`;
   if found within tolerance, trim it; otherwise discard the read.

2. Expose in Python:
   ```python
   dada2.trim_primers(fwd_primer="GTGYCAGCMGCCGCGGTAA",
                      rev_primer="GGACTACNVGGGTWTCTAAT",
                      input_path="raw.fastq", output_path="trimmed.fastq",
                      max_mismatches=1)
   ```

3. Integrate into `run_pipeline` as an optional pre-filter step.

**Files:** `primer.rs` (new), `lib.rs`, `dada2-py/src/lib.rs`.

---

## P3 — Error model graceful fallback

**What is broken today:**
`learn_errors` returns `Dada2Error::Convergence` when the gradient descent fails
(demonstrated with `n_reads=1` in the test suite). Callers must handle this themselves.
In a pipeline context this crashes the whole run.

**Plan:**

1. In `error_model.rs`, change `fit_logistic_row` to return the best-so-far
   parameters rather than `Err(Convergence)` when `max_iter` is exhausted — log a
   `warn!` but do not error.

2. Only return `Err(Convergence)` if the very first EM step produces `NaN` or `±Inf`
   (true numerical failure), not just slow convergence.

3. Remove the `try/except RuntimeError` workaround from `test_learn_errors_n_reads_1`
   in the Python test suite and assert the call always succeeds.

**Files:** `error_model.rs`, `dada2_test/tests/test_error_model.py`.

---

## P3 — Progress reporting

**What is missing:**
`learn_errors` on 1 M reads and pooled DADA across many samples can take minutes with
no feedback. Long-running operations should emit structured log messages.

**Plan:**

1. Add `log::info!` calls at the start and end of each major stage, and after each
   EM iteration in `learn_errors` and `run_dada` (gated behind `log::debug!`).
   Example: `log::info!("learn_errors: processed {n} reads, iter {i}/{max}, ΔlogL={delta:.2e}")`.

2. In `dada2-py`, initialise `env_logger` from Python:
   ```python
   dada2.init_logging(level="INFO")   # sets RUST_LOG and calls env_logger::init
   ```
   This lets Python users opt in to Rust-level log output without touching env vars.

3. No external dependency needed — `log` and `env_logger` are already in the
   dependency graph.

**Files:** all stage files (add log calls), `dada2-py/src/lib.rs` (add `init_logging`).

---

## P3 — Parity test automation

**What is missing:**
`tests/parity/fixtures/r_output.json` is a hand-written static fixture, not generated
from real R dada2. The parity test cannot detect regressions against the reference
implementation.

**Plan:**

1. Add `tests/parity/generate_r_output.R` — a self-contained R script that:
   - Generates the same 1 000-read synthetic FASTQ (same seed/parameters as the Rust
     integration test).
   - Runs the R dada2 pipeline on it.
   - Writes the ASV table to `fixtures/r_output.json` in the expected format.

2. Document in README how to regenerate: `Rscript tests/parity/generate_r_output.R`.

3. Commit the regenerated `r_output.json` so CI can run `compare_r_output.py` against
   a real R baseline rather than the toy two-entry fixture.

**Files:** `tests/parity/generate_r_output.R` (new), `tests/parity/fixtures/r_output.json`.

---

## Dependency summary

```
P1 pool=TRUE          ← depends on: pool.rs (done)
P1 paired filtering   ← no dependencies
P1 lineage TSV        ← no dependencies
P2 quality profile    ← no dependencies
P2 sequence table     ← depends on: P1 pool=TRUE (for run_pipeline multi-sample)
P2 primer trimming    ← depends on: align.rs (done)
P3 error fallback     ← no dependencies
P3 progress logging   ← no dependencies
P3 parity automation  ← no dependencies
```

## Suggested implementation order

| Step | Item | Reason |
|------|------|--------|
| 1 | P3 error fallback | Tiny change; removes a known crash path; unblocks other tests |
| 2 | P1 paired filtering | Self-contained; unlocks the most common real-world use case |
| 3 | P1 lineage TSV | Self-contained; makes taxonomy usable with real reference DBs |
| 4 | P2 quality profile | No dependencies; naturally precedes filtering in the workflow |
| 5 | P1 pool=TRUE | Builds on existing pool.rs; unlocks large-dataset runs |
| 6 | P2 sequence table | Depends on pooled results being correct |
| 7 | P2 primer trimming | Builds on align.rs; completes the pre-filter stage |
| 8 | P3 progress logging | Polish; safe to do at any point |
| 9 | P3 parity automation | Final validation step; best done last |
