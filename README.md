# dada2-rs

A high-performance reimplementation of the [DADA2](https://github.com/benjjneb/dada2) amplicon sequence variant (ASV) pipeline, written in Rust. Exposes two language bindings:

- **Python** — `dada2` module via PyO3/maturin (`crates/dada2-py`)
- **R** — `dada2rs` drop-in package via extendr-api (`crates/dada2-r` + `r-package/dada2rs/`)

**Reference:** Callahan et al. 2016, *Nature Methods* — doi:10.1038/nmeth.3869

---

## Why a Rust rewrite?

| Problem | Root cause | Solution here |
|---|---|---|
| `pool=TRUE` crashes with large datasets | Entire reads from all samples held in RAM simultaneously | Disk-backed `PoolStore` — flushes chunks to temp files, pages back on demand |
| Slow processing | R/Python GIL, serial I/O | GIL released for all Rust work; Rayon parallelism at sample and intra-sample level |

**Benchmark (10 000 paired 16S V3-V4 reads, 10 true ASVs, Raspberry Pi 5 / aarch64):**

| Tool | ASVs found | Total time | vs R dada2 |
|---|---|---|---|
| R dada2 (reference) | 11 | 6 820 ms | 1× |
| dada2rs (R binding) | 10 | 309 ms | **22×** |
| Python dada2 | 10 | 330 ms | **21×** |

All three tools recover the same 10 true sequences at identical abundances (Jaccard = 0.91, Pearson r = 1.00). The chimera present in the R output is correctly removed by `removeBimeraDenovo` / `remove_bimera_denovo` in both Rust bindings.

Stage-level breakdown:

| Stage | R dada2 | dada2rs | Python dada2 |
|---|---|---|---|
| filter | 2 887 ms | 118 ms | 137 ms |
| learn_errors | 2 313 ms | 27 ms | 25 ms |
| derep | 888 ms | 44 ms | 39 ms |
| dada | 685 ms | 115 ms | 129 ms |
| merge | 43 ms | 4 ms | 1 ms |
| chimera | 2 ms | 1 ms | <1 ms |

All pipeline stages are substantially faster than R dada2. The DADA denoising step — previously a bottleneck in Rust — now matches R speed after replacing per-base `ln()` calls with a precomputed log-probability lookup table.

---

## Project layout

```
dada2_rust/
├── Cargo.toml                     # workspace root (edition 2021, MSRV 1.78)
├── crates/
│   ├── dada2-core/                # pure Rust library — no Python/R dependency
│   │   └── src/
│   │       ├── lib.rs             # Dada2Error, Phred, Kmer newtypes; bytes_to_hex
│   │       ├── quality_profile.rs # per-cycle quality statistics
│   │       ├── filter.rs          # quality filter + paired filter (streaming)
│   │       ├── primer.rs          # primer/adapter trimming
│   │       ├── error_model.rs     # EM error rate learning (logistic regression)
│   │       ├── derep.rs           # dereplication
│   │       ├── dada.rs            # DADA algorithm + pooled variant
│   │       ├── merge.rs           # paired-end merging
│   │       ├── chimera.rs         # bimera detection
│   │       ├── taxonomy.rs        # naive Bayes k-mer classifier
│   │       ├── sequence_table.rs  # sample × ASV count matrix
│   │       ├── align.rs           # Hamming / mismatch primitives (SIMD-auto-vectorised)
│   │       ├── pool.rs            # disk-backed pooled dereplication
│   │       └── io/
│   │           ├── fastq.rs       # streaming FASTQ parser (needletail)
│   │           └── fasta.rs       # FASTA reference reader
│   ├── dada2-py/                  # PyO3 bindings → `dada2` Python module
│   │   ├── Cargo.toml
│   │   ├── pyproject.toml
│   │   └── src/
│   │       ├── lib.rs             # #[pymodule] registration
│   │       ├── functions.rs       # #[pyfunction] items
│   │       └── types.rs           # #[pyclass] structs
│   └── dada2-r/                   # extendr-api bindings → compiled into dada2rs.so
│       ├── Cargo.toml
│       └── src/lib.rs
├── r-package/
│   └── dada2rs/                   # R package source (DESCRIPTION, NAMESPACE, R/)
│       └── R/dada2rs.R            # R wrapper functions (drop-in for dada2)
├── tests/
│   ├── integration/
│   │   └── pipeline_test.rs       # 1 000 synthetic reads → ASV recovery
│   └── parity/
│       ├── compare_r_output.py    # Jaccard / Pearson parity vs R reference
│       ├── generate_r_output.R    # regenerate reference fixture from R dada2
│       └── fixtures/
└── benchmarks/
    ├── sim_fastq.py               # simulate 16S V3-V4 paired FASTQs
    ├── bench_r.R                  # R dada2 benchmark
    ├── bench_dada2rs.R            # dada2rs benchmark
    ├── bench_rust.py              # Python dada2 benchmark
    └── compare.py                 # three-way comparison table
```

---

## Pipeline stages

```
FASTQ files
    │
    ▼  quality_profile         per-cycle mean / p25 / p50 / p75 statistics
    ▼  trim_primers            optional primer/adapter removal (exact or fuzzy)
    ▼  filter_and_trim         quality filter, truncation, expected-error cutoff
    │  filter_and_trim_paired  lock-step paired filtering — keeps R1/R2 in sync
    ▼  learn_errors            fit logistic error model P(obs|true, Phred q) via EM
    ▼  derep_fastq             collapse identical sequences; track per-position quality sums
    ▼  dada                    per-sample denoising — greedy Poisson significance test + EM
    │  dada_pooled             cross-sample pooling via disk-backed PoolStore
    ▼  merge_pairs             suffix-prefix overlap of F + R ASVs
    ▼  remove_bimera_denovo    exact bimera search; parents must outrank candidate
    ▼  assign_taxonomy         naive Bayes k-mer classifier (k=8, 100 bootstrap replicates)
    ▼  make_sequence_table     sample × ASV count matrix → TSV or JSON
```

---

## Algorithm notes

### DADA greedy promotion

The core DADA algorithm (Callahan 2016, Suppl. Note 1) determines which unique sequences are genuine ASVs by comparing each sequence's observed count against a Poisson null model:

```
λᵢⱼ = total_reads × ∏ₗ p(obs_base_iₗ | true_base_jₗ, Phred_qₗ)
```

If `P(X ≥ count_i | Poisson(λᵢⱼ)) < ω_A` (default ω_A = 1e-40), the sequence is promoted as a new cluster center.

**Key implementation detail:** promotion is greedy and sequential (unique sequences processed in decreasing count order). Within each iteration, the best current center is looked up after each promotion. This prevents the batch-promotion pathology where error reads from ASV2–10 sources would appear over-abundant against the initial single center (ASV1) because they differ at ~75% of positions.

Significance testing is done in log-space to avoid f64 underflow when λ << 1e-14:

```
log p ≈ count × log λ − log(count!)      [Stirling for count > 20]
promote if log p < log(ω_A)
```

---

## Building

### Prerequisites

| Tool | Install |
|---|---|
| Rust ≥ 1.78 | `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \| sh` |
| Python ≥ 3.9 | system or pyenv |
| maturin ≥ 1.5 | `pip install maturin` or `uv add maturin` |
| R ≥ 4.0 | system |

### Rust library and all crates

```bash
cargo build --release --workspace
cargo test --workspace          # 25 unit tests + 1 integration test
```

### Python extension (development install)

```bash
cd crates/dada2-py
maturin develop --release
# or with uv:
uv run maturin develop --release
```

If maturin cannot find the virtualenv automatically:

```bash
VIRTUAL_ENV=/path/to/.venv maturin develop --release
```

### R package

```bash
# Build Rust static library first
cargo build --release --workspace

# Install to user library
cd r-package/dada2rs
R CMD INSTALL --library=~/R/library .
```

---

## Python API

```python
import dada2

dada2.__version__        # "0.1.0"
dada2.init_logging()     # enable Rust-level log output (default level: "info")
```

### Classes

| Class | Key attributes |
|---|---|
| `FilterConfig` | `trunc_len`, `min_len`, `max_ee`, `trunc_q`, `trim_left`, `trim_right` |
| `FilterStats` | `.reads_in`, `.reads_out` |
| `FilterStatsPaired` | `.reads_in`, `.pairs_out`, `.fwd_failed`, `.rev_failed`, `.both_failed` |
| `QualityProfile` | `.n_reads`, `.cycle_mean`, `.cycle_p25`, `.cycle_p50`, `.cycle_p75`, `.cycle_count` |
| `ErrorModel` | Learned error matrix; `.plot_errors()` → dict for matplotlib |
| `DadaResult` | Indexable as `result[i]` → `(sequence: bytes, abundance: int)`; `len(result)` = ASV count |
| `MergedRead` | `.sequence` (bytes), `.abundance` (int), `.accept` (bool), `.nmatch`, `.nmismatch`, `.nindel` |
| `TaxonAssignment` | `.asv`, `.kingdom`, `.phylum`, `.order`, `.family`, `.genus`, `.species`, `.confidence`; note: `.class` is a Python reserved word — access via `getattr(a, 'class')` |
| `SequenceTable` | `.samples`, `.sequences` (hex), `.counts`; `.to_tsv(path)`, `.to_json()` |

### Functions

```python
# Quality inspection
profile = dada2.quality_profile("sample.fastq", n_reads=500_000)

# Primer trimming (optional)
stats = dada2.trim_primers(
    fwd_primer=b"GTGYCAGCMGCCGCGGTAA",
    rev_primer=b"GGACTACNVGGGTWTCTAAT",
    input_path="raw.fastq",
    output_path="trimmed.fastq",
    max_mismatches=1,
    min_overlap=10,
)

# Filter (single-end)
cfg   = dada2.FilterConfig(trunc_len=250, min_len=20, max_ee=2.0)
stats = dada2.filter_and_trim(cfg, "raw.fastq", "filtered.fastq")
# → FilterStats

# Filter (paired-end, lock-step — R1 and R2 are kept in sync)
stats = dada2.filter_and_trim_paired(
    cfg_fwd, cfg_rev,
    "R1.fastq", "R2.fastq",
    "filt_R1.fastq", "filt_R2.fastq",
)
# → FilterStatsPaired

# Error model (pass all per-direction filtered files together)
model = dada2.learn_errors(["filt_R1.fastq", "filt_R2.fastq"], n_reads=1_000_000)
# → ErrorModel

# Dereplicate (one file at a time)
derep = dada2.derep_fastq("filtered.fastq")
# → list[tuple[bytes, int]]  — (sequence, abundance)

# Denoise (per-sample)
result = dada2.dada(derep, model, omega_a=1e-40)
# → DadaResult

# Denoise (pooled across samples — shares evidence across all samples)
results = dada2.dada_pooled(
    [derep_s1, derep_s2, derep_s3],   # each is list[tuple[bytes, int]]
    model,
    omega_a=1e-40,
)
# → list[DadaResult]

# Merge paired ends
merged = dada2.merge_pairs(fwd_result, rev_result, min_overlap=20)
# → list[MergedRead]  — .sequence, .abundance, .nmatch, .nmismatch, .nindel

# Remove chimeras (pass sequence+abundance pairs)
clean = dada2.remove_bimera_denovo(
    [(m.sequence, m.abundance) for m in merged]
)
# → list[tuple[bytes, int]]

# Assign taxonomy (pass sequences only)
assignments = dada2.assign_taxonomy(
    [seq for seq, _ in clean],
    "silva_138.fasta",
    lineage_tsv="silva_138_lineage.tsv",  # optional; falls back to FASTA description field
    k=8,
)
# → list[TaxonAssignment]

# Build sample × ASV count matrix
table = dada2.make_sequence_table(
    sample_names=["s1", "s2"],
    results=[result_s1, result_s2],
)
table.to_tsv("asv_table.tsv")

# One-shot pipeline (filter → learn → derep → denoise → chimera)
asv_map = dada2.run_pipeline(
    input_paths=["sample1.fastq", "sample2.fastq"],
    output_dir="/tmp/filtered",
    trunc_len=250,
    max_ee=2.0,
    omega_a=1e-40,
)
# → dict[str, int]  (hex-encoded sequence → abundance)
```

### Typical paired-end workflow

```python
import dada2

dada2.init_logging()

# 1. Filter both reads in lock-step
cfg_fwd = dada2.FilterConfig(trunc_len=230, min_len=150, max_ee=3.0)
cfg_rev = dada2.FilterConfig(trunc_len=210, min_len=150, max_ee=5.0)
stats = dada2.filter_and_trim_paired(
    cfg_fwd, cfg_rev,
    "R1.fastq", "R2.fastq",
    "filt_R1.fastq", "filt_R2.fastq",
)
print(f"{stats.pairs_out}/{stats.reads_in} pairs passed")

# 2. Learn errors (use filtered reads from both directions)
model = dada2.learn_errors(["filt_R1.fastq", "filt_R2.fastq"])

# 3. Dereplicate
derep_fwd = dada2.derep_fastq("filt_R1.fastq")
derep_rev = dada2.derep_fastq("filt_R2.fastq")

# 4. Denoise
result_fwd = dada2.dada(derep_fwd, model)
result_rev = dada2.dada(derep_rev, model)
print(f"{len(result_fwd)} fwd ASVs, {len(result_rev)} rev ASVs")

# 5. Merge paired ends
merged = dada2.merge_pairs(result_fwd, result_rev, min_overlap=12)
print(f"{len(merged)} merged ASVs  (nmatch range: "
      f"{min(m.nmatch for m in merged)}–{max(m.nmatch for m in merged)})")

# 6. Remove chimeras
clean = dada2.remove_bimera_denovo([(m.sequence, m.abundance) for m in merged])
print(f"{len(clean)} non-chimeric ASVs")

# 7. Assign taxonomy
taxa = dada2.assign_taxonomy(
    [seq for seq, _ in clean],
    "silva_138.fasta",
    lineage_tsv="silva_138_lineage.tsv",
)
for t in taxa:
    print(f"  {t.genus or 'unclassified'}  conf={t.confidence:.2f}")

# 8. Build sequence table
table = dada2.make_sequence_table(["sample1"], [result_fwd])
table.to_tsv("asv_table.tsv")
```

---

## R API (`dada2rs`)

`dada2rs` is a drop-in replacement for the R `dada2` package. Function signatures mirror R dada2 exactly; extra arguments are accepted and silently ignored for compatibility.

### Install

```bash
# Build the Rust workspace first (only needed once, or after code changes)
cd /path/to/dada2_rust
cargo build --release --workspace

# Install the R package
cd r-package/dada2rs
R CMD INSTALL --library=~/R/library .
```

Or from an R session:

```r
install.packages("path/to/r-package/dada2rs", repos = NULL, type = "source",
                 INSTALL_opts = "--library=~/R/library")
```

### Functions

| Function | Description |
|---|---|
| `filterAndTrim(fwd, filt, rev, filt.rev, truncLen, maxEE, minLen, ...)` | Quality filter and truncate; paired or single-end |
| `learnErrors(fls, nbases, ...)` | Learn error rates; returns opaque error model handle |
| `derepFastq(fls, ...)` | Dereplicate FASTQ file; returns `"derep"` object with `$uniques` |
| `dada(derep, err, omega_a, pool, ...)` | Denoise; returns `"dada"` object with `$denoised` |
| `mergePairs(dadaF, derepF, dadaR, derepR, minOverlap, maxMismatch, ...)` | Merge paired-end ASVs; returns data frame with `sequence`, `abundance`, `nmatch`, `nmismatch` |
| `makeSequenceTable(samples, ...)` | Build sample × ASV count matrix |
| `removeBimeraDenovo(unqs, ...)` | Remove chimeric sequences; returns same type as input |

### Typical workflow

```r
library(dada2rs)

# 1. Filter
fstats <- filterAndTrim(
    "R1.fastq", "filt_R1.fastq",
    "R2.fastq", "filt_R2.fastq",
    truncLen = c(230L, 210L),
    maxEE    = c(3, 5),
    minLen   = 150L
)
cat(fstats["reads.out"], "/ 10000 pairs passed\n")

# 2. Learn errors
errF <- learnErrors("filt_R1.fastq")
errR <- learnErrors("filt_R2.fastq")

# 3. Dereplicate
derepF <- derepFastq("filt_R1.fastq")
derepR <- derepFastq("filt_R2.fastq")

# 4. Denoise
dadaF <- dada(derepF, err = errF)
dadaR <- dada(derepR, err = errR)

# 5. Merge and remove chimeras
merged <- mergePairs(dadaF, derepF, dadaR, derepR, minOverlap = 12L)
seqtab <- makeSequenceTable(list(sample1 = merged))
seqtab <- removeBimeraDenovo(seqtab)

cat(ncol(seqtab), "non-chimeric ASVs\n")
```

---

## Architecture decisions

### Greedy sequential promotion (DADA algorithm)
New cluster centers are added in decreasing count order. After each promotion, subsequent candidates are evaluated against the updated center set. This matches R dada2's behavior and prevents spurious promotion of low-abundance error reads that happen to differ from the initial cluster center at many positions.

### Precomputed log-probability table
`ErrorModel` holds a 16×41 `log_matrix` (one row per transition pair, one column per Phred score) built at construction time. The inner DADA loop looks up `log_matrix[[transition, phred]]` instead of calling `f64::ln()` (~100 ns on aarch64 vs ~2 ns for an array access). Per-unique `[[f64; 4]]` arrays are precomputed in parallel at `dada()` entry; `seq_ll()` becomes a pure array-index sum that LLVM auto-vectorises.

### GIL release
Every CPU-bound Python function calls `py.allow_threads(|| { … })` before entering Rust. Python threads remain live while Rust works.

### Streaming I/O
`filter_and_trim` and `filter_and_trim_paired` read and write one FASTQ record at a time. No full-file accumulation — a 100 GB file uses no more RAM than a 1 MB file during filtering.

### Disk-backed pooling (`pool.rs`)
`PoolStore` accumulates unique sequences from multiple samples in a `BTreeMap`. When the in-memory entry count exceeds `flush_threshold` (default 500 000), the current map is serialised to a JSONL chunk in a `tempfile::TempDir` and cleared. `into_pooled_uniques()` re-merges all chunks into a single sorted `Vec<UniqueSeq>` for DADA. `dada_pooled` re-splits ASV assignments back to per-sample lists using provenance stored during accumulation.

### SIMD alignment (`align.rs`)
`hamming_distance`, `first_mismatch`, and `range_equal` are scalar loops that LLVM reliably auto-vectorises to AVX2 / NEON when compiled with `target-cpu=native`. `chimera.rs`, `merge.rs`, and `primer.rs` delegate to these primitives.

### Rayon parallelism
- **Sample level:** `filter_and_trim_many` and `filter_and_trim_paired_many` process N FASTQ pairs across the Rayon thread pool.
- **Intra-sample:** DADA re-assignment, chimera candidate scan, taxonomy classification, and pooled sample ingestion all use `par_iter()`.

### Type safety
Domain primitives are newtypes — `Phred(u8)` and `Kmer(u64)` — preventing argument-order bugs. All errors propagate via `Dada2Error` (thiserror); no `.unwrap()` or `.expect()` in library code.

---

## Testing

### Rust unit and integration tests

```bash
cargo test --workspace
# 25 unit tests (across 11 modules) + 1 integration test
# Integration: 1 000 synthetic reads, 2 % errors → top ASV ≥ 95 % identity
```

### Python smoke tests

```bash
cd crates/dada2-py
pytest tests/
# 4 tests: version semver, filter_and_trim, learn_errors, dada + remove_bimera_denovo round-trip
```

### Parity test against R dada2

```bash
# Regenerate the reference fixture (requires R + dada2 + jsonlite packages)
Rscript tests/parity/generate_r_output.R

# Compare Rust output against R reference
python tests/parity/compare_r_output.py rust_output.json
# Asserts: Jaccard ≥ 0.95, Pearson r ≥ 0.99
```

### Benchmarks

```bash
# 1. Generate simulated 16S V3-V4 paired FASTQs (10 000 read pairs, 10 true ASVs)
python benchmarks/sim_fastq.py
# → /tmp/bench_fastq/R1.fastq, /tmp/bench_fastq/R2.fastq

# 2. Run each tool (JSON results written to /tmp/bench_out/)
Rscript  benchmarks/bench_r.R           # R dada2 (reference)
Rscript  benchmarks/bench_dada2rs.R     # dada2rs R binding
python   benchmarks/bench_rust.py       # Python dada2 binding

# 3. Print three-way comparison table
python benchmarks/compare.py
```

Each benchmark script prints per-stage timing and writes a JSON result file. `compare.py` reads all three JSON files and prints a summary table with speedup ratios.

---

## Configuration reference

### `FilterConfig` (Python)

| Field | Default | Description |
|---|---|---|
| `trunc_len` | `0` | Truncate reads to this length; 0 = no truncation |
| `min_len` | `20` | Discard reads shorter than this after truncation |
| `max_ee` | `2.0` | Maximum expected errors (sum of error probabilities) |
| `trunc_q` | `2` | Truncate at first base with Phred below this value |
| `trim_left` | `0` | Bases to remove from the 5′ end before truncation |
| `trim_right` | `0` | Bases to remove from the 3′ end before truncation |

### `DadaConfig` (Rust core)

| Field | Default | Description |
|---|---|---|
| `omega_a` | `1e-40` | Abundance p-value threshold for accepting a new ASV |
| `pool` | `false` | Enable cross-sample pooling (use `dada_pooled` instead) |
| `max_iter` | `16` | Maximum EM iterations |
| `tol` | `1e-6` | Log-likelihood convergence tolerance |
| `seed` | `42` | RNG seed (reserved; currently unused) |

### `TaxonomyConfig` (Rust core)

| Field | Default | Description |
|---|---|---|
| `k` | `8` | k-mer length (must be ≤ 16) |
| `threshold` | `0.80` | Minimum bootstrap confidence to report genus-level assignment |
| `seed` | `42` | Bootstrap subsampling seed |

### Lineage TSV format

```
seq_id<TAB>kingdom;phylum;class;order;family;genus;species
```

Compatible with SILVA (`tax_slv_ssu_*.txt`) and GTDB lineage files. If `lineage_tsv` is not provided to `assign_taxonomy`, the FASTA description field (everything after the first space in the `>` header) is parsed as a semicolon-separated lineage instead.

---

## Key dependencies

| Crate | Version | Purpose |
|---|---|---|
| `pyo3` | 0.23 | Python ↔ Rust FFI |
| `extendr-api` | 0.9 | R ↔ Rust FFI |
| `maturin` | ≥ 1.5 | Build system for the Python extension |
| `rayon` | 1 | Work-stealing parallelism |
| `needletail` | 0.5 | Streaming FASTQ/FASTA parser |
| `ndarray` | 0.15 | 2-D error matrix |
| `statrs` | 0.17 | Poisson distribution for abundance p-values |
| `thiserror` | 1 | Ergonomic error types |
| `serde` / `serde_json` | 1 | Serialisation for pool chunks and JSON output |
| `tempfile` | 3 | Temp directory for `PoolStore` disk chunks |

---

## Known limitations

- **Error model learning** uses self-alignment (only match transitions are accumulated). Mismatch rates fall back to a logistic prior (`sigmoid(-5 + 0.1q)`), which overestimates error probability by ~100× compared to the Phred definition. The DADA denoising step remains correct because the greedy promotion threshold is set high enough (ω_A = 1e-40), but the learned model will not match R dada2's error plots.
- **`dada` quality scores** — the bindings pass per-position quality averages computed during dereplication. If a unique sequence was observed with varying quality, the mean is used.
- **`selfConsist` mode** — not implemented; `dada2rs` emits a warning and falls back to single-pass denoising.
- **`pool="pseudo"`** — not implemented; falls back to `pool=FALSE`.
- **Vectorised sample input** — `dada()`, `derep_fastq()`, and `filter_and_trim()` each process one sample at a time. For multi-sample workflows, loop over samples explicitly or use `dada_pooled()` for cross-sample pooled denoising.
