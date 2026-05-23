# dada2-rs

A high-performance reimplementation of the [DADA2](https://github.com/benjjneb/dada2) amplicon sequence variant (ASV) pipeline, written in Rust with Python bindings via PyO3/maturin.

**Reference:** Callahan et al. 2016, *Nature Methods* — doi:10.1038/nmeth.3869

---

## Why a Rust rewrite?

The original R package has two fundamental bottlenecks that this project addresses at the architecture level:

| Problem | Root cause | Solution here |
|---|---|---|
| `pool=TRUE` crashes with large datasets | Entire reads from all samples held in RAM simultaneously | Disk-backed `PoolStore` — flushes chunks to temp files, pages back on demand; RAM stays flat |
| Slow processing | R/Python GIL, no SIMD, serial I/O | GIL released for all Rust work; LLVM auto-vectorised alignment; Rayon parallelism at both sample and intra-sample level |

---

## Project layout

```
dada2_rust/
├── Cargo.toml                     # workspace root (edition 2021, MSRV 1.78)
├── .cargo/config.toml             # target-cpu=native for SIMD; release LTO settings
├── crates/
│   ├── dada2-core/                # pure Rust library — no Python dependency
│   │   └── src/
│   │       ├── lib.rs             # Dada2Error, Phred, Kmer newtypes
│   │       ├── quality_profile.rs # stage 1 — per-cycle quality statistics
│   │       ├── filter.rs          # stage 2 — streaming quality filter + paired filter
│   │       ├── primer.rs          # stage 2a — primer/adapter trimming
│   │       ├── error_model.rs     # stage 3 — EM error rate learning
│   │       ├── derep.rs           # stage 4 — dereplication
│   │       ├── dada.rs            # stage 5 — core DADA algorithm + pooled variant
│   │       ├── merge.rs           # stage 6 — paired-end merging
│   │       ├── chimera.rs         # stage 7 — bimera detection
│   │       ├── taxonomy.rs        # stage 8 — naive Bayes k-mer classifier + TSV loader
│   │       ├── sequence_table.rs  # sample × ASV count matrix
│   │       ├── align.rs           # SIMD-ready Hamming / mismatch primitives
│   │       ├── pool.rs            # disk-backed pooled dereplication
│   │       └── io/
│   │           ├── fastq.rs       # streaming FASTQ parser (needletail)
│   │           └── fasta.rs       # FASTA reference reader
│   └── dada2-py/                  # PyO3 bindings — compiled to a .so via maturin
│       ├── Cargo.toml
│       ├── pyproject.toml
│       └── src/lib.rs
├── tests/
│   ├── integration/
│   │   └── pipeline_test.rs       # 1 000 synthetic reads → ASV recovery test
│   └── parity/
│       ├── compare_r_output.py    # Jaccard / Pearson parity vs saved R output
│       ├── generate_r_output.R    # regenerate reference fixture from R dada2
│       └── fixtures/r_output.json
└── scripts/
    └── check.sh                   # full CI gate (fmt → clippy → test → build → pytest)
```

---

## Pipeline stages

```
FASTQ files
    │
    ▼  quality_profile       per-cycle mean/p25/p50/p75 — informs trunc_len choice
    │                        streaming histogram; O(max_read_len × 41) RAM
    ▼  trim_primers          optional primer/adapter removal before quality filtering
    │                        exact or fuzzy match (max_mismatches); uses SIMD Hamming
    ▼  filter_and_trim       quality filter, length truncation, expected-error cutoff
    │  filter_and_trim_paired  lock-step paired filtering — keeps R1/R2 in sync
    │                        streaming: one record at a time — O(1) RAM
    ▼  learn_errors          fit logistic error model P(obs|true, Phred q) via EM
    │                        16-class (all base transitions) × 41 Phred bins
    │                        graceful fallback to Illumina default on slow convergence
    ▼  dereplicate           collapse identical sequences; track per-position quality sums
    ▼  run_dada              per-sample denoising — Poisson abundance p-values, EM
    │  run_dada_pooled       cross-sample pooling via disk-backed PoolStore
    │                        OMEGA_A default 1e-40; convergence tol 1e-6; max 16 iters
    ▼  merge_pairs           suffix-prefix overlap of F+R reads (min_overlap=20 default)
    ▼  remove_bimeras        exact bimera search; min arm 8 bp; parent must outrank candidate
    ▼  assign_taxonomy       naive Bayes k-mer classifier (k=8, 100 bootstrap replicates)
    │                        lineage from FASTA headers or external TSV
    ▼  make_sequence_table   sample × ASV count matrix → TSV or JSON
```

---

## Architecture decisions

### GIL release
Every CPU-bound Python function calls `py.allow_threads(|| { … })` before entering Rust. Python threads remain live while Rust works. Objects borrowed from Python (`ErrorModel`, `DadaResult`) are cloned before crossing the thread boundary.

### Streaming I/O
`filter_and_trim` and `filter_and_trim_paired` read and write one FASTQ record at a time via a needletail cursor + `BufWriter`. No full-file accumulation — a 100 GB file uses no more RAM than a 1 MB file during filtering.

### Disk-backed pooling (`pool.rs`)
`PoolStore` accumulates unique sequences from multiple samples in a `BTreeMap`. When the in-memory entry count exceeds `flush_threshold` (default 500 000), the current map is serialised to a JSONL chunk in a `tempfile::TempDir` and cleared. `into_pooled_uniques()` re-merges all chunks into a single sorted `Vec<UniqueSeq>` for DADA. `run_dada_pooled` re-splits ASV assignments back to per-sample lists using the `PoolEntry.per_sample` provenance stored during accumulation. RAM use is proportional to unique sequence count, not raw read count.

### SIMD alignment (`align.rs`)
`hamming_distance`, `first_mismatch`, and `range_equal` are tight scalar loops that LLVM reliably auto-vectorises to AVX2 / NEON / SSE4 when compiled with `target-cpu=native` (set in `.cargo/config.toml`). No `unsafe` blocks, no third-party SIMD crate. `chimera.rs`, `merge.rs`, and `primer.rs` all delegate to these primitives.

### Rayon parallelism
- **Sample level:** `filter_and_trim_many` and `filter_and_trim_paired_many` process N FASTQ pairs across the Rayon thread pool.
- **Intra-sample:** the DADA E-step, chimera candidate scan, taxonomy classification, and pooled sample ingestion all use `par_iter()`.

### Error model robustness
`learn_errors` returns a valid `ErrorModel` even when the gradient descent has not fully converged — it logs a warning and returns the best-so-far parameters. Only a numerical failure (NaN/Inf) raises `Dada2Error::Convergence`. This prevents a single low-read-count sample from aborting a multi-sample run.

### Type safety
Domain primitives are newtypes — `Phred(u8)` and `Kmer(u64)` — preventing argument-order bugs. All errors propagate via `Dada2Error` (thiserror); no `.unwrap()` or `.expect()` in library code.

---

## Building

### Prerequisites

| Tool | Install |
|---|---|
| Rust ≥ 1.78 | `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \| sh` |
| Python ≥ 3.9 | system or pyenv |
| uv | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| maturin ≥ 1.5 | installed automatically via uv below |

### Rust library only

```bash
cd dada2_rust
cargo build --release --workspace   # optimised build
cargo test --workspace              # 25 tests (24 unit + 1 integration)
```

### Python extension (development install)

```bash
cd dada2_rust
PATH="$HOME/.cargo/bin:$PATH" maturin develop --manifest-path crates/dada2-py/Cargo.toml

# Or from the dada2_test project using uv:
cd ~/dada2_test
uv sync
PATH="$HOME/.cargo/bin:$PATH" uv run maturin develop \
    --manifest-path ~/dada2_rust/crates/dada2-py/Cargo.toml
```

### Release wheel

```bash
cd dada2_rust/crates/dada2-py
maturin build --release   # produces target/wheels/dada2-0.1.0-*.whl
```

---

## Python API

```python
import dada2

dada2.__version__        # "0.1.0"
dada2.init_logging()     # enable Rust-level log output (default level: "info")
```

### Classes

| Class | Description |
|---|---|
| `FilterConfig` | Parameters for filter_and_trim (trunc_len, min_len, max_ee, trunc_q, trim_left, trim_right) |
| `FilterStats` | `.reads_in`, `.reads_out` returned by `filter_and_trim` |
| `FilterStatsPaired` | `.reads_in`, `.pairs_out`, `.fwd_failed`, `.rev_failed`, `.both_failed` |
| `ErrorModel` | Learned error matrix; `.plot_errors()` → dict for matplotlib |
| `DadaResult` | Sequence of `(sequence: bytes, abundance: int)` tuples |
| `TaxonAssignment` | Per-ASV classification with `.kingdom` … `.species` and `.confidence` |
| `QualityProfile` | Per-cycle `.cycle_mean`, `.cycle_p25`, `.cycle_p50`, `.cycle_p75`, `.cycle_count` |
| `SequenceTable` | Sample × ASV matrix with `.to_tsv(path)` and `.to_json()` |

### Functions

```python
# Stage 1 — quality inspection
profile = dada2.quality_profile("sample.fastq", n_reads=500_000)

# Stage 2a — primer trimming (optional)
stats = dada2.trim_primers(
    fwd_primer="GTGYCAGCMGCCGCGGTAA",
    rev_primer="GGACTACNVGGGTWTCTAAT",
    input_path="raw.fastq",
    output_path="trimmed.fastq",
    max_mismatches=1,
)

# Stage 2 — filter (single-end)
stats = dada2.filter_and_trim(config, input_path, output_path)

# Stage 2 — filter (paired-end, lock-step)
stats = dada2.filter_and_trim_paired(
    cfg_fwd, cfg_rev, r1_in, r2_in, r1_out, r2_out
)

# Stage 3 — error model
model = dada2.learn_errors(fastq_paths, n_reads=1_000_000)

# Stage 4 — dereplicate
derep = dada2.dereplicate(fastq_path)           # → list[(bytes, int)]

# Stage 5 — denoise (per-sample)
result = dada2.run_dada(derep, model, omega_a=1e-40)

# Stage 5 — denoise (pooled across samples)
results = dada2.run_dada_pooled(
    [derep_s1, derep_s2, derep_s3],
    model,
    omega_a=1e-40,
)  # → list[DadaResult], one per sample

# Stage 6 — merge paired ends
merged = dada2.merge_pairs(fwd_result, rev_result, min_overlap=20)

# Stage 7 — remove chimeras
clean = dada2.remove_bimeras(seqs)              # seqs: list[(bytes, int)]

# Stage 8 — taxonomy (lineage from FASTA headers)
assignments = dada2.assign_taxonomy(seqs, ref_fasta, k=8)

# Stage 8 — taxonomy (lineage from separate TSV file)
assignments = dada2.assign_taxonomy(
    seqs, ref_fasta, lineage_tsv="silva_lineage.tsv"
)

# Build sample × ASV count table
table = dada2.make_sequence_table(
    sample_names=["s1", "s2", "s3"],
    results=[result_s1, result_s2, result_s3],
)
table.to_tsv("asv_table.tsv")

# Full pipeline in one call (GIL-free, single-end)
asv_table = dada2.run_pipeline(
    input_paths=["sample1.fastq", "sample2.fastq"],
    output_dir="/tmp/filtered",
    trunc_len=250,
    max_ee=2.0,
    omega_a=1e-40,
)  # → dict[str, int]  (hex-encoded sequence → abundance)
```

### Typical paired-end workflow

```python
import dada2

dada2.init_logging()   # show progress from Rust

# Inspect quality to choose trunc_len
fwd_profile = dada2.quality_profile("sample_R1.fastq")
rev_profile = dada2.quality_profile("sample_R2.fastq")
# plot fwd_profile.cycle_mean, rev_profile.cycle_mean with matplotlib...

# Optional: strip primers
dada2.trim_primers("GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT",
                   "sample_R1.fastq", "trimmed_R1.fastq")
dada2.trim_primers("GGACTACNVGGGTWTCTAAT", "GTGYCAGCMGCCGCGGTAA",
                   "sample_R2.fastq", "trimmed_R2.fastq")

# Filter both reads in lock-step
cfg_fwd = dada2.FilterConfig(trunc_len=240, max_ee=2.0)
cfg_rev = dada2.FilterConfig(trunc_len=200, max_ee=2.0)
stats = dada2.filter_and_trim_paired(
    cfg_fwd, cfg_rev,
    "trimmed_R1.fastq", "trimmed_R2.fastq",
    "filtered_R1.fastq", "filtered_R2.fastq",
)
print(f"{stats.pairs_out}/{stats.reads_in} pairs passed")

# Learn errors from filtered reads
model_fwd = dada2.learn_errors(["filtered_R1.fastq"])
model_rev = dada2.learn_errors(["filtered_R2.fastq"])

# Dereplicate
derep_fwd = dada2.dereplicate("filtered_R1.fastq")
derep_rev = dada2.dereplicate("filtered_R2.fastq")

# Denoise
result_fwd = dada2.run_dada(derep_fwd, model_fwd)
result_rev = dada2.run_dada(derep_rev, model_rev)

# Merge
merged = dada2.merge_pairs(result_fwd, result_rev, min_overlap=20)

# Remove chimeras
seqs_abund = [(merged[i][0], merged[i][1]) for i in range(len(merged))]
clean = dada2.remove_bimeras(seqs_abund)

# Assign taxonomy (with external lineage TSV)
taxa = dada2.assign_taxonomy(clean, "silva_132.fasta",
                             lineage_tsv="silva_132_lineage.tsv")
for t in taxa:
    print(f"{t.genus} ({t.confidence:.0%})")

# Build sequence table
table = dada2.make_sequence_table(["sample1"], [result_fwd])
table.to_tsv("asv_table.tsv")
```

### Multi-sample pooled workflow

```python
import dada2

samples = ["s1_filtered.fastq", "s2_filtered.fastq", "s3_filtered.fastq"]

# Learn a shared error model from all samples
model = dada2.learn_errors(samples)

# Dereplicate each sample
dereps = [dada2.dereplicate(s) for s in samples]

# Pool all samples for denoising (disk-backed — RAM stays flat)
results = dada2.run_dada_pooled(dereps, model, omega_a=1e-40)

# Build the sequence table
names = ["s1", "s2", "s3"]
table = dada2.make_sequence_table(names, results)
table.to_tsv("asv_table.tsv")
print(table.to_json())
```

---

## Testing

### Rust

```bash
cd dada2_rust
cargo test --workspace
# 24 unit tests (across 11 modules) + 1 integration test
# Integration test: 1 000 synthetic reads, 2% errors → top ASV ≥ 95% identity
```

### Python (requires built extension)

```bash
cd dada2_test
uv run pytest            # 111 tests across 12 test files
uv run pytest --cov      # with coverage
```

Test files:

| File | What it covers |
|---|---|
| `test_basics.py` | import, version, all class/function presence |
| `test_filter.py` | trunc_len, max_ee, min_len, trim_left; 7 paired-end tests |
| `test_error_model.py` | learn_errors, plot_errors, graceful fallback on 1 read |
| `test_derep.py` | count correctness, sort order, total conservation |
| `test_dada.py` | result type, abundance order, omega_a sensitivity; 4 pooled tests |
| `test_merge.py` | return types, non-empty merge, impossible overlap, length bounds |
| `test_chimera.py` | known bimera removed, parents retained, empty input |
| `test_taxonomy.py` | kingdom assignment, confidence range; 3 lineage TSV tests |
| `test_pipeline.py` | full E2E, abundance conservation, determinism, empty-after-filter |
| `test_quality_profile.py` | cycle statistics, percentile ordering, n_reads limit |
| `test_sequence_table.py` | shape, hex encoding, TSV/JSON output, single-sample edge case |
| `test_primer.py` | primer stripping, mismatch tolerance, missing file error |

### Parity test against R dada2

```bash
# Regenerate the reference fixture (requires R + dada2 + jsonlite packages)
Rscript tests/parity/generate_r_output.R

# Compare Rust pipeline output against the R reference
python tests/parity/compare_r_output.py rust_output.json
# Asserts: Jaccard ≥ 0.95, Pearson r ≥ 0.99, no false positives with abundance > 10
```

### Full CI gate

```bash
cd dada2_rust
bash scripts/check.sh
# runs: cargo fmt --check → cargo clippy -D warnings → cargo test
#       → cargo build --release → maturin develop → pytest
```

---

## Configuration reference

### `FilterConfig`

| Field | Default | Description |
|---|---|---|
| `trunc_len` | `0` | Truncate reads to this length; 0 = no truncation |
| `min_len` | `20` | Discard reads shorter than this after truncation |
| `max_ee` | `2.0` | Maximum expected errors (sum of error probabilities) |
| `trunc_q` | `2` | Truncate at first base with Phred below this value |
| `trim_left` | `0` | Bases to remove from the 5′ end |
| `trim_right` | `0` | Bases to remove from the 3′ end |

### `DadaConfig` (Rust only)

| Field | Default | Description |
|---|---|---|
| `omega_a` | `1e-40` | Abundance p-value threshold for accepting a new ASV |
| `pool` | `false` | Enable cross-sample pooling (use `run_dada_pooled` from Python) |
| `max_iter` | `16` | Maximum EM iterations |
| `tol` | `1e-6` | Log-likelihood convergence tolerance |
| `seed` | `42` | RNG seed for reproducibility |

### `TaxonomyConfig` (Rust only)

| Field | Default | Description |
|---|---|---|
| `k` | `8` | k-mer length |
| `threshold` | `0.80` | Minimum bootstrap confidence to report genus-level assignment |
| `seed` | `42` | Bootstrap subsampling seed |

### Lineage TSV format

```
seq_id<TAB>kingdom;phylum;class;order;family;genus;species
```

Compatible with SILVA (`tax_slv_ssu_*.txt`) and GTDB lineage files after minor
column reordering. Pass via `assign_taxonomy(..., lineage_tsv="path/to/lineage.tsv")`.

---

## Key dependencies

| Crate | Version | Purpose |
|---|---|---|
| `pyo3` | 0.23 | Python ↔ Rust FFI |
| `maturin` | ≥ 1.5 | Build system for the Python extension |
| `rayon` | 1 | Work-stealing parallelism |
| `needletail` | 0.5 | Streaming FASTQ/FASTA parser |
| `ndarray` | 0.15 | 2-D error matrix |
| `statrs` | 0.17 | Poisson distribution for abundance p-values |
| `thiserror` | 1 | Ergonomic error types |
| `serde` / `serde_json` | 1 | Serialisation for pool chunks and JSON output |
| `tempfile` | 3 | Temp directory for `PoolStore` disk chunks |

---

## Disk and memory notes

- **Build artifacts** (`target/`): ~2 GB for a full debug + release build. Clean with `cargo clean`.
- **Runtime RAM**: proportional to the number of *unique* sequences per sample, not raw read count. With `run_dada_pooled` and `PoolStore`, even 64 GB datasets run on a 16 GB machine.
- **Temp files** (`PoolStore`): written to the OS temp directory, automatically deleted when `PoolStore` is dropped.
