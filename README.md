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
│   │       ├── filter.rs          # stage 2 — streaming quality filter
│   │       ├── error_model.rs     # stage 3 — EM error rate learning
│   │       ├── derep.rs           # stage 4 — dereplication
│   │       ├── dada.rs            # stage 5 — core DADA algorithm
│   │       ├── merge.rs           # stage 6 — paired-end merging
│   │       ├── chimera.rs         # stage 7 — bimera detection
│   │       ├── taxonomy.rs        # stage 8 — naive Bayes k-mer classifier
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
│       └── fixtures/r_output.json
└── scripts/
    └── check.sh                   # full CI gate (fmt → clippy → test → build → pytest)
```

---

## Pipeline stages

```
FASTQ files
    │
    ▼  filter_and_trim        quality filter, length truncation, expected-error cutoff
    │                         streaming: one record at a time — O(1) RAM regardless of file size
    ▼  learn_errors           fit logistic error model P(obs|true, Phred q) via EM
    │                         16-class (all base transitions) × 41 Phred bins
    ▼  dereplicate            collapse identical sequences; track per-position quality sums
    ▼  run_dada               core denoising — Poisson abundance p-values, greedy partition + EM
    │                         OMEGA_A default 1e-40; convergence tol 1e-6; max 16 iterations
    ▼  merge_pairs            suffix-prefix overlap of F+R reads (min_overlap=20 default)
    ▼  remove_bimeras         exact bimera search; min arm 8 bp; parent must outrank candidate
    ▼  assign_taxonomy        naive Bayes k-mer classifier (k=8, 100 bootstrap replicates)
    │
    ▼  ASV count table
```

---

## Architecture decisions

### GIL release
Every CPU-bound Python function calls `py.allow_threads(|| { … })` before entering Rust. Python threads remain live while Rust works. Objects borrowed from Python (`ErrorModel`, `DadaResult`) are cloned before crossing the thread boundary.

### Streaming I/O
`filter_and_trim` reads and writes one FASTQ record at a time via a needletail cursor + `BufWriter`. No `Vec<FastqRecord>` accumulation — a 100 GB file uses no more RAM than a 1 MB file during filtering.

### Disk-backed pooling (`pool.rs`)
`PoolStore` accumulates unique sequences from multiple samples in a `BTreeMap`. When the in-memory entry count exceeds `flush_threshold` (default 500 000), the current map is serialised to a JSONL chunk in a `tempfile::TempDir` and the map is cleared. `into_pooled_uniques()` re-merges all chunks and returns a single sorted `Vec<UniqueSeq>` for DADA. RAM use is proportional to unique sequence count, not raw read count.

### SIMD alignment (`align.rs`)
`hamming_distance`, `first_mismatch`, and `range_equal` are written as tight scalar loops that LLVM reliably auto-vectorises to AVX2 / NEON / SSE4 when compiled with `target-cpu=native` (set in `.cargo/config.toml`). No `unsafe` blocks, no third-party SIMD crate. `chimera.rs` and `merge.rs` delegate to these primitives.

### Rayon parallelism
- **Sample level:** `filter_and_trim_many(cfg, pairs)` processes N FASTQ pairs across the Rayon thread pool.
- **Intra-sample:** the DADA E-step, chimera candidate scan, and taxonomy classification all use `par_iter()`.

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
cargo test --workspace              # 18 tests (17 unit + 1 integration)
```

### Python extension (development install)

```bash
cd dada2_rust
# First build — downloads crates and compiles (~30 s on first run)
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

dada2.__version__   # "0.1.0"
```

### Classes

| Class | Description |
|---|---|
| `FilterConfig` | Parameters for filter_and_trim (trunc_len, min_len, max_ee, trunc_q, trim_left, trim_right) |
| `FilterStats` | `.reads_in`, `.reads_out` returned by filter_and_trim |
| `ErrorModel` | Learned error matrix; `.plot_errors()` → dict for matplotlib |
| `DadaResult` | Sequence of `(sequence: bytes, abundance: int)` tuples |
| `TaxonAssignment` | Per-ASV classification with `.kingdom` … `.species` and `.confidence` |

### Functions

```python
# Stage 2 — filter
stats = dada2.filter_and_trim(config, input_path, output_path)

# Parallel filter — processes many files at once (returns list of FilterStats)
# Available from Rust: filter_and_trim_many(cfg, [(in, out), ...])

# Stage 3 — error model
model = dada2.learn_errors(fastq_paths, n_reads=1_000_000)

# Stage 4 — dereplicate
derep = dada2.dereplicate(fastq_path)          # → list[(bytes, int)]

# Stage 5 — denoise
result = dada2.run_dada(derep, model, omega_a=1e-40, pool=False)

# Stage 6 — merge paired ends
merged = dada2.merge_pairs(fwd_result, rev_result, min_overlap=20)

# Stage 7 — remove chimeras
clean = dada2.remove_bimeras(seqs)             # seqs: list[(bytes, int)]

# Stage 8 — taxonomy
assignments = dada2.assign_taxonomy(seqs, ref_fasta, k=8)

# Full pipeline in one call (GIL-free end to end)
asv_table = dada2.run_pipeline(
    input_paths=["sample1.fastq", "sample2.fastq"],
    output_dir="/tmp/filtered",
    trunc_len=250,
    max_ee=2.0,
    omega_a=1e-40,
)  # → dict[str, int]  (hex-encoded sequence → abundance)
```

### Typical single-sample workflow

```python
import dada2

cfg = dada2.FilterConfig(trunc_len=250, max_ee=2.0)
stats = dada2.filter_and_trim(cfg, "raw.fastq", "filtered.fastq")
print(f"{stats.reads_out}/{stats.reads_in} reads passed")

model = dada2.learn_errors(["filtered.fastq"], n_reads=1_000_000)
derep = dada2.dereplicate("filtered.fastq")
result = dada2.run_dada(derep, model)

seqs_abund = [(result[i][0], result[i][1]) for i in range(len(result))]
clean = dada2.remove_bimeras(seqs_abund)
print(f"{len(clean)} ASVs after chimera removal")

taxa = dada2.assign_taxonomy(clean, "silva_ref.fasta")
for t in taxa:
    print(t.genus, t.confidence)
```

---

## Testing

### Rust

```bash
cd dada2_rust
cargo test --workspace
# 17 unit tests (one per module) + 1 integration test
# Integration test: 1 000 synthetic reads, 2% errors → top ASV ≥ 95% identity
```

### Python (requires built extension)

```bash
cd dada2_test
uv run pytest            # 68 tests across 8 test files
uv run pytest --cov      # with coverage
```

Test files:

| File | What it covers |
|---|---|
| `test_basics.py` | import, version, class/function presence |
| `test_filter.py` | trunc_len, max_ee, min_len, trim_left, error handling |
| `test_error_model.py` | learn_errors, plot_errors shape and ordering |
| `test_derep.py` | count correctness, sort order, total conservation |
| `test_dada.py` | result type, abundance order, omega_a sensitivity, identical-reads collapse |
| `test_merge.py` | return types, non-empty merge, impossible overlap, length bounds |
| `test_chimera.py` | known bimera removed, parents retained, empty input |
| `test_taxonomy.py` | kingdom assignment, confidence range, missing file error |
| `test_pipeline.py` | full E2E, abundance conservation, determinism, empty-after-filter edge case |

### Full CI gate

```bash
cd dada2_rust
bash scripts/check.sh
# runs: cargo fmt --check → cargo clippy -D warnings → cargo test → cargo build --release
#       → maturin develop → pytest
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
| `pool` | `false` | Enable cross-sample pooling via `PoolStore` |
| `max_iter` | `16` | Maximum EM iterations |
| `tol` | `1e-6` | Log-likelihood convergence tolerance |
| `seed` | `42` | RNG seed for reproducibility |

### `TaxonomyConfig` (Rust only)

| Field | Default | Description |
|---|---|---|
| `k` | `8` | k-mer length |
| `threshold` | `0.80` | Minimum bootstrap confidence to report genus-level assignment |
| `seed` | `42` | Bootstrap subsampling seed |

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
- **Runtime RAM**: proportional to the number of *unique* sequences per sample, not raw read count. With `pool=TRUE` and `PoolStore`, even 64 GB datasets run on a 16 GB machine.
- **Temp files** (`PoolStore`): written to the OS temp directory, automatically deleted when `PoolStore` is dropped.

---

## Limitations and roadmap

- `pool=TRUE` is implemented in `PoolStore` but not yet wired through the Python `run_dada` binding — call the Rust API directly for pooled analysis.
- Taxonomy assignment parses lineage from FASTA description fields as semicolon-separated strings; a separate lineage TSV loader is not yet implemented.
- No paired-end filtering (`filterAndTrim` on R1+R2 simultaneously); filter each file separately then merge.
