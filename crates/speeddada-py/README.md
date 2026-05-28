# speeddada

[![PyPI](https://img.shields.io/pypi/v/speeddada.svg)](https://pypi.org/project/speeddada/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**High-performance DADA2 amplicon sequence variant pipeline — Rust core with Python bindings.**

A drop-in replacement for the Python ecosystem's call to DADA2
(Callahan et al. 2016, [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869))
implemented in Rust. Typically **10–20× faster** than R `dada2` and uses
**~10× less peak memory** through streaming dereplication, disk-backed
pooling, and compact bitset taxonomy profiles.

## Installation

Distributed via **GitHub Releases** for now. PyPI publishing
(`pip install speeddada`) will follow in a later release once the
Trusted Publisher is set up. Pre-built wheels are attached to every
tagged release for Linux (x86_64 / aarch64 / i686), macOS
(universal2: Intel + Apple Silicon), and Windows x64.

**Pre-built wheel (no Rust toolchain needed):** pick the matching wheel
from the [latest release](https://github.com/Genedance/speed-dada/releases/latest)
and install it directly:

```bash
# example: linux x86_64, cpython 3.12 — substitute your platform tag
pip install "https://github.com/Genedance/speed-dada/releases/download/v0.99.1/speeddada-0.99.1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
```

**Source install (any platform, needs Rust):**

```bash
pip install "git+https://github.com/Genedance/speed-dada.git@v0.99.1#subdirectory=crates/speeddada-py"
```

Building from source requires a Rust toolchain (install via
[rustup](https://rustup.rs)) and Python 3.9+.

## Quick start

```python
import speeddada

# 1. (optional) configure rayon thread pool
speeddada.configure_runtime(n_threads=8)

# 2. quality filtering — paired-end
cfg_fwd = speeddada.FilterConfig(trunc_len=240, max_ee=2.0, min_len=50)
cfg_rev = speeddada.FilterConfig(trunc_len=200, max_ee=5.0, min_len=50)
stats = speeddada.filter_and_trim_paired(
    cfg_fwd, cfg_rev,
    "samples/R1.fastq.gz", "samples/R2.fastq.gz",
    "filtered/R1.fastq.gz", "filtered/R2.fastq.gz",
)
print(f"reads_in={stats.reads_in}, pairs_out={stats.pairs_out}")

# 3. learn error model
err = speeddada.learn_errors(["filtered/R1.fastq.gz"], n_reads=1_000_000)

# 4. dereplicate + DADA denoising (per sample)
derep_fwd = speeddada.derep_fastq("filtered/R1.fastq.gz")
asvs_fwd  = speeddada.dada(derep_fwd, err, omega_a=1e-40)

# 5. merge pairs + remove chimeras
derep_rev = speeddada.derep_fastq("filtered/R2.fastq.gz")
asvs_rev  = speeddada.dada(derep_rev, err, omega_a=1e-40)
merged    = speeddada.merge_pairs(asvs_fwd, asvs_rev, min_overlap=20)
clean     = speeddada.remove_bimera_denovo([(m.sequence, m.abundance) for m in merged])
```

The shorter `run_pipeline` / `run_pipeline_samples` helpers run the full
pipeline end-to-end if you prefer a single call.

## API reference

Every function and class is type-annotated; your editor / type-checker
will pick up the bundled `.pyi` stubs automatically. The full API
reference is at <https://genedance.github.io/speed-dada/>.

## Cross-platform support

| OS / Arch         | Wheel | Source build |
|-------------------|:-----:|:------------:|
| Linux x86_64      | ✅    | ✅           |
| Linux aarch64     | ✅    | ✅           |
| macOS x86_64      | ✅    | ✅           |
| macOS arm64       | ✅    | ✅           |
| Windows x86_64    | ✅    | ✅           |
| Raspberry Pi (aarch64) | ✅ (Linux wheel) | ✅ |

## License

MIT © 2026 Genedance GmbH. Author: Alexandre Jousset
(<info@genedance.com>).

## Citation

If you use speeddada in a publication, please cite the original DADA2 paper
(Callahan et al. 2016, *Nature Methods*) and reference this package via its
[`CITATION.cff`](https://github.com/Genedance/speed-dada/blob/main/CITATION.cff).
