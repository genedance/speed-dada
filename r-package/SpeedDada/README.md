# SpeedDada

**High-performance drop-in replacement for the DADA2 amplicon pipeline,
backed by a Rust core.**

SpeedDada exposes the seven core dada2 functions
(`filterAndTrim`, `learnErrors`, `derepFastq`, `dada`, `mergePairs`,
`makeSequenceTable`, `removeBimeraDenovo`) with the same API as the
original [dada2](https://github.com/benjjneb/dada2) package
(Callahan et al. 2016, [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869))
but powered by a Rust implementation that is typically **10-20x faster** and
uses **~10x less peak memory**.

## Installation

Prebuilt binary packages are attached to each tagged release for the
common desktop platforms. Installing them does **not** invoke the
Rust compiler on your machine.

### Prebuilt binaries (recommended — no Rust required)

```r
# macOS Apple Silicon (M-series Macs)
install.packages(
  "https://github.com/Genedance/speed-dada/releases/download/v0.99.1/SpeedDada_0.99.1.tgz",
  repos = NULL
)

# Windows x64
install.packages(
  "https://github.com/Genedance/speed-dada/releases/download/v0.99.1/SpeedDada_0.99.1.zip",
  repos = NULL
)

# Linux x86_64
install.packages(
  "https://github.com/Genedance/speed-dada/releases/download/v0.99.1/SpeedDada_0.99.1_R_x86_64-pc-linux-gnu.tar.gz",
  repos = NULL
)
```

Requires R ≥ 4.1.

### From source (Intel Mac, aarch64 Linux, or any other platform)

The source tarball compiles the Rust core during installation, so you
need:

* R >= 4.1
* A Rust toolchain (install via [rustup](https://rustup.rs))
* `libbz2-dev` (Linux only — `apt install libbz2-dev` / `dnf install bzip2-devel`)

```r
install.packages(
  "https://github.com/Genedance/speed-dada/releases/download/v0.99.1/SpeedDada_0.99.1.tar.gz",
  repos = NULL, type = "source"
)

# Or from GitHub at a tag
# install.packages("remotes")
remotes::install_github("Genedance/speed-dada",
                        ref    = "v0.99.1",
                        subdir = "r-package/SpeedDada")
```

### Future channels (planned)

```r
# Bioconductor — after initial testing window
BiocManager::install("SpeedDada")
```

## Quick start

```r
library(SpeedDada)

fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
rev <- system.file("extdata", "sample_R2.fastq", package = "SpeedDada")

# 1. quality filter
fwd_filt <- tempfile(fileext = ".fastq")
rev_filt <- tempfile(fileext = ".fastq")
filterAndTrim(fwd, fwd_filt, rev, rev_filt,
              truncLen = c(30, 30), maxEE = c(5, 5), minLen = 10)

# 2. learn errors
err <- learnErrors(c(fwd_filt, rev_filt), nbases = 1e4)

# 3. dereplicate + denoise
dF <- derepFastq(fwd_filt); dR <- derepFastq(rev_filt)
aF <- dada(dF, err, omega_a = 1e-5)
aR <- dada(dR, err, omega_a = 1e-5)

# 4. merge pairs, build table, remove chimeras
merged <- mergePairs(aF, dF, aR, dR, minOverlap = 5)
seqtab <- makeSequenceTable(list(s1 = merged))
seqtab_nochim <- removeBimeraDenovo(seqtab)
```

The full end-to-end workflow is in the vignette:

```r
vignette("SpeedDada-pipeline", "SpeedDada")
```

## Cross-platform support

SpeedDada builds out of the box on:

* Linux x86_64 / aarch64 (including Raspberry Pi 5)
* macOS x86_64 / arm64 (Apple Silicon)
* Windows x86_64 (Rtools 4.x)

## Citing

If you use SpeedDada in a publication, please cite the original DADA2
paper *and* this package:

```r
citation("SpeedDada")
```

## License

MIT (c) 2026 Genedance GmbH. Author: Alexandre Jousset
(<info@genedance.com>).
