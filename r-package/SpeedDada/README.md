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

SpeedDada is currently a source-only package. You need a Rust toolchain
(install via [rustup](https://rustup.rs)) and R >= 4.1.

### From Bioconductor (planned)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SpeedDada")
```

### From GitHub (current)

```r
# install.packages("remotes")
remotes::install_github("Genedance/speed-dada",
                        subdir = "r-package/SpeedDada")
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
