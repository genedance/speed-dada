#!/usr/bin/env Rscript
# Minimal end-to-end test for SpeedDada drop-in compatibility.
# Mirrors the bench_r.R pipeline using the Rust-backed package.

.libPaths(c(path.expand("~/R/library"), .libPaths()))
library(SpeedDada)

# ── Generate synthetic paired FASTQs ─────────────────────────────────────────

TRUE_SEQ <- paste(rep(c("A","C","G","T"), 30), collapse="")  # 120 bp

rc <- function(s) {
  comp <- chartr("ACGT", "TGCA", s)
  paste(rev(strsplit(comp, "")[[1]]), collapse="")
}

cycle_err <- function(pos, len, base=0.001, tail=0.05)
  base + (tail - base) * (pos / (len - 1))^2

phred_char <- function(err) {
  q <- pmin(40, pmax(2, round(-10 * log10(pmax(err, 1e-10)))))
  rawToChar(as.raw(q + 33L))
}

make_reads <- function(seq, n, seed) {
  set.seed(seed)
  L    <- nchar(seq)
  errs <- sapply(seq_len(L) - 1, cycle_err, len = L)
  qual <- paste(sapply(errs, phred_char), collapse = "")
  bases_mat <- strsplit(rep(seq, n), "")
  lines <- character(4 * n)
  for (i in seq_len(n)) {
    bases <- bases_mat[[i]]
    for (j in seq_len(L))
      if (runif(1) < errs[j])
        bases[j] <- sample(setdiff(c("A","C","G","T"), bases[j]), 1L)
    lines[4*(i-1)+1] <- paste0("@read_", i-1)
    lines[4*(i-1)+2] <- paste(bases, collapse="")
    lines[4*(i-1)+3] <- "+"
    lines[4*(i-1)+4] <- qual
  }
  lines
}

tmpdir <- tempdir()
r1 <- file.path(tmpdir, "R1.fastq")
r2 <- file.path(tmpdir, "R2.fastq")
writeLines(make_reads(TRUE_SEQ, 500L, 42L), r1)
writeLines(make_reads(rc(TRUE_SEQ), 500L, 43L), r2)

fr1 <- file.path(tmpdir, "filt_R1.fastq")
fr2 <- file.path(tmpdir, "filt_R2.fastq")

# ── 1. filterAndTrim ─────────────────────────────────────────────────────────

cat("── filterAndTrim\n")
stats <- filterAndTrim(r1, fr1, r2, fr2,
                       truncLen = c(100L, 100L), maxEE = c(3, 3), minLen = 50L)
stopifnot(is.matrix(stats))
stopifnot(identical(colnames(stats), c("reads.in", "reads.out")))
stopifnot(stats[1, "reads.in"] == 500L)
stopifnot(stats[1, "reads.out"] > 0L)
cat("  reads.in =", stats[1,1], " reads.out =", stats[1,2], "\n")

# ── 2. learnErrors ───────────────────────────────────────────────────────────

cat("── learnErrors\n")
errF <- learnErrors(fr1)
errR <- learnErrors(fr2)
stopifnot(!is.null(errF))
cat("  errF class =", class(errF), "\n")

# ── 3. derepFastq ────────────────────────────────────────────────────────────

cat("── derepFastq\n")
derepF <- derepFastq(fr1)
derepR <- derepFastq(fr2)
stopifnot(inherits(derepF, "derep"))
stopifnot(is.integer(derepF$uniques) && !is.null(names(derepF$uniques)))
cat("  fwd unique seqs =", length(derepF$uniques),
    " rev =", length(derepR$uniques), "\n")

# ── 4. dada ──────────────────────────────────────────────────────────────────

cat("── dada\n")
dadaF <- dada(derepF, err = errF)
dadaR <- dada(derepR, err = errR)
stopifnot(inherits(dadaF, "dada"))
stopifnot(is.integer(dadaF$denoised) && !is.null(names(dadaF$denoised)))
cat("  fwd ASVs =", length(dadaF$denoised),
    " rev ASVs =", length(dadaR$denoised), "\n")

# ── 5. mergePairs ─────────────────────────────────────────────────────────────

cat("── mergePairs\n")
merged <- mergePairs(dadaF, derepF, dadaR, derepR)
stopifnot(is.data.frame(merged))
stopifnot(all(c("sequence","abundance","accept") %in% names(merged)))
cat("  merged pairs =", sum(merged$accept), "\n")

# ── 6. makeSequenceTable ─────────────────────────────────────────────────────

cat("── makeSequenceTable\n")
seqtab <- makeSequenceTable(list(s1 = merged))
stopifnot(is.matrix(seqtab) && is.integer(seqtab))
stopifnot(nrow(seqtab) == 1L)
cat("  samples =", nrow(seqtab), " ASVs =", ncol(seqtab), "\n")
cat("  top ASV:", substr(colnames(seqtab)[1], 1, 30), "...\n")

# ── 7. removeBimeraDenovo ────────────────────────────────────────────────────

cat("── removeBimeraDenovo\n")
seqtab_clean <- removeBimeraDenovo(seqtab)
stopifnot(is.matrix(seqtab_clean) && is.integer(seqtab_clean))
stopifnot(ncol(seqtab_clean) <= ncol(seqtab))
cat("  ASVs after chimera removal =", ncol(seqtab_clean), "\n")

# ── Verify top ASV matches true sequence ──────────────────────────────────────

top <- colnames(seqtab_clean)[1]
ref <- substr(TRUE_SEQ, 1L, nchar(top))
matches <- sum(strsplit(top,"")[[1]] == strsplit(ref,"")[[1]])
identity <- matches / nchar(top)
cat("\nTop ASV identity to true sequence:", round(identity, 3), "\n")
stopifnot(identity >= 0.90)

cat("\nAll tests passed.\n")
