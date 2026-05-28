#!/usr/bin/env Rscript
# Synthetic FASTQ fixture generator for cross-platform parity tests.
#
# Generates deterministic paired-end FASTQ files for each supported
# sequencing platform from a small reference of 16S V4-like sequences,
# applying platform-specific quality binning and substitution error rates.
#
# Outputs land under inst/extdata/platforms/<platform>_R{1,2}.fastq.gz.
# Re-running with the same seed yields byte-identical files, so the
# committed fixtures stay reproducible.

set.seed(20250528L)  # date of the speeddada parity work

# --------------------------------------------------------------------------- #
# Reference sequences (5 mock ASVs, 300 bp each)
# --------------------------------------------------------------------------- #
.mkref <- function(seed, len = 300L) {
  set.seed(seed)
  paste(sample(c("A","C","G","T"), len, replace = TRUE), collapse = "")
}

REFS <- list(
  asv1 = .mkref(101L),
  asv2 = .mkref(102L),
  asv3 = .mkref(103L),
  asv4 = .mkref(104L),
  asv5 = .mkref(105L)
)
ABUNDANCES <- c(50L, 40L, 30L, 20L, 10L)  # total reads per ASV

# --------------------------------------------------------------------------- #
# Platform profiles
#
#   quality_bins: vector of integer Phred values the platform reports
#   bin_prob:     probability vector over quality_bins (same length)
#   read_len:     forward / reverse read length (paired-end)
#   notes:        free-text for documentation
# --------------------------------------------------------------------------- #
PLATFORMS <- list(
  miseq = list(
    quality_bins = 0:41,
    bin_prob     = dnorm(0:41, mean = 35, sd = 4),
    read_len     = 250L,
    notes        = "Illumina MiSeq v3, 2x300 mode (250 bp post-trim)"
  ),
  hiseq = list(
    quality_bins = 0:41,
    bin_prob     = dnorm(0:41, mean = 32, sd = 5),
    read_len     = 150L,
    notes        = "Illumina HiSeq 2500, 2x150 mode"
  ),
  novaseq = list(
    # Illumina NovaSeq 6000 default 4-bin Q-table: 2, 12, 23, 37
    quality_bins = c(2L, 12L, 23L, 37L),
    bin_prob     = c(0.02, 0.10, 0.20, 0.68),
    read_len     = 150L,
    notes        = "Illumina NovaSeq 6000, 4-bin binned quality"
  ),
  nextseq = list(
    # NextSeq 500/550/2000 with 8-bin Q-table
    quality_bins = c(2L, 12L, 18L, 25L, 32L, 36L, 38L, 40L),
    bin_prob     = c(0.02, 0.04, 0.06, 0.12, 0.20, 0.26, 0.18, 0.12),
    read_len     = 150L,
    notes        = "Illumina NextSeq 1000/2000, 8-bin binned quality"
  ),
  mgi = list(
    # MGI DNBSEQ-G400 / T7 typical binned profile (12 bins)
    quality_bins = c(2L, 8L, 14L, 18L, 22L, 26L, 30L, 33L, 36L, 38L, 40L, 41L),
    bin_prob     = c(0.02, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.14,
                     0.14, 0.12, 0.07, 0.03),
    read_len     = 200L,
    notes        = "MGI DNBSEQ-G400/T7, binned quality"
  ),
  pacbio_ccs = list(
    quality_bins = c(35L, 38L, 40L, 41L),
    bin_prob     = c(0.05, 0.15, 0.40, 0.40),
    read_len     = 250L,
    notes        = "PacBio CCS / HiFi, near-Q40 (clamped to MAX_QUAL=41 for now)"
  ),
  nanopore = list(
    # ONT R10.4 typical: Phred 5-20, mean ~12
    quality_bins = 5:20,
    bin_prob     = dnorm(5:20, mean = 12, sd = 3),
    read_len     = 300L,  # kept short for fixture size; real ONT is kb+
    notes        = "Oxford Nanopore R10.4 — NOT SUPPORTED, here only for fail-mode tests"
  )
)

# --------------------------------------------------------------------------- #
# Quality / error simulators
# --------------------------------------------------------------------------- #
.sample_qual <- function(n, profile) {
  prob <- profile$bin_prob / sum(profile$bin_prob)
  sample(profile$quality_bins, n, replace = TRUE, prob = prob)
}

.phred_to_ascii <- function(q) {
  intToUtf8(pmin(pmax(as.integer(q), 0L), 93L) + 33L, multiple = TRUE)
}

.mutate_seq <- function(seq_chars, qual_phred) {
  # P(error | q) = 10^(-q/10). On error, draw a uniform non-match base.
  p_err <- 10^(-qual_phred / 10)
  errors <- runif(length(seq_chars)) < p_err
  if (any(errors)) {
    bases <- c("A","C","G","T")
    for (i in which(errors)) {
      alt <- setdiff(bases, seq_chars[i])
      seq_chars[i] <- alt[sample.int(3L, 1L)]
    }
  }
  seq_chars
}

.revcomp <- function(s) {
  bs <- strsplit(s, "", fixed = TRUE)[[1L]]
  comp <- c(A="T", C="G", G="C", T="A", N="N")
  paste(rev(comp[bs]), collapse = "")
}

# --------------------------------------------------------------------------- #
# One platform → paired-end FASTQ files
# --------------------------------------------------------------------------- #
gen_platform <- function(name, profile, refs, abundances, out_dir) {
  rl <- profile$read_len
  total_reads <- sum(abundances)
  fwd_lines <- character(4L * total_reads)
  rev_lines <- character(4L * total_reads)
  idx <- 1L

  for (asv_i in seq_along(refs)) {
    n <- abundances[asv_i]
    ref <- refs[[asv_i]]
    # Forward read = first rl bases. Reverse read = revcomp(last rl bases).
    fwd_template <- substr(ref, 1L, min(rl, nchar(ref)))
    rev_template <- .revcomp(substr(ref, max(1L, nchar(ref) - rl + 1L),
                                    nchar(ref)))
    if (nchar(fwd_template) < rl) {
      fwd_template <- paste0(fwd_template,
                             paste(rep("A", rl - nchar(fwd_template)),
                                   collapse = ""))
    }
    if (nchar(rev_template) < rl) {
      rev_template <- paste0(rev_template,
                             paste(rep("A", rl - nchar(rev_template)),
                                   collapse = ""))
    }
    fwd_chars0 <- strsplit(fwd_template, "", fixed = TRUE)[[1L]]
    rev_chars0 <- strsplit(rev_template, "", fixed = TRUE)[[1L]]

    for (k in seq_len(n)) {
      qF <- .sample_qual(rl, profile)
      qR <- .sample_qual(rl, profile)
      fwd_chars <- .mutate_seq(fwd_chars0, qF)
      rev_chars <- .mutate_seq(rev_chars0, qR)
      id <- sprintf("@%s_asv%d_r%05d", name, asv_i, k)
      fwd_lines[idx]     <- id
      fwd_lines[idx + 1L] <- paste(fwd_chars, collapse = "")
      fwd_lines[idx + 2L] <- "+"
      fwd_lines[idx + 3L] <- paste(.phred_to_ascii(qF), collapse = "")
      rev_lines[idx]     <- id
      rev_lines[idx + 1L] <- paste(rev_chars, collapse = "")
      rev_lines[idx + 2L] <- "+"
      rev_lines[idx + 3L] <- paste(.phred_to_ascii(qR), collapse = "")
      idx <- idx + 4L
    }
  }

  fwd_out <- file.path(out_dir, sprintf("%s_R1.fastq.gz", name))
  rev_out <- file.path(out_dir, sprintf("%s_R2.fastq.gz", name))
  con <- gzfile(fwd_out, "w"); writeLines(fwd_lines, con); close(con)
  con <- gzfile(rev_out, "w"); writeLines(rev_lines, con); close(con)
  invisible(list(fwd = fwd_out, rev = rev_out, reads = total_reads))
}

# --------------------------------------------------------------------------- #
# Driver: regenerate every platform fixture
# --------------------------------------------------------------------------- #
out_dir <- file.path(
  if (nzchar(Sys.getenv("SPEEDDADA_FIXTURE_DIR")))
    Sys.getenv("SPEEDDADA_FIXTURE_DIR")
  else
    "inst/extdata/platforms"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (name in names(PLATFORMS)) {
  res <- gen_platform(name, PLATFORMS[[name]], REFS, ABUNDANCES, out_dir)
  cat(sprintf("[ok] %-12s %5d reads -> %s\n", name, res$reads,
              dirname(res$fwd)))
}

# Drop a README in the platforms dir explaining provenance.
writeLines(c(
  "# Synthetic platform fixtures",
  "",
  "Generated by `inst/scripts/gen-fixture.R` from a fixed seed.",
  "Re-run with `Rscript inst/scripts/gen-fixture.R` from the package root.",
  "",
  "Each *_R{1,2}.fastq.gz pair contains 150 paired reads drawn from 5",
  "mock 16S V4-like reference ASVs with platform-specific quality binning",
  "and substitution-rate-only error injection. The `nanopore` fixture is",
  "included only to verify SpeedDada's failure modes — ONT is not yet",
  "supported (indel-aware alignment is a separate body of work)."
),
con = file.path(out_dir, "README.md"))
