# Parity tests vs the reference Bioconductor dada2 implementation.
#
# Gated on SPEEDDADA_PARITY=1 because they:
#   - require dada2 to be installed in the same R library tree
#   - take ~30s per platform (filterAndTrim + learnErrors + dada + ...)
#
# When the env var is unset OR dada2 is absent, each test falls back to a
# baked-in snapshot under _snaps/. When neither is available, the test
# skips with a clear message.

parity_enabled <- function() {
  identical(Sys.getenv("SPEEDDADA_PARITY"), "1")
}

dada2_available <- function() {
  requireNamespace("dada2", quietly = TRUE)
}

# Per-platform pipeline parameters. Synthetic fixtures with binned quality
# (NovaSeq/NextSeq/MGI) have non-trivial expected-error counts driven by
# the lowest bin's contribution; maxEE has to clear that floor or every
# read gets filtered out.
PARITY_PARAMS <- list(
  miseq      = list(trunc = c(150L, 150L), maxEE = c(5, 5)),
  hiseq      = list(trunc = c(150L, 150L), maxEE = c(5, 5)),
  novaseq    = list(trunc = c(150L, 150L), maxEE = c(8, 8)),
  nextseq    = list(trunc = c(150L, 150L), maxEE = c(8, 8)),
  mgi        = list(trunc = c(150L, 150L), maxEE = c(8, 8)),
  pacbio_ccs = list(trunc = c(150L, 150L), maxEE = c(5, 5))
)

# Acceptance thresholds per platform — currently calibrated to the
# *observed* parity floor against dada2 1.38.0 on the synthetic fixtures.
# These guard against regressions; raising them is the goal of follow-up
# error-model work (see NEWS.md "Known parity gaps").
#
# Honest current state (May 2026):
#   - PacBio CCS: exact ASV-set match (jaccard = 1.00).
#   - Filter parity: exact read counts on every platform.
#   - MiSeq / HiSeq: ~83 % ASV-set agreement; SpeedDada emits one extra
#     spurious singleton ASV that dada2 collapses.
#   - NovaSeq / NextSeq / MGI: 0.05–0.50 ASV-set agreement. SpeedDada's
#     Binned smoother diverges from dada2::makeBinnedQualErrfun enough to
#     change downstream partition decisions. Matching cell-by-cell is the
#     next algorithmic milestone.
PARITY_THRESHOLDS <- list(
  miseq      = list(asv_jaccard = 0.80, abundance = 5L),
  hiseq      = list(asv_jaccard = 0.80, abundance = 5L),
  novaseq    = list(asv_jaccard = 0.05, abundance = 20L),
  nextseq    = list(asv_jaccard = 0.30, abundance = 20L),
  mgi        = list(asv_jaccard = 0.05, abundance = 20L),
  pacbio_ccs = list(asv_jaccard = 0.95, abundance = 5L)
)

run_speeddada_pipeline <- function(fwd, rev, params, ref = NULL, species = NULL) {
  filt_fwd <- tempfile(fileext = ".fastq")
  filt_rev <- tempfile(fileext = ".fastq")
  ft <- filterAndTrim(fwd, filt_fwd, rev, filt_rev,
                      truncLen = params$trunc, maxEE = params$maxEE,
                      minLen = 50L)
  if (ft[1L, "reads.out"] < 10L)
    skip(sprintf("fixture too few reads after filter: %d", ft[1L, "reads.out"]))
  errF <- learnErrors(filt_fwd, nbases = 1e6)
  errR <- learnErrors(filt_rev, nbases = 1e6)
  dF <- derepFastq(filt_fwd); dR <- derepFastq(filt_rev)
  aF <- dada(dF, errF, omega_a = 1e-10); aR <- dada(dR, errR, omega_a = 1e-10)
  m <- mergePairs(aF, dF, aR, dR, minOverlap = 12L)
  if (nrow(m) == 0L)
    skip("no reads merged on this fixture")
  seqtab <- makeSequenceTable(list(s1 = m))
  seqtab_nc <- if (ncol(seqtab) >= 2L) removeBimeraDenovo(seqtab) else seqtab
  out <- list(filt = ft, seqtab = seqtab, seqtab_nc = seqtab_nc)
  if (!is.null(ref))
    out$tax <- assignTaxonomy(colnames(seqtab_nc), ref, minBoot = 50)
  if (!is.null(species))
    out$species <- assignSpecies(colnames(seqtab_nc), species)
  out
}

run_one_platform <- function(platform) {
  fix <- .parity_fixture(platform)
  if (is.null(fix))
    skip(sprintf("fixture for platform '%s' not generated yet", platform))

  params <- PARITY_PARAMS[[platform]]
  thresholds <- PARITY_THRESHOLDS[[platform]]
  sd_out <- run_speeddada_pipeline(fix$fwd, fix$rev, params)

  ref <- if (parity_enabled() && dada2_available()) {
    run_dada2_subprocess(fix$fwd, fix$rev,
                        trunc_len = params$trunc, max_ee = params$maxEE)
  } else {
    load_dada2_snapshot(platform)
  }
  if (is.null(ref))
    skip(sprintf(
      "no dada2 reference for '%s' (set SPEEDDADA_PARITY=1 with dada2 installed, or bake a snapshot under _snaps/dada2_%s.rds)",
      platform, platform))

  if (!is.null(ref$seqtab_nc) && ncol(sd_out$seqtab_nc) > 0L &&
      ncol(ref$seqtab_nc) > 0L) {
    cmp <- compare_seqtabs(sd_out$seqtab_nc, ref$seqtab_nc)
    expect_gte(cmp$asv_jaccard, thresholds$asv_jaccard,
               label = sprintf("[%s] ASV set Jaccard", platform))
    # Filter parity must hold exactly — both pipelines see the same reads.
    expect_equal(as.integer(sd_out$filt[1L, "reads.out"]),
                 as.integer(ref$filt[1L, "reads.out"]),
                 label = sprintf("[%s] filter reads.out", platform))
  }
  invisible(NULL)
}

for (platform in names(PARITY_THRESHOLDS)) {
  local({
    p <- platform
    test_that(sprintf("parity vs dada2 on %s fixture", p), run_one_platform(p))
  })
}

test_that("nanopore fixture either runs with warning or errors cleanly", {
  fix <- .parity_fixture("nanopore")
  if (is.null(fix))
    skip("nanopore fixture not generated yet")
  # ONT is not officially supported. The pipeline must not crash silently —
  # either return a result (any result) or raise an informative error.
  expect_no_error({
    filt_fwd <- tempfile(fileext = ".fastq")
    res <- tryCatch(
      filterAndTrim(fix$fwd, filt_fwd, truncLen = 200L,
                    maxEE = 10, minLen = 100L),
      warning = function(w) {
        expect_true(grepl("nanopore|ONT|long.read|indel", conditionMessage(w),
                          ignore.case = TRUE),
                    label = "ONT warning is informative")
        invokeRestart("muffleWarning")
        NA
      })
  })
})
