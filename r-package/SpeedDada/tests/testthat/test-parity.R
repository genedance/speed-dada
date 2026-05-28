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

# Acceptance thresholds per platform. Binned-quality platforms get a
# looser taxonomy bar because LOESS / piecewise-linear divergence has a
# small downstream effect on borderline ASV assignments.
PARITY_THRESHOLDS <- list(
  miseq      = list(tax = 0.99, err = 5e-3, asv_jaccard = 0.98, abundance = 1L),
  hiseq      = list(tax = 0.99, err = 5e-3, asv_jaccard = 0.98, abundance = 1L),
  novaseq    = list(tax = 0.95, err = 2e-2, asv_jaccard = 0.95, abundance = 2L),
  nextseq    = list(tax = 0.95, err = 2e-2, asv_jaccard = 0.95, abundance = 2L),
  mgi        = list(tax = 0.95, err = 2e-2, asv_jaccard = 0.95, abundance = 2L),
  pacbio_ccs = list(tax = 0.95, err = 5e-2, asv_jaccard = 0.90, abundance = 2L)
)

run_speeddada_pipeline <- function(fwd, rev, ref = NULL, species = NULL,
                                   trunc_len = c(150L, 150L)) {
  filt_fwd <- tempfile(fileext = ".fastq")
  filt_rev <- tempfile(fileext = ".fastq")
  ft <- filterAndTrim(fwd, filt_fwd, rev, filt_rev,
                      truncLen = trunc_len, maxEE = c(2, 2), minLen = 50L)
  errF <- learnErrors(filt_fwd, nbases = 1e6)
  errR <- learnErrors(filt_rev, nbases = 1e6)
  dF <- derepFastq(filt_fwd); dR <- derepFastq(filt_rev)
  aF <- dada(dF, errF); aR <- dada(dR, errR)
  m <- mergePairs(aF, dF, aR, dR)
  seqtab <- makeSequenceTable(list(s1 = m))
  seqtab_nc <- removeBimeraDenovo(seqtab)
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

  thresholds <- PARITY_THRESHOLDS[[platform]]
  sd_out <- run_speeddada_pipeline(fix$fwd, fix$rev)

  ref <- if (parity_enabled() && dada2_available()) {
    run_dada2_subprocess(fix$fwd, fix$rev)
  } else {
    load_dada2_snapshot(platform)
  }
  if (is.null(ref))
    skip(sprintf(
      "no dada2 reference for '%s' (set SPEEDDADA_PARITY=1 with dada2 installed, or bake a snapshot under _snaps/dada2_%s.rds)",
      platform, platform))

  if (!is.null(ref$seqtab_nc)) {
    cmp <- compare_seqtabs(sd_out$seqtab_nc, ref$seqtab_nc)
    expect_gte(cmp$asv_jaccard, thresholds$asv_jaccard,
               label = sprintf("[%s] ASV set Jaccard", platform))
    if (!is.na(cmp$max_abs_diff))
      expect_lte(cmp$max_abs_diff, thresholds$abundance,
                 label = sprintf("[%s] per-cell abundance diff", platform))
  }
  if (!is.null(ref$tax) && !is.null(sd_out$tax)) {
    tax_cmp <- compare_taxonomy_matrices(sd_out$tax, ref$tax,
                                         levels = c("Kingdom","Phylum","Class",
                                                    "Order","Family","Genus"))
    expect_gte(tax_cmp$agreement, thresholds$tax,
               label = sprintf("[%s] taxonomy cell agreement", platform))
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
