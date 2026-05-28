#!/usr/bin/env Rscript
# Bake reference dada2 outputs as committed RDS snapshots.
#
# Run this once on a machine where Bioconductor `dada2` is installed in
# the same library tree. The resulting fixtures-dada2/*.rds files are
# committed alongside the synthetic FASTQ fixtures, so future parity test
# runs (without dada2 installed) still get meaningful diffs against a
# fixed reference.
#
# Usage:
#   cd r-package/SpeedDada
#   Rscript inst/scripts/bake-dada2-snapshots.R

if (!requireNamespace("dada2", quietly = TRUE))
  stop("dada2 is required to bake snapshots; install via BiocManager::install(\"dada2\")")
if (!requireNamespace("SpeedDada", quietly = TRUE))
  stop("install SpeedDada first (R CMD INSTALL SpeedDada)")

source("tests/testthat/helper-parity.R")

PARITY_PARAMS <- list(
  miseq      = list(trunc = c(150L, 150L), maxEE = c(5, 5)),
  hiseq      = list(trunc = c(150L, 150L), maxEE = c(5, 5)),
  novaseq    = list(trunc = c(150L, 150L), maxEE = c(8, 8)),
  nextseq    = list(trunc = c(150L, 150L), maxEE = c(8, 8)),
  mgi        = list(trunc = c(150L, 150L), maxEE = c(8, 8)),
  pacbio_ccs = list(trunc = c(150L, 150L), maxEE = c(5, 5))
)

out_dir <- "tests/testthat/fixtures-dada2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (p in names(PARITY_PARAMS)) {
  fix <- .parity_fixture(p)
  if (is.null(fix)) {
    message(sprintf("[skip] %s: fixture not found", p)); next
  }
  params <- PARITY_PARAMS[[p]]
  cat(sprintf("[bake] %s ...\n", p))
  res <- run_dada2_subprocess(fix$fwd, fix$rev,
                              trunc_len = params$trunc, max_ee = params$maxEE)
  saveRDS(res, file.path(out_dir, sprintf("dada2_%s.rds", p)))
}

writeLines(c(
  "# dada2 reference snapshots",
  "",
  sprintf("Baked %s against `dada2` %s by `inst/scripts/bake-dada2-snapshots.R`.",
          format(Sys.Date()),
          as.character(packageVersion("dada2"))),
  "",
  "Re-bake whenever the synthetic fixtures change (`inst/scripts/gen-fixture.R`)",
  "or whenever a new dada2 release lands."
),
con = file.path(out_dir, "README.md"))

cat("done.\n")
