# Parity-test helpers — shared by the SPEEDDADA_PARITY=1-gated tests.
#
# These functions run the reference dada2 implementation in a subprocess
# (so SpeedDada and dada2 don't fight over their respective namespaces in
# the same R session) and compare the results to SpeedDada's output.

# Resolve a per-platform fixture pair. Looks first in inst/extdata/platforms
# (committed), then in the user cache (regenerated on demand).
.parity_fixture <- function(platform) {
  base <- system.file("extdata", "platforms", package = "SpeedDada")
  fwd <- file.path(base, sprintf("%s_R1.fastq.gz", platform))
  rev <- file.path(base, sprintf("%s_R2.fastq.gz", platform))
  if (file.exists(fwd) && file.exists(rev))
    return(list(fwd = fwd, rev = rev))

  cache <- tools::R_user_dir("SpeedDada", "cache")
  fwd <- file.path(cache, sprintf("%s_R1.fastq.gz", platform))
  rev <- file.path(cache, sprintf("%s_R2.fastq.gz", platform))
  if (file.exists(fwd) && file.exists(rev))
    return(list(fwd = fwd, rev = rev))

  NULL
}

# Run a dada2 subprocess on the given fixture, capture results to an RDS,
# return the parsed list. The subprocess receives a small R script via
# stdin so SpeedDada doesn't have to ship per-test driver scripts.
run_dada2_subprocess <- function(fwd, rev, ref_fasta = NULL,
                                 species_fasta = NULL,
                                 trunc_len = c(150L, 150L),
                                 max_ee = c(2, 2),
                                 timeout_sec = 600) {
  out_rds <- tempfile(fileext = ".rds")
  script <- sprintf('
suppressPackageStartupMessages({
  library(dada2)
})
fwd <- "%s"; rev <- "%s"
ref <- %s; sp <- %s
trunc <- c(%d, %d); maxEE <- c(%g, %g)
filt_fwd <- tempfile(fileext = ".fastq.gz")
filt_rev <- tempfile(fileext = ".fastq.gz")
ft <- filterAndTrim(fwd, filt_fwd, rev, filt_rev,
                    truncLen = trunc, maxEE = maxEE,
                    multithread = FALSE, verbose = FALSE)
errF <- learnErrors(filt_fwd, nbases = 1e6, multithread = FALSE, verbose = FALSE)
errR <- learnErrors(filt_rev, nbases = 1e6, multithread = FALSE, verbose = FALSE)
dF <- derepFastq(filt_fwd); dR <- derepFastq(filt_rev)
aF <- dada(dF, err = errF, multithread = FALSE, verbose = FALSE)
aR <- dada(dR, err = errR, multithread = FALSE, verbose = FALSE)
m <- mergePairs(aF, dF, aR, dR, verbose = FALSE)
seqtab <- makeSequenceTable(list(s1 = m))
seqtab_nc <- removeBimeraDenovo(seqtab, method = "consensus",
                                multithread = FALSE, verbose = FALSE)
res <- list(filt = ft,
            errF = errF$err_out, errR = errR$err_out,
            seqtab = seqtab, seqtab_nc = seqtab_nc)
if (!is.null(ref)) {
  res$tax <- assignTaxonomy(colnames(seqtab_nc), ref,
                            multithread = FALSE, verbose = FALSE)
}
if (!is.null(sp)) {
  res$species <- assignSpecies(colnames(seqtab_nc), sp, verbose = FALSE)
}
saveRDS(res, "%s")
',
    fwd, rev,
    if (is.null(ref_fasta)) "NULL" else sprintf('"%s"', ref_fasta),
    if (is.null(species_fasta)) "NULL" else sprintf('"%s"', species_fasta),
    trunc_len[1], trunc_len[2], max_ee[1], max_ee[2],
    out_rds)

  result <- tryCatch(
    system2("Rscript", args = c("--vanilla", "-"),
            input = script, stdout = TRUE, stderr = TRUE,
            timeout = timeout_sec),
    error = function(e) e
  )
  if (inherits(result, "error") || !file.exists(out_rds)) {
    msg <- if (inherits(result, "error")) conditionMessage(result)
           else paste(tail(result, 20), collapse = "\n")
    stop("dada2 subprocess failed: ", msg)
  }
  readRDS(out_rds)
}

# Cell-by-cell agreement of two taxonomy character matrices.
# NA == NA counts as agreement; NA vs. value counts as disagreement.
compare_taxonomy_matrices <- function(a, b, levels = NULL) {
  if (is.null(levels)) levels <- intersect(colnames(a), colnames(b))
  shared_asvs <- intersect(rownames(a), rownames(b))
  if (length(shared_asvs) == 0L)
    return(list(agreement = 0, n = 0, level_agreement = setNames(numeric(0), levels)))
  a <- a[shared_asvs, levels, drop = FALSE]
  b <- b[shared_asvs, levels, drop = FALSE]
  matches <- (a == b) | (is.na(a) & is.na(b))
  list(
    agreement = mean(matches, na.rm = FALSE),
    n = length(matches),
    level_agreement = vapply(levels,
      function(L) mean((a[, L] == b[, L]) | (is.na(a[, L]) & is.na(b[, L])),
                       na.rm = FALSE),
      numeric(1))
  )
}

# Max abs diff in two 16xMAX_QUAL matrices (rows: 16 transitions; cols: quality).
compare_error_matrices <- function(a, b) {
  n <- min(nrow(a), nrow(b))
  k <- min(ncol(a), ncol(b))
  max(abs(a[seq_len(n), seq_len(k)] - b[seq_len(n), seq_len(k)]),
      na.rm = TRUE)
}

# Compare two sample x ASV count matrices. Returns:
#   asv_jaccard: |A∩B| / |A∪B| on the ASV column-name sets
#   max_abs_diff: per-cell max difference on the intersection of columns
compare_seqtabs <- function(a, b) {
  asvs_a <- colnames(a); asvs_b <- colnames(b)
  shared <- intersect(asvs_a, asvs_b)
  jaccard <- length(shared) / max(1L, length(union(asvs_a, asvs_b)))
  shared_rows <- intersect(rownames(a), rownames(b))
  if (length(shared) == 0L || length(shared_rows) == 0L)
    return(list(asv_jaccard = jaccard, max_abs_diff = NA_real_,
                n_shared_asvs = length(shared)))
  diff <- abs(a[shared_rows, shared, drop = FALSE] -
              b[shared_rows, shared, drop = FALSE])
  list(asv_jaccard = jaccard, max_abs_diff = max(diff),
       n_shared_asvs = length(shared))
}

# Load a baked-in dada2 snapshot if it exists. Snapshots are produced once
# (on a machine with dada2 installed) and committed so the harness can
# report meaningful diffs even when dada2 isn't available locally.
load_dada2_snapshot <- function(platform) {
  path <- file.path(testthat::test_path("_snaps"),
                    sprintf("dada2_%s.rds", platform))
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}
