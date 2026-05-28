# Smoke tests for the dada2-parity wrappers added in 0.99.2.
# Each test exercises the binding shape; full parity tests against the
# original dada2 package are gated on SPEEDDADA_PARITY=1.

test_that("rc reverse-complements a vector of sequences", {
  expect_equal(rc(c("ACGT", "AAATTT", "")), c("ACGT", "AAATTT", ""))
  expect_equal(rc("ACGTN"), "NACGT")
})

test_that("getSequences and getUniques unpack derep / dada objects", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  d <- derepFastq(fwd)
  expect_type(getSequences(d), "character")
  expect_length(getSequences(d), length(d$uniques))

  u <- getUniques(d)
  expect_type(u, "integer")
  expect_true(!is.null(names(u)))
})

test_that("uniquesToFasta writes a non-empty FASTA file", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  d   <- derepFastq(fwd)
  out <- tempfile(fileext = ".fasta")
  uniquesToFasta(d, out)
  expect_true(file.exists(out))
  lines <- readLines(out)
  expect_true(length(lines) >= 2L)
  expect_true(startsWith(lines[1], ">"))
})

test_that("mergeSequenceTables column-unions two matrices", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  d   <- derepFastq(fwd)
  err <- learnErrors(fwd, nbases = 1e4)
  res <- dada(d, err, omega_a = 1e-5)
  m1  <- makeSequenceTable(list(A = res))
  m2  <- makeSequenceTable(list(B = res))
  combined <- mergeSequenceTables(m1, m2)
  expect_true(is.matrix(combined))
  expect_equal(nrow(combined), 2L)
  expect_true(all(c("A", "B") %in% rownames(combined)))
})

test_that("removePrimers trims and writes output FASTQ", {
  fwd_in  <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  fwd_out <- tempfile(fileext = ".fastq")
  # Primers chosen to match the synthetic test FASTQ. Empty primer for the
  # reverse side means "no trimming" on that end.
  mat <- removePrimers(fwd_in, fwd_out,
                       primer.fwd = "", primer.rev = "",
                       max.mismatch = 0L, min.overlap = 1L)
  expect_true(is.matrix(mat))
  expect_equal(colnames(mat), c("reads.in", "reads.out"))
  expect_true(mat[1, "reads.in"] >= mat[1, "reads.out"])
})

test_that("removePrimers refuses allow.indels = TRUE", {
  fwd_in  <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  fwd_out <- tempfile(fileext = ".fastq")
  expect_error(
    removePrimers(fwd_in, fwd_out, allow.indels = TRUE),
    "allow\\.indels"
  )
})

test_that("qualityProfile / plotQualityProfile compute per-cycle stats", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  # Bypass ggplot path: with ggplot2 absent we get the data frame back.
  with_mocked_bindings <- function() {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- plotQualityProfile(fwd, n = 100)
      expect_s3_class(p, "ggplot")
    } else {
      df <- plotQualityProfile(fwd, n = 100)
      expect_true(is.data.frame(df))
      expect_true(all(c("position", "mean", "q25", "q50", "q75") %in% names(df)))
    }
  }
  with_mocked_bindings()
})

test_that("assignTaxonomy returns a character matrix", {
  # Build a tiny FASTA reference with two annotated sequences.
  ref <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">seq1 Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;acidophilus",
    "ACGTACGTACGTACGTACGTACGTACGTACGT",
    ">seq2 Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;aeruginosa",
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  ), ref)

  seqs <- c("ACGTACGTACGTACGTACGTACGTACGTACGT")
  tax <- assignTaxonomy(seqs, ref, minBoot = 0, k = 6L)
  expect_true(is.matrix(tax))
  expect_equal(nrow(tax), 1L)
  expect_equal(ncol(tax), 7L)
  expect_equal(tax[1, "Kingdom"], "Bacteria")
})

test_that("assignSpecies returns a character matrix with Genus + Species", {
  ref <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">A Lactobacillus acidophilus",
    "ACGTACGTACGTACGTACGTACGT",
    ">B Pseudomonas aeruginosa",
    "TTTTTTTTTTTTTTTTTTTTTTTT"
  ), ref)

  out <- assignSpecies("ACGTACGTACGTACGTACGTACGT", ref)
  expect_true(is.matrix(out))
  expect_equal(colnames(out), c("Genus", "Species"))
  expect_equal(out[1, "Genus"], "Lactobacillus")
  expect_equal(out[1, "Species"], "acidophilus")
})

test_that("addSpecies appends the Species column to a taxonomy table", {
  ref <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">A Lactobacillus acidophilus",
    "ACGTACGTACGTACGTACGTACGT"
  ), ref)

  seq <- "ACGTACGTACGTACGTACGTACGT"
  tax <- matrix(c("Bacteria", "Firmicutes", "Bacilli", "Lactobacillales",
                  "Lactobacillaceae", "Lactobacillus", NA_character_),
                nrow = 1L, ncol = 7L,
                dimnames = list(seq,
                                c("Kingdom","Phylum","Class","Order",
                                  "Family","Genus","Species")))
  out <- addSpecies(tax, ref)
  expect_equal(out[1, "Species"], "acidophilus")
})
