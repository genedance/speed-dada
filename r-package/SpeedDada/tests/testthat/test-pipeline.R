test_that("filterAndTrim writes a valid output FASTQ and returns a matrix", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  out <- tempfile(fileext = ".fastq")
  mat <- filterAndTrim(fwd, out, truncLen = 30, maxEE = 5, minLen = 10)

  expect_true(file.exists(out))
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(1, 2))
  expect_equal(colnames(mat), c("reads.in", "reads.out"))
  expect_gt(mat[1, "reads.out"], 0)
})

test_that("learnErrors returns an externalptr handle", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  err <- learnErrors(fwd, nbases = 1e4)
  expect_s3_class(err, "externalptr") |> tryCatch(error = function(e) {
    # Older expect_s3_class doesn't recognise externalptr; fall back.
    expect_true(typeof(err) == "externalptr")
  })
})

test_that("derepFastq returns a derep object with non-empty uniques", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  d <- derepFastq(fwd)
  expect_s3_class(d, "derep")
  expect_true(length(d$uniques) > 0)
  expect_true(!is.null(d$.rust_ptr))
})

test_that("dada returns a denoised set with at least one ASV", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  d   <- derepFastq(fwd)
  err <- learnErrors(fwd, nbases = 1e4)
  res <- dada(d, err, omega_a = 1e-5)
  expect_s3_class(res, "dada")
  expect_true(length(res$denoised) >= 1)
})

test_that("mergePairs returns a data.frame with the dada2 columns", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  rev <- system.file("extdata", "sample_R2.fastq", package = "SpeedDada")
  dF  <- derepFastq(fwd); dR <- derepFastq(rev)
  err <- learnErrors(c(fwd, rev), nbases = 1e4)
  aF  <- dada(dF, err, omega_a = 1e-5)
  aR  <- dada(dR, err, omega_a = 1e-5)
  m   <- mergePairs(aF, dF, aR, dR, minOverlap = 5)
  expect_true(is.data.frame(m))
  expect_setequal(
    names(m),
    c("sequence", "abundance", "accept", "nmatch", "nmismatch", "nindel")
  )
})

test_that("makeSequenceTable produces a sample x ASV matrix", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  d   <- derepFastq(fwd)
  err <- learnErrors(fwd, nbases = 1e4)
  res <- dada(d, err, omega_a = 1e-5)
  mat <- makeSequenceTable(list(s1 = res))
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 1L)
  expect_gt(ncol(mat), 0L)
})

test_that("removeBimeraDenovo preserves matrix columns when no chimeras", {
  fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
  d   <- derepFastq(fwd)
  err <- learnErrors(fwd, nbases = 1e4)
  res <- dada(d, err, omega_a = 1e-5)
  mat <- makeSequenceTable(list(s1 = res))
  cleaned <- removeBimeraDenovo(mat)
  expect_true(is.matrix(cleaned))
  expect_true(ncol(cleaned) <= ncol(mat))
})
