# SpeedDada — drop-in R wrappers over Rust extendr bindings.
# Each function mirrors the R dada2 API so existing scripts work unchanged.

#' @useDynLib SpeedDada, .registration = TRUE
NULL

# ── Internal helpers ──────────────────────────────────────────────────────────

# Accept and ignore multithread: Rust uses Rayon automatically.
.ignore_multithread <- function(multithread) invisible(NULL)

# Extract a named integer vector (seq → count) from dada/derep/data.frame.
.get_uniques <- function(x) {
  if (inherits(x, "dada"))  return(x$denoised)
  if (inherits(x, "derep")) return(x$uniques)
  if (is.data.frame(x) && all(c("sequence", "abundance") %in% names(x))) {
    v <- as.integer(x$abundance[x$accept])
    names(v) <- x$sequence[x$accept]
    return(v)
  }
  if (is.integer(x) && !is.null(names(x))) return(x)
  stop("SpeedDada: unrecognised sample type in .get_uniques()")
}

# Resolve a bundled FASTQ fixture for examples / tests.
.fixture <- function(name) {
  system.file("extdata", name, package = "SpeedDada", mustWork = TRUE)
}

# ── 1. filterAndTrim ─────────────────────────────────────────────────────────

#' Filter and Trim FASTQ Reads
#'
#' Drop-in replacement for \code{dada2::filterAndTrim}. Filters and trims
#' single- or paired-end FASTQ files using the high-performance Rust core.
#'
#' @param fwd Character vector of forward-read input paths.
#' @param filt Character vector of forward-read output paths.
#' @param rev Character vector of reverse-read input paths, or \code{NULL}
#'   for single-end.
#' @param filt.rev Character vector of reverse-read output paths, or \code{NULL}.
#' @param truncLen Integer of length 1 or 2: truncation length per direction.
#' @param trimLeft Integer of length 1 or 2: bases to trim from the 5' end.
#' @param maxEE Numeric of length 1 or 2: maximum expected errors.
#' @param truncQ Integer: truncate at first base below this Phred score.
#' @param minLen Integer: discard reads shorter than this after truncation.
#' @param compress Accepted for compatibility; output is plain FASTQ.
#' @param multithread Accepted for compatibility; Rust uses Rayon automatically.
#' @param ... Extra arguments ignored for compatibility.
#'
#' @return Integer matrix of dimensions [files x 2] with column names
#'   \code{c("reads.in", "reads.out")}.
#'
#' @examples
#' fwd_in  <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
#' fwd_out <- tempfile(fileext = ".fastq")
#' filterAndTrim(fwd_in, fwd_out, truncLen = 30, maxEE = 5, minLen = 10)
#'
#' @export
filterAndTrim <- function(fwd, filt, rev = NULL, filt.rev = NULL,
                          truncLen = 0L, trimLeft = 0L, maxEE = Inf,
                          truncQ = 2L, minLen = 20L,
                          compress = TRUE, multithread = FALSE, ...) {
  .ignore_multithread(multithread)

  # Vectorise per-direction parameters
  to2 <- function(x) { x <- rep_len(as.numeric(x), 2L); x[!is.finite(x)] <- 1e9; x }
  truncLen <- to2(truncLen)
  trimLeft <- to2(trimLeft)
  maxEE    <- to2(maxEE)

  paired <- !is.null(rev)
  raw <- .Call("wrap__filterAndTrim",
    as.character(fwd),
    as.character(filt),
    if (paired) as.character(rev)      else NULL,
    if (paired) as.character(filt.rev) else NULL,
    as.integer(truncLen[1L]),
    as.integer(truncLen[2L]),
    as.integer(trimLeft[1L]),
    as.integer(trimLeft[2L]),
    as.double(maxEE[1L]),
    as.double(maxEE[2L]),
    as.integer(truncQ),
    as.integer(minLen))

  # Build [n × 2] integer matrix matching dada2's return type.
  mat <- matrix(
    c(as.integer(raw$reads_in), as.integer(raw$reads_out)),
    nrow = length(raw$reads_in),
    ncol = 2L,
    dimnames = list(raw$rownames, c("reads.in", "reads.out"))
  )
  invisible(mat)
}

# ── 1b. removePrimers ────────────────────────────────────────────────────────

#' Remove PCR Primers from FASTQ Reads
#'
#' Drop-in replacement for \code{dada2::removePrimers}. Locates the forward
#' and reverse primers in each read and writes the trimmed insert to
#' \code{fout}. Reads where either primer cannot be matched within
#' \code{max.mismatch} are discarded.
#'
#' @param fn Character vector of input FASTQ paths.
#' @param fout Character vector of output FASTQ paths.
#' @param primer.fwd Forward primer sequence (5'->3'). Empty string skips.
#' @param primer.rev Reverse primer sequence (5'->3' on the same strand
#'   as the read; SpeedDada matches it near the 3' end). Empty string skips.
#' @param max.mismatch Maximum mismatches when locating each primer.
#' @param min.overlap Minimum primer bases required to call a match.
#' @param orient Accepted for compatibility; SpeedDada always searches the
#'   forward orientation only. Reverse-orientation reads should be reverse-
#'   complemented before calling this function.
#' @param allow.indels Not supported; raises an error if \code{TRUE}.
#' @param compress Accepted for compatibility; output is plain FASTQ.
#' @param verbose Logical: ignored.
#' @param ... Extra arguments ignored.
#'
#' @return Integer matrix [files x 2] with column names
#'   \code{c("reads.in", "reads.out")}.
#'
#' @export
removePrimers <- function(fn, fout, primer.fwd = "", primer.rev = "",
                          max.mismatch = 2L, min.overlap = 4L,
                          orient = TRUE, allow.indels = FALSE,
                          compress = TRUE, verbose = FALSE, ...) {
  if (isTRUE(allow.indels))
    stop("SpeedDada: removePrimers(allow.indels = TRUE) is not supported; ",
         "the current primer matcher uses Hamming distance only")
  if (length(fn) != length(fout))
    stop("SpeedDada: length(fn) must equal length(fout)")

  raw <- .Call("wrap__removePrimers",
    as.character(fn),
    as.character(fout),
    as.character(primer.fwd %||% ""),
    as.character(primer.rev %||% ""),
    as.integer(max.mismatch),
    as.integer(min.overlap),
    isTRUE(orient))

  mat <- matrix(
    c(as.integer(raw$reads_in), as.integer(raw$reads_out)),
    nrow = length(raw$reads_in),
    ncol = 2L,
    dimnames = list(raw$rownames, c("reads.in", "reads.out"))
  )
  invisible(mat)
}

# Local NULL-coalesce so users don't need to depend on rlang.
`%||%` <- function(a, b) if (is.null(a)) b else a

# ── 2. learnErrors ───────────────────────────────────────────────────────────

#' Learn Error Rates from FASTQ Files
#'
#' Drop-in replacement for \code{dada2::learnErrors}. Estimates a substitution
#' error model from FASTQ data using the Rust EM implementation.
#'
#' The default \code{errFun = "auto"} sniffs the quality profile and picks
#' \code{loess} for full-range Illumina (MiSeq / HiSeq full quality),
#' \code{binned} for binned-quality platforms (NovaSeq, NextSeq, MGI
#' DNBSEQ), or \code{pacbio} for near-Q40 long reads (PacBio CCS / HiFi).
#' Oxford Nanopore data is detected and warned about — the substitution-
#' dominant model used here will not track ONT's indel-heavy errors well;
#' proper ONT support requires indel-aware alignment.
#'
#' @param fls Character vector of FASTQ paths, or a directory containing them.
#' @param nbases Total bases to use for learning (converted to approx. read
#'   count internally).
#' @param errFun One of \code{"auto"} (default), \code{"loess"},
#'   \code{"binned"}, \code{"pacbio"}, or \code{"logistic"}. Maps to dada2's
#'   \code{loessErrfun} / \code{makeBinnedQualErrfun} / \code{PacBioErrfun}.
#' @param multithread Accepted for compatibility; Rust uses Rayon automatically.
#' @param verbose Logical: if \code{TRUE}, print the detected platform / errFun.
#' @param ... Extra arguments ignored for compatibility.
#'
#' @return Opaque error-model handle (an R \code{externalptr}) consumed by
#'   \code{\link{dada}}.
#'
#' @examples
#' fastq <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
#' err <- learnErrors(fastq, nbases = 1e4)
#' err
#'
#' @export
learnErrors <- function(fls, nbases = 1e8, errFun = "auto",
                        multithread = FALSE, verbose = FALSE, ...) {
  .ignore_multithread(multithread)
  if (length(fls) == 1L && dir.exists(fls)) {
    fls <- list.files(fls, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE)
  }
  errFun <- tolower(as.character(errFun))[1L]
  if (identical(errFun, "auto")) {
    info <- .Call("wrap__detectErrFun", as.character(fls), as.double(nbases))
    if (isTRUE(verbose)) {
      message(sprintf(
        "SpeedDada::learnErrors auto-detect: platform=%s, distinct_q=%d, mean_q=%.1f, max_q=%d, mean_len=%.0f",
        info$kind, as.integer(info$distinct_q),
        as.numeric(info$mean_q), as.integer(info$max_q),
        as.numeric(info$mean_len)))
    }
    if (as.numeric(info$mean_q) < 20 && as.numeric(info$mean_len) > 500) {
      warning(
        "SpeedDada: input looks like Oxford Nanopore (mean Phred ",
        sprintf("%.1f", as.numeric(info$mean_q)),
        ", mean read length ", sprintf("%.0f", as.numeric(info$mean_len)),
        " bp). SpeedDada's substitution-dominant error model and single-",
        "indel gap correction will not track ONT data faithfully — ASVs ",
        "produced from this run will be inaccurate. Proper ONT support ",
        "requires banded indel-aware alignment (not yet implemented).")
    }
  }
  .Call("wrap__learnErrors",
        as.character(fls), as.double(nbases), as.character(errFun))
}

# ── 3. derepFastq ────────────────────────────────────────────────────────────

#' Dereplicate FASTQ File(s)
#'
#' Drop-in replacement for \code{dada2::derepFastq}. Streams the FASTQ file
#' through the Rust dereplicator, returning a \code{"derep"} object.
#'
#' @param fls Character vector of FASTQ paths.
#' @param verbose Logical: ignored (accepted for compatibility).
#' @param ... Extra arguments ignored.
#'
#' @return For one file: a \code{"derep"} list with \code{$uniques} (named
#'   integer vector of sequence → count). For multiple files: a named list
#'   of such objects. The internal \code{.rust_ptr} field carries the full
#'   per-position quality scores across the FFI boundary so that downstream
#'   \code{\link{dada}} calls use real qualities, not a placeholder.
#'
#' @examples
#' fastq <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
#' d <- derepFastq(fastq)
#' head(d$uniques)
#'
#' @export
derepFastq <- function(fls, verbose = FALSE, ...) {
  make_derep <- function(path) {
    raw <- .Call("wrap__derepFastq", as.character(path))
    uniq <- as.integer(raw$counts)
    names(uniq) <- raw$seqs
    # `.rust_ptr` carries per-position quality across the FFI for `dada()`.
    # The user-visible `$uniques` (seq → count) stays for back-compat.
    structure(list(uniques = uniq, quals = NULL, map = NULL,
                   .rust_ptr = raw$ptr),
              class = "derep")
  }
  if (length(fls) == 1L) {
    make_derep(fls)
  } else {
    nms <- if (!is.null(names(fls))) names(fls) else basename(fls)
    setNames(lapply(fls, make_derep), nms)
  }
}

# ── 4. dada ──────────────────────────────────────────────────────────────────

#' Denoise Amplicon Reads
#'
#' Drop-in replacement for \code{dada2::dada}. Runs the DADA partition + EM
#' refinement algorithm on dereplicated samples.
#'
#' @param derep A \code{"derep"} object from \code{\link{derepFastq}}, or a
#'   list of them for multi-sample denoising.
#' @param err Error model from \code{\link{learnErrors}}.
#' @param selfConsist Logical: not implemented; emits a warning if \code{TRUE}.
#' @param pool Logical or \code{"pseudo"}: \code{TRUE} enables cross-sample
#'   pooling; \code{"pseudo"} uses Callahan's two-pass pseudo-pooling.
#' @param omega_a Numeric: DADA significance threshold (default \code{1e-40}).
#' @param multithread Accepted for compatibility; Rust uses Rayon automatically.
#' @param verbose Logical: ignored.
#' @param ... Extra arguments ignored.
#'
#' @return A \code{"dada"} object with \code{$denoised} (named integer vector
#'   of sequence → count), or a list of such objects if \code{derep} is a list.
#'
#' @examples
#' fastq <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
#' d <- derepFastq(fastq)
#' err <- learnErrors(fastq, nbases = 1e4)
#' res <- dada(d, err, omega_a = 1e-5)
#' head(res$denoised)
#'
#' @export
dada <- function(derep, err, selfConsist = FALSE, pool = FALSE,
                 omega_a = 1e-40, multithread = FALSE, verbose = TRUE, ...) {
  .ignore_multithread(multithread)
  if (!isFALSE(selfConsist))
    warning("SpeedDada: selfConsist is not implemented and was ignored")

  pool_flag   <- isTRUE(pool)
  pseudo_flag <- identical(pool, "pseudo")

  # ── Cross-sample multi-sample path ────────────────────────────────────────
  if (is.list(derep) && !inherits(derep, "derep") && length(derep) >= 2) {
    ptrs <- vector("list", length(derep))
    for (i in seq_along(derep)) {
      d <- derep[[i]]
      if (!inherits(d, "derep") || is.null(d$.rust_ptr))
        stop("SpeedDada: each element of derep must be a 'derep' object with a .rust_ptr")
      ptrs[[i]] <- d$.rust_ptr
    }
    entry <- if (pseudo_flag)    "wrap__dada_pseudo"
             else if (pool_flag) "wrap__dada_pooled"
             else                "wrap__dada_many"
    raw <- .Call(entry, ptrs, err, as.double(omega_a))

    out <- vector("list", length(derep))
    names(out) <- if (!is.null(names(derep))) names(derep)
                  else paste0("Sample", seq_along(derep))
    for (i in seq_along(derep)) {
      mask <- raw$sample_idx == (i - 1L)
      denoised <- as.integer(raw$counts[mask])
      names(denoised) <- raw$seqs[mask]
      out[[i]] <- structure(list(denoised = denoised), class = "dada")
    }
    return(out)
  }

  # ── Single-sample path ────────────────────────────────────────────────────
  run_one <- function(d) {
    if (!inherits(d, "derep") || is.null(d$.rust_ptr))
      stop("SpeedDada: derep must be a 'derep' object with a .rust_ptr")
    raw <- .Call("wrap__dada",
      d$.rust_ptr,
      err,
      as.double(omega_a),
      pool_flag)
    denoised <- as.integer(raw$counts)
    names(denoised) <- raw$seqs
    structure(list(denoised = denoised), class = "dada")
  }

  if (inherits(derep, "derep")) run_one(derep)
  else if (is.list(derep))      lapply(derep, run_one)
  else stop("SpeedDada: derep must be a 'derep' object or list thereof")
}

# ── 5. mergePairs ────────────────────────────────────────────────────────────

#' Merge Paired-End Reads
#'
#' Drop-in replacement for \code{dada2::mergePairs}.
#'
#' @param dadaF \code{"dada"} object for forward reads.
#' @param derepF \code{"derep"} object for forward reads (accepted for compatibility).
#' @param dadaR \code{"dada"} object for reverse reads.
#' @param derepR \code{"derep"} object for reverse reads (accepted for compatibility).
#' @param minOverlap Integer: minimum overlap length (default 12).
#' @param maxMismatch Integer: maximum mismatches in overlap (default 0).
#' @param justConcatenate Logical: concatenate instead of merging.
#' @param verbose Logical: ignored.
#' @param ... Extra arguments ignored.
#'
#' @return Data frame with columns \code{sequence}, \code{abundance},
#'   \code{accept}, \code{nmatch}, \code{nmismatch}, \code{nindel}.
#'
#' @examples
#' fwd <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
#' rev <- system.file("extdata", "sample_R2.fastq", package = "SpeedDada")
#' dF <- derepFastq(fwd); dR <- derepFastq(rev)
#' err <- learnErrors(c(fwd, rev), nbases = 1e4)
#' aF <- dada(dF, err, omega_a = 1e-5); aR <- dada(dR, err, omega_a = 1e-5)
#' mergePairs(aF, dF, aR, dR, minOverlap = 5)
#'
#' @export
mergePairs <- function(dadaF, derepF, dadaR, derepR,
                       minOverlap = 12L, maxMismatch = 0L,
                       justConcatenate = FALSE, verbose = FALSE, ...) {
  if (!inherits(dadaF, "dada") || !inherits(dadaR, "dada"))
    stop("SpeedDada: dadaF and dadaR must be 'dada' objects from dada()")

  fwd <- dadaF$denoised
  rev <- dadaR$denoised

  raw <- .Call("wrap__mergePairs",
    as.character(names(fwd)), as.integer(fwd),
    as.character(names(rev)),  as.integer(rev),
    as.integer(minOverlap),
    as.integer(maxMismatch),
    isTRUE(justConcatenate))

  data.frame(
    sequence  = raw$sequence,
    abundance = raw$abundance,
    accept    = raw$accept,
    nmatch    = raw$nmatch,
    nmismatch = raw$nmismatch,
    nindel    = raw$nindel,
    stringsAsFactors = FALSE
  )
}

# ── 6. makeSequenceTable ─────────────────────────────────────────────────────

#' Build Sample x ASV Count Matrix
#'
#' Drop-in replacement for \code{dada2::makeSequenceTable}.
#'
#' @param samples List of \code{dada} / \code{derep} / \code{mergePairs}
#'   data frames, optionally named.
#' @param orderBy \code{"abundance"} (default) sorts columns by total count;
#'   \code{""} keeps insertion order.
#'
#' @return Integer matrix [samples x ASV sequences].
#'
#' @examples
#' fastq <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
#' d <- derepFastq(fastq)
#' err <- learnErrors(fastq, nbases = 1e4)
#' res <- dada(d, err, omega_a = 1e-5)
#' makeSequenceTable(list(s1 = res))
#'
#' @export
makeSequenceTable <- function(samples, orderBy = "abundance") {
  if (!is.list(samples) ||
        inherits(samples, "dada") ||
        inherits(samples, "derep") ||
        is.data.frame(samples))
    samples <- list(samples)

  snames <- if (!is.null(names(samples))) names(samples)
            else paste0("Sample", seq_along(samples))

  all_seqs    <- character(0L)
  all_counts  <- integer(0L)
  all_idx     <- integer(0L)

  for (i in seq_along(samples)) {
    uniq <- .get_uniques(samples[[i]])
    all_seqs   <- c(all_seqs,   names(uniq))
    all_counts <- c(all_counts, as.integer(uniq))
    all_idx    <- c(all_idx,    rep(i - 1L, length(uniq)))
  }

  raw <- .Call("wrap__makeSequenceTable",
    as.character(snames),
    as.character(all_seqs),
    as.integer(all_counts),
    as.integer(all_idx),
    identical(orderBy, "abundance"))

  matrix(
    as.integer(raw$data),
    nrow = length(raw$samples),
    ncol = length(raw$seqs),
    dimnames = list(raw$samples, raw$seqs)
  )
}

# ── 7. removeBimeraDenovo ────────────────────────────────────────────────────

#' Remove Chimeric Sequences (Bimera Detection)
#'
#' Drop-in replacement for \code{dada2::removeBimeraDenovo}.
#'
#' @param unqs Integer matrix (sequence table), named integer vector, dada
#'   object, derep object, or mergePairs data frame.
#' @param method Accepted for compatibility (the Rust core uses a single,
#'   consensus-style scoring path).
#' @param multithread Accepted for compatibility; Rust uses Rayon automatically.
#' @param verbose Logical: ignored.
#' @param ... Extra arguments ignored.
#'
#' @return Same type as input with chimeric sequences removed.
#'
#' @examples
#' fastq <- system.file("extdata", "sample_R1.fastq", package = "SpeedDada")
#' d <- derepFastq(fastq)
#' err <- learnErrors(fastq, nbases = 1e4)
#' res <- dada(d, err, omega_a = 1e-5)
#' mat <- makeSequenceTable(list(s1 = res))
#' removeBimeraDenovo(mat)
#'
#' @export
removeBimeraDenovo <- function(unqs, method = "consensus",
                               multithread = FALSE, verbose = FALSE, ...) {
  .ignore_multithread(multithread)

  if (is.matrix(unqs)) {
    seqs   <- colnames(unqs)
    counts <- as.integer(colSums(unqs))
    mask   <- .Call("wrap__removeBimeraDenovo",
                    as.character(seqs), as.integer(counts))
    return(unqs[, mask == 1L, drop = FALSE])
  }

  uniq <- .get_uniques(unqs)
  mask <- .Call("wrap__removeBimeraDenovo",
                as.character(names(uniq)), as.integer(uniq))
  uniq[mask == 1L]
}
