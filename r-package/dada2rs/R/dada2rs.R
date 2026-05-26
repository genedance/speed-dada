# dada2rs — drop-in R wrappers over Rust extendr bindings.
# Each function mirrors the R dada2 API so existing scripts work unchanged.

#' @useDynLib dada2rs, .registration = TRUE
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
  stop("dada2rs: unrecognised sample type in .get_uniques()")
}

# ── 1. filterAndTrim ─────────────────────────────────────────────────────────

#' Filter and trim FASTQ reads (drop-in for dada2::filterAndTrim)
#'
#' @param fwd    character vector of forward-read input paths.
#' @param filt   character vector of forward-read output paths.
#' @param rev    character vector of reverse-read input paths, or NULL.
#' @param filt.rev character vector of reverse-read output paths, or NULL.
#' @param truncLen integer(1 or 2): truncation length per direction.
#' @param trimLeft integer(1 or 2): bases to trim from the 5' end.
#' @param maxEE  numeric(1 or 2): maximum expected errors.
#' @param truncQ integer: truncate at first base below this Phred score.
#' @param minLen integer: discard reads shorter than this after truncation.
#' @param compress logical: accepted for compatibility, output is plain FASTQ.
#' @param multithread accepted for compatibility, Rust uses Rayon automatically.
#' @param ...    extra arguments ignored for compatibility.
#' @return integer matrix [files × c("reads.in", "reads.out")].
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

# ── 2. learnErrors ───────────────────────────────────────────────────────────

#' Learn error rates from FASTQ files (drop-in for dada2::learnErrors)
#'
#' @param fls character vector of FASTQ paths, or a directory containing them.
#' @param nbases total bases to use for learning (converted to approx read count).
#' @param multithread accepted for compatibility.
#' @param ...   extra arguments ignored.
#' @return opaque error model handle (R externalptr).
learnErrors <- function(fls, nbases = 1e8, multithread = FALSE, ...) {
  .ignore_multithread(multithread)
  if (length(fls) == 1L && dir.exists(fls)) {
    fls <- list.files(fls, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE)
  }
  .Call("wrap__learnErrors", as.character(fls), as.double(nbases))
}

# ── 3. derepFastq ────────────────────────────────────────────────────────────

#' Dereplicate FASTQ file(s) (drop-in for dada2::derepFastq)
#'
#' @param fls character vector of FASTQ paths.
#' @param verbose logical: ignored.
#' @param ...  extra arguments ignored.
#' @return For one file: a "derep" list with \code{$uniques} (named integer vector).
#'         For multiple files: a named list of such objects.
derepFastq <- function(fls, verbose = FALSE, ...) {
  make_derep <- function(path) {
    raw <- .Call("wrap__derepFastq", as.character(path))
    uniq <- as.integer(raw$counts)
    names(uniq) <- raw$seqs
    structure(list(uniques = uniq, quals = NULL, map = NULL), class = "derep")
  }
  if (length(fls) == 1L) {
    make_derep(fls)
  } else {
    nms <- if (!is.null(names(fls))) names(fls) else basename(fls)
    setNames(lapply(fls, make_derep), nms)
  }
}

# ── 4. dada ──────────────────────────────────────────────────────────────────

#' Denoise amplicon reads (drop-in for dada2::dada)
#'
#' @param derep  a "derep" object from \code{derepFastq}, or a list of them.
#' @param err    error model from \code{learnErrors}.
#' @param selfConsist logical: not implemented; emits a warning if TRUE.
#' @param pool   logical or "pseudo": TRUE enables cross-sample pooling.
#' @param omega_a numeric: DADA significance threshold (default 1e-40).
#' @param multithread accepted for compatibility.
#' @param verbose logical: ignored.
#' @param ...    extra arguments ignored.
#' @return "dada" object with \code{$denoised} (named integer vector seq→count),
#'         or a list of such objects if \code{derep} is a list.
dada <- function(derep, err, selfConsist = FALSE, pool = FALSE,
                 omega_a = 1e-40, multithread = FALSE, verbose = TRUE, ...) {
  .ignore_multithread(multithread)
  if (!isFALSE(selfConsist))
    warning("dada2rs: selfConsist is not implemented and was ignored")

  pool_flag   <- isTRUE(pool)
  pseudo_flag <- identical(pool, "pseudo")

  # ── Cross-sample multi-sample path ────────────────────────────────────────
  # Triggered for any list-of-derep input with >=2 samples. Dispatches to one
  # of three Rust orchestrators (all parallelised across Rayon):
  #   pool="pseudo" → wrap__dada_pseudo  (two-pass with priors)
  #   pool=TRUE     → wrap__dada_pooled  (true cross-sample pool)
  #   pool=FALSE    → wrap__dada_many    (independent per-sample, parallel)
  # All three flatten the per-sample dereps into parallel vectors and re-split
  # ASVs back to a per-sample list of "dada" objects.
  if (is.list(derep) && !inherits(derep, "derep") && length(derep) >= 2) {
    sample_idx <- integer(0)
    seqs       <- character(0)
    counts     <- integer(0)
    for (i in seq_along(derep)) {
      d <- derep[[i]]
      if (!inherits(d, "derep"))
        stop("dada2rs: each element of derep must be a 'derep' object")
      uniq <- d$uniques
      sample_idx <- c(sample_idx, rep.int(i - 1L, length(uniq)))
      seqs       <- c(seqs,       names(uniq))
      counts     <- c(counts,     as.integer(uniq))
    }
    entry <- if (pseudo_flag)    "wrap__dada_pseudo"
             else if (pool_flag) "wrap__dada_pooled"
             else                "wrap__dada_many"
    raw <- .Call(entry,
      as.integer(sample_idx),
      as.character(seqs),
      as.integer(counts),
      as.integer(length(derep)),
      err,
      as.double(omega_a))

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

  # ── Per-sample path ────────────────────────────────────────────────────────
  run_one <- function(d) {
    if (!inherits(d, "derep")) stop("dada2rs: each element of derep must be a 'derep' object")
    uniq <- d$uniques
    raw <- .Call("wrap__dada",
      as.character(names(uniq)),
      as.integer(uniq),
      err,
      as.double(omega_a),
      pool_flag)
    denoised <- as.integer(raw$counts)
    names(denoised) <- raw$seqs
    structure(list(denoised = denoised), class = "dada")
  }

  if (inherits(derep, "derep")) run_one(derep)
  else if (is.list(derep))      lapply(derep, run_one)
  else stop("dada2rs: derep must be a 'derep' object or list thereof")
}

# ── 5. mergePairs ────────────────────────────────────────────────────────────

#' Merge paired-end reads (drop-in for dada2::mergePairs)
#'
#' @param dadaF  "dada" object for forward reads.
#' @param derepF "derep" object for forward reads (accepted for compatibility).
#' @param dadaR  "dada" object for reverse reads.
#' @param derepR "derep" object for reverse reads (accepted for compatibility).
#' @param minOverlap integer: minimum overlap length (default 12).
#' @param maxMismatch integer: maximum mismatches in overlap (default 0).
#' @param justConcatenate logical: concatenate instead of merging.
#' @param verbose logical: ignored.
#' @param ...    extra arguments ignored.
#' @return data.frame with columns sequence, abundance, accept, nmatch, nmismatch, nindel.
mergePairs <- function(dadaF, derepF, dadaR, derepR,
                       minOverlap = 12L, maxMismatch = 0L,
                       justConcatenate = FALSE, verbose = FALSE, ...) {
  if (!inherits(dadaF, "dada") || !inherits(dadaR, "dada"))
    stop("dada2rs: dadaF and dadaR must be 'dada' objects from dada()")

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

#' Build sample × ASV count matrix (drop-in for dada2::makeSequenceTable)
#'
#' @param samples list of dada / derep / mergePairs data.frames, optionally named.
#' @param orderBy "abundance" to sort columns by total count (default), or "" for
#'        insertion order.
#' @return integer matrix [samples × ASV sequences].
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

#' Remove chimeric sequences (drop-in for dada2::removeBimeraDenovo)
#'
#' @param unqs    integer matrix (sequence table), named integer vector, dada,
#'                derep, or mergePairs data.frame.
#' @param method  accepted for compatibility (always uses "consensus").
#' @param multithread accepted for compatibility.
#' @param verbose logical: ignored.
#' @param ...     extra arguments ignored.
#' @return Same type as input with chimeric sequences removed.
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
