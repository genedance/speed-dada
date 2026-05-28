# SpeedDada utility wrappers — mirror dada2 utility functions.

# ── plotQualityProfile ────────────────────────────────────────────────────────

#' Plot Per-Cycle Quality Profile of FASTQ File(s)
#'
#' Drop-in replacement for \code{dada2::plotQualityProfile}. Computes per-cycle
#' mean and quartile Phred scores for one or more FASTQ files (using the
#' fast Rust streaming implementation) and renders them as a ggplot2 panel.
#'
#' @param fl Character vector of FASTQ paths.
#' @param n Number of reads per file to sample (default 5e5).
#' @param aggregate If \code{TRUE}, pool reads across all files into a single
#'   profile. Currently a no-op (per-file profiles are returned and faceted).
#'
#' @return A \code{ggplot} object (or invisibly a data frame if \pkg{ggplot2}
#'   is not installed).
#'
#' @export
plotQualityProfile <- function(fl, n = 5e5, aggregate = FALSE) {
  fl <- as.character(fl)
  if (length(fl) == 0L) stop("SpeedDada: plotQualityProfile requires >= 1 file")

  df <- do.call(rbind, lapply(fl, function(path) {
    raw <- .Call("wrap__qualityProfile", as.character(path), as.integer(n))
    data.frame(
      file = basename(path),
      position = as.integer(raw$position),
      mean = as.numeric(raw$mean),
      q25 = as.numeric(raw$q25),
      q50 = as.numeric(raw$q50),
      q75 = as.numeric(raw$q75),
      count = as.integer(raw$count),
      stringsAsFactors = FALSE
    )
  }))

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("SpeedDada: install 'ggplot2' for plotting; returning data frame.")
    return(invisible(df))
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data$position)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q25, ymax = .data$q75),
                         alpha = 0.2, fill = "steelblue") +
    ggplot2::geom_line(ggplot2::aes(y = .data$q50), colour = "steelblue") +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean), colour = "black",
                       linetype = "dashed") +
    ggplot2::facet_wrap(~ .data$file) +
    ggplot2::labs(x = "Cycle", y = "Quality score",
                  title = "Per-cycle quality profile") +
    ggplot2::theme_minimal()
}


# ── rc — reverse complement ───────────────────────────────────────────────────

#' Reverse Complement DNA Sequences
#'
#' Drop-in replacement for \code{dada2::rc}. Accepts a character vector of
#' DNA sequences and returns their reverse complement. Non-ACGT bases are
#' mapped to \code{N}.
#'
#' @param seq Character vector of DNA sequences.
#'
#' @return Character vector of reverse-complemented sequences.
#'
#' @examples
#' rc(c("ACGT", "AAATTT"))
#'
#' @export
rc <- function(seq) {
  if (length(seq) == 0L) return(character(0L))
  .Call("wrap__rc", as.character(seq))
}

# ── getSequences ──────────────────────────────────────────────────────────────

#' Extract Sequences from a dada2 Object
#'
#' Drop-in replacement for \code{dada2::getSequences}. Returns the sequence
#' strings underlying a \code{derep}, \code{dada}, \code{mergePairs} data
#' frame, named integer vector, or sequence table matrix.
#'
#' @param object A \code{derep}/\code{dada} object, \code{mergePairs} data
#'   frame, named integer vector, or sequence table matrix.
#' @param collapse Accepted for compatibility; ignored.
#' @param silence Accepted for compatibility; ignored.
#'
#' @return Character vector of sequences.
#'
#' @export
getSequences <- function(object, collapse = TRUE, silence = FALSE) {
  if (is.matrix(object)) return(colnames(object))
  if (is.character(object)) return(object)
  uniq <- .get_uniques(object)
  names(uniq)
}

# ── getUniques ────────────────────────────────────────────────────────────────

#' Extract a Named Count Vector from a dada2 Object
#'
#' Drop-in replacement for \code{dada2::getUniques}. Returns the named
#' integer vector (sequence -> count) underlying a \code{derep}, \code{dada},
#' or \code{mergePairs} result.
#'
#' @param object A \code{derep}/\code{dada} object, \code{mergePairs} data
#'   frame, or named integer vector.
#' @param collapse Accepted for compatibility; ignored.
#' @param silence Accepted for compatibility; ignored.
#'
#' @return Named integer vector.
#'
#' @export
getUniques <- function(object, collapse = TRUE, silence = FALSE) {
  .get_uniques(object)
}

# ── uniquesToFasta ────────────────────────────────────────────────────────────

#' Write a Uniques Object to a FASTA File
#'
#' Drop-in replacement for \code{dada2::uniquesToFasta}.
#'
#' @param unqs A \code{derep}/\code{dada} object or named integer vector.
#' @param fout Output FASTA file path.
#' @param ids Optional character vector of FASTA IDs (default: \code{seqN}).
#' @param mode File mode: \code{"w"} (default) or \code{"a"} for append.
#' @param width Line width for sequence wrapping (default 20000 = no wrap).
#'
#' @return Invisibly returns \code{fout}.
#'
#' @export
uniquesToFasta <- function(unqs, fout, ids = NULL, mode = "w", width = 20000L) {
  uniq <- .get_uniques(unqs)
  seqs <- names(uniq)
  n <- length(seqs)
  if (is.null(ids)) ids <- sprintf("seq%d", seq_len(n))
  if (length(ids) != n) stop("SpeedDada: length(ids) must equal length(unqs)")

  con <- file(fout, open = mode)
  on.exit(close(con), add = TRUE)
  width <- as.integer(width)
  for (i in seq_len(n)) {
    writeLines(paste0(">", ids[i], ";size=", as.integer(uniq[i]), ";"), con)
    s <- seqs[i]
    if (nchar(s) <= width) {
      writeLines(s, con)
    } else {
      starts <- seq(1L, nchar(s), by = width)
      writeLines(substring(s, starts, starts + width - 1L), con)
    }
  }
  invisible(fout)
}

# ── mergeSequenceTables ───────────────────────────────────────────────────────

#' Merge Two or More Sequence Tables
#'
#' Drop-in replacement for \code{dada2::mergeSequenceTables}. Combines
#' sequence-table matrices over the union of their ASV columns. Samples are
#' row-stacked; duplicate sample names trigger an error.
#'
#' @param table1 First sequence-table matrix.
#' @param table2 Second sequence-table matrix.
#' @param ... Additional sequence-table matrices.
#' @param repeats How to resolve duplicate ASV columns within a single input
#'   table: \code{"error"} (default), \code{"sum"}, or \code{"first"}.
#' @param orderBy Column order in the result: \code{"abundance"} (default)
#'   sorts by total count, \code{NULL} keeps insertion order.
#' @param tryRC If \code{TRUE}, also match each ASV to the reverse complement
#'   of every other table's ASVs before adding it as a new column.
#'
#' @return Combined sample x ASV integer matrix.
#'
#' @export
mergeSequenceTables <- function(table1, table2, ..., repeats = "error",
                                orderBy = "abundance", tryRC = FALSE) {
  tables <- c(list(table1, table2), list(...))
  for (i in seq_along(tables)) {
    if (!is.matrix(tables[[i]]))
      stop("SpeedDada: each input must be a sequence-table matrix")
  }
  # Collapse duplicate columns within each input per `repeats` policy.
  tables <- lapply(tables, function(m) {
    dups <- duplicated(colnames(m))
    if (!any(dups)) return(m)
    if (identical(repeats, "error"))
      stop("SpeedDada: duplicate ASV column in input; set repeats='sum' to combine")
    if (identical(repeats, "first"))
      return(m[, !dups, drop = FALSE])
    if (identical(repeats, "sum")) {
      uniq <- unique(colnames(m))
      out <- matrix(0L, nrow = nrow(m), ncol = length(uniq),
                    dimnames = list(rownames(m), uniq))
      for (j in seq_along(uniq)) {
        cols <- which(colnames(m) == uniq[j])
        out[, j] <- as.integer(rowSums(m[, cols, drop = FALSE]))
      }
      return(out)
    }
    stop("SpeedDada: unknown repeats option: ", repeats)
  })

  # Optionally fold reverse-complement matches across tables.
  if (isTRUE(tryRC)) {
    rc_lookup <- rc(colnames(tables[[1]]))
    names(rc_lookup) <- colnames(tables[[1]])
    for (i in 2:length(tables)) {
      cn <- colnames(tables[[i]])
      hits <- match(cn, rc_lookup)
      keep <- !is.na(hits)
      if (any(keep)) {
        # Map this table's RC-matching columns to table 1's column names.
        cn[keep] <- names(rc_lookup)[hits[keep]]
        colnames(tables[[i]]) <- cn
      }
    }
  }

  # Sample names must be unique across tables.
  all_samples <- unlist(lapply(tables, rownames))
  if (any(duplicated(all_samples)))
    stop("SpeedDada: duplicate sample names across tables")

  all_asvs <- unique(unlist(lapply(tables, colnames)))
  out <- matrix(0L, nrow = length(all_samples), ncol = length(all_asvs),
                dimnames = list(all_samples, all_asvs))
  row <- 1L
  for (m in tables) {
    nr <- nrow(m)
    out[row:(row + nr - 1L), colnames(m)] <- as.integer(m)
    row <- row + nr
  }

  if (identical(orderBy, "abundance")) {
    ord <- order(colSums(out), decreasing = TRUE)
    out <- out[, ord, drop = FALSE]
  }
  out
}
