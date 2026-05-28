# SpeedDada taxonomy wrappers — mirror dada2::assignTaxonomy / assignSpecies /
# addSpecies. Implements the Wang et al. 2007 naive-Bayes k-mer classifier
# with 100-rep bootstrap confidence, and (assignSpecies) exact-match species
# DB lookup with optional reverse-complement search.

#' Assign Taxonomy to ASV Sequences
#'
#' Drop-in replacement for \code{dada2::assignTaxonomy}. Implements the Wang
#' et al. (2007) naive-Bayes k-mer classifier (RDP-style) with 100-replicate
#' bootstrap confidence estimation against a reference database.
#'
#' @param seqs Character vector of ASV sequences, OR a sequence-table matrix
#'   (column names are the ASV sequences).
#' @param refFasta Path to a reference FASTA file. Headers may carry the
#'   semicolon-separated lineage directly (SILVA / GTDB style); alternatively,
#'   pass \code{lineageTsv} to specify lineages out-of-band.
#' @param minBoot Bootstrap confidence threshold (0..100) below which a level
#'   is reported as \code{NA}. Default 50 matches dada2.
#' @param tryRC If \code{TRUE}, also classify the reverse complement of each
#'   ASV and keep whichever orientation yields higher confidence.
#' @param outputBootstraps If \code{TRUE}, return a list with the taxonomy
#'   matrix \code{tax} and a per-level bootstrap matrix \code{boot}.
#' @param taxLevels Column names for the returned taxonomy matrix.
#' @param multithread Accepted for compatibility; Rust uses Rayon automatically.
#' @param verbose Logical: ignored.
#' @param lineageTsv Optional path to a TSV with
#'   \code{seq_id\\tkingdom;phylum;...} lineages (overrides FASTA headers).
#' @param k k-mer length used to build the bitset profile (default 8).
#'
#' @return Character matrix \code{[n × length(taxLevels)]} with one row per
#'   input ASV and rownames set to the ASV sequence. If
#'   \code{outputBootstraps = TRUE}, a list of two matrices (\code{tax} and
#'   \code{boot}) is returned instead.
#'
#' @export
assignTaxonomy <- function(seqs, refFasta, minBoot = 50, tryRC = FALSE,
                           outputBootstraps = FALSE,
                           taxLevels = c("Kingdom", "Phylum", "Class", "Order",
                                         "Family", "Genus", "Species"),
                           multithread = FALSE, verbose = FALSE,
                           lineageTsv = "", k = 8L) {
  .ignore_multithread(multithread)
  if (is.matrix(seqs)) seqs <- colnames(seqs)
  if (inherits(seqs, "dada") || inherits(seqs, "derep"))
    seqs <- names(.get_uniques(seqs))
  seqs <- as.character(seqs)
  if (length(seqs) == 0L)
    stop("SpeedDada: assignTaxonomy received zero sequences")

  db <- .Call("wrap__buildTaxonomyDb",
              as.character(refFasta),
              as.character(lineageTsv %||% ""),
              as.integer(k),
              as.double(0.0),      # threshold applied at classify time
              as.double(42))
  raw <- .Call("wrap__assignTaxonomy",
               seqs, db,
               as.double(minBoot),
               isTRUE(tryRC))

  cols <- list(
    Kingdom = raw$Kingdom, Phylum = raw$Phylum, Class = raw$Class,
    Order   = raw$Order,   Family = raw$Family, Genus = raw$Genus,
    Species = raw$Species
  )
  tax <- do.call(cbind, lapply(cols, function(v) {
    v[v == "" | v == "NA"] <- NA_character_
    v
  }))
  colnames(tax) <- names(cols)
  rownames(tax) <- seqs
  tax <- tax[, taxLevels, drop = FALSE]

  if (isTRUE(outputBootstraps)) {
    boot <- matrix(as.numeric(raw$bootstrap),
                   nrow = length(seqs), ncol = 7L, byrow = TRUE,
                   dimnames = list(seqs,
                                   c("Kingdom","Phylum","Class","Order",
                                     "Family","Genus","Species")))
    boot <- boot[, taxLevels, drop = FALSE]
    return(list(tax = tax, boot = boot))
  }
  tax
}

#' Build a Reusable Taxonomy Reference Database
#'
#' Pre-builds the bitset k-mer profile database from a reference FASTA so
#' that subsequent \code{\link{assignTaxonomy}} calls can reuse the index
#' (the build dominates wall time on large databases like SILVA).
#'
#' @param refFasta Path to the reference FASTA.
#' @param lineageTsv Optional path to a lineage TSV (see \code{assignTaxonomy}).
#' @param k k-mer length (default 8).
#'
#' @return An opaque \code{externalptr} consumable by \code{\link{assignTaxonomyDb}}.
#'
#' @export
buildTaxonomyDb <- function(refFasta, lineageTsv = "", k = 8L) {
  .Call("wrap__buildTaxonomyDb",
        as.character(refFasta),
        as.character(lineageTsv %||% ""),
        as.integer(k),
        as.double(0.0),
        as.double(42))
}

#' Assign Taxonomy Against a Pre-built Database
#'
#' Same as \code{\link{assignTaxonomy}} but takes a database handle from
#' \code{\link{buildTaxonomyDb}} instead of a FASTA path, avoiding the
#' index-build cost on repeated calls.
#'
#' @inheritParams assignTaxonomy
#' @param db Database handle from \code{\link{buildTaxonomyDb}}.
#'
#' @export
assignTaxonomyDb <- function(seqs, db, minBoot = 50, tryRC = FALSE,
                             outputBootstraps = FALSE,
                             taxLevels = c("Kingdom","Phylum","Class","Order",
                                           "Family","Genus","Species"),
                             multithread = FALSE) {
  .ignore_multithread(multithread)
  if (is.matrix(seqs)) seqs <- colnames(seqs)
  seqs <- as.character(seqs)
  raw <- .Call("wrap__assignTaxonomy", seqs, db,
               as.double(minBoot), isTRUE(tryRC))
  cols <- list(
    Kingdom = raw$Kingdom, Phylum = raw$Phylum, Class = raw$Class,
    Order   = raw$Order,   Family = raw$Family, Genus = raw$Genus,
    Species = raw$Species
  )
  tax <- do.call(cbind, lapply(cols, function(v) {
    v[v == "" | v == "NA"] <- NA_character_
    v
  }))
  colnames(tax) <- names(cols)
  rownames(tax) <- seqs
  tax <- tax[, taxLevels, drop = FALSE]

  if (isTRUE(outputBootstraps)) {
    boot <- matrix(as.numeric(raw$bootstrap),
                   nrow = length(seqs), ncol = 7L, byrow = TRUE,
                   dimnames = list(seqs,
                                   c("Kingdom","Phylum","Class","Order",
                                     "Family","Genus","Species")))
    return(list(tax = tax, boot = boot[, taxLevels, drop = FALSE]))
  }
  tax
}

# ── assignSpecies ────────────────────────────────────────────────────────────

#' Assign Species via Exact Sequence Match
#'
#' Drop-in replacement for \code{dada2::assignSpecies}. Matches each ASV to
#' a species-annotated reference FASTA by exact full-sequence equality
#' (and optionally a small Hamming tolerance). The reference headers must
#' carry genus and species as the first two whitespace-separated tokens of
#' the description: \code{>SeqID Genus species}.
#'
#' @param seqs Character vector of ASV sequences, or a sequence-table matrix.
#' @param refFasta Path to a species-annotated reference FASTA (e.g.,
#'   \code{silva_species_assignment_v138.fa.gz}).
#' @param allowMultiple If \code{TRUE}, ambiguous matches are joined with
#'   \code{/} on the Species column rather than dropped.
#' @param tryRC If \code{TRUE}, also try the reverse complement of each ASV.
#' @param n Maximum allowed Hamming mismatches when no exact match is found.
#'   Default 0 matches dada2.
#' @param verbose Logical: ignored.
#'
#' @return Character matrix [n x 2] with columns \code{Genus} and
#'   \code{Species}. Cells with no match are \code{NA}.
#'
#' @export
assignSpecies <- function(seqs, refFasta, allowMultiple = FALSE, tryRC = FALSE,
                          n = 0L, verbose = FALSE) {
  if (is.matrix(seqs)) seqs <- colnames(seqs)
  if (inherits(seqs, "dada") || inherits(seqs, "derep"))
    seqs <- names(.get_uniques(seqs))
  seqs <- as.character(seqs)
  if (length(seqs) == 0L)
    stop("SpeedDada: assignSpecies received zero sequences")
  raw <- .Call("wrap__assignSpecies",
               seqs, as.character(refFasta),
               isTRUE(allowMultiple), isTRUE(tryRC),
               as.integer(n))
  genus <- raw$Genus
  species <- raw$Species
  genus[genus == ""] <- NA_character_
  species[species == ""] <- NA_character_
  mat <- cbind(Genus = genus, Species = species)
  rownames(mat) <- seqs
  mat
}

# ── addSpecies ───────────────────────────────────────────────────────────────

#' Append a Species Column to a Taxonomy Table
#'
#' Drop-in replacement for \code{dada2::addSpecies}. Calls
#' \code{\link{assignSpecies}} on the ASV sequences (taken from
#' \code{rownames(taxtab)}) and overwrites or appends the \code{Genus} /
#' \code{Species} columns. The Genus from \code{assignSpecies} is only
#' written when it is consistent with the existing \code{Genus} column;
#' otherwise that row's species is dropped (matches dada2 semantics).
#'
#' @param taxtab Taxonomy character matrix with ASV sequences in
#'   \code{rownames}, as returned by \code{\link{assignTaxonomy}}.
#' @param refFasta Path to a species-annotated reference FASTA.
#' @param allowMultiple,tryRC,n,verbose See \code{\link{assignSpecies}}.
#'
#' @return The taxonomy matrix with an updated \code{Species} column (and
#'   a new \code{Genus} column if it was absent).
#'
#' @export
addSpecies <- function(taxtab, refFasta, allowMultiple = FALSE, tryRC = FALSE,
                       n = 0L, verbose = FALSE) {
  if (!is.matrix(taxtab) || is.null(rownames(taxtab)))
    stop("SpeedDada: addSpecies expects a taxonomy matrix with ASV rownames")
  seqs <- rownames(taxtab)
  sp <- assignSpecies(seqs, refFasta, allowMultiple = allowMultiple,
                      tryRC = tryRC, n = n, verbose = verbose)

  # If taxtab has a Genus column, only accept species rows whose genus matches.
  existing_genus <- if ("Genus" %in% colnames(taxtab)) taxtab[, "Genus"] else NULL
  new_species <- sp[, "Species"]
  if (!is.null(existing_genus)) {
    keep <- is.na(existing_genus) | is.na(sp[, "Genus"]) |
            existing_genus == sp[, "Genus"]
    new_species[!keep] <- NA_character_
  }

  out <- taxtab
  if ("Species" %in% colnames(out)) {
    out[, "Species"] <- new_species
  } else {
    out <- cbind(out, Species = new_species)
  }
  out
}
