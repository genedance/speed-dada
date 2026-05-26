#!/usr/bin/env Rscript
# Benchmark: R dada2 (Bioconductor) on 6 paired manure samples.
# Multi-sample workflow:
#   filterAndTrim → learnErrors (fwd+rev pooled) → derepFastq (per sample)
#   → dada (per sample) → mergePairs → makeSequenceTable → removeBimeraDenovo

.libPaths(c(path.expand("~/R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(dada2)
  library(jsonlite)
})

args     <- commandArgs(trailingOnly = TRUE)
n_threads <- if (length(args) >= 1) as.integer(args[1]) else 16L
in_dir    <- if (length(args) >= 2) args[2] else "/Users/alex/Downloads/raw_data_FIELD"
out_dir   <- if (length(args) >= 3) args[3] else "/tmp/bench_field_out/r"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- sprintf("T0-Manure-rep%d", 1:6)
fwd_in  <- file.path(in_dir, paste0(samples, "_R1.fastq.gz"))
rev_in  <- file.path(in_dir, paste0(samples, "_R2.fastq.gz"))
fwd_filt <- file.path(out_dir, paste0(samples, "_R1_filt.fastq.gz"))
rev_filt <- file.path(out_dir, paste0(samples, "_R2_filt.fastq.gz"))
names(fwd_filt) <- samples
names(rev_filt) <- samples

cat(sprintf("[R dada2] n_threads=%d  samples=%d\n", n_threads, length(samples)))

# Modified loess error-rate estimator for NextSeq/NovaSeq binned Phred scores.
# Without this, learnErrors fails on binned-quality FASTQs (only ~4 distinct Q values).
# Source: https://github.com/benjjneb/dada2/issues/1307 (community-standard fix).
loessErrfun_binned <- function(trans) {
  qq <- as.numeric(colnames(trans))
  err <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs  <- trans[paste0(nti, "2", ntj), ]
        tot   <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1) / tot)
        rlogp[is.infinite(rlogp)] <- NA
        df    <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
        mod   <- suppressWarnings(loess(rlogp ~ q, df,
                                        weights = log10(pmax(tot, 1))))
        pred  <- predict(mod, qq)
        pred[!is.finite(pred)] <- 0
        err   <- rbind(err, 10 ^ pred)
      } else {
        err <- rbind(err, rep(0, length(qq)))
      }
    }
  }
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2",
                          c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  # Enforce monotonic decrease of error rate with increasing Q.
  for (r in seq_len(nrow(err))) {
    for (c in seq.int(ncol(err), 2)) {
      if (err[r, c] > err[r, c - 1]) err[r, c - 1] <- err[r, c]
    }
  }
  err[err < 0] <- 0
  err
}

t_total <- proc.time()

# 1. Filter (paired, lock-step, per-sample)
cat("[filterAndTrim]\n")
t <- proc.time()
fstats <- filterAndTrim(
  fwd_in,  fwd_filt,
  rev_in,  rev_filt,
  truncLen   = c(240L, 180L),
  maxEE      = c(2,    4),
  truncQ     = 2L,
  minLen     = 50L,
  compress   = TRUE,
  multithread = n_threads
)
t_filter <- (proc.time() - t)[3] * 1000
cat(sprintf("  total_in=%d  total_out=%d  (%.1f ms)\n",
            sum(fstats[, 1]), sum(fstats[, 2]), t_filter))

# 2. Learn errors (one model fwd, one rev)
# Try binned-loess first (quality-aware); fall back to dada2's noqualErrfun
# if EM diverges (happens on extremely binned NovaSeq rev reads).
learn_with_fallback <- function(fls, label) {
  res <- tryCatch(
    learnErrors(fls, multithread = n_threads, verbose = FALSE,
                MAX_CONSIST = 20L,
                errorEstimationFunction = loessErrfun_binned),
    error = function(e) {
      cat(sprintf("  %s: binned loess failed (%s); falling back to noqualErrfun\n",
                  label, conditionMessage(e)))
      learnErrors(fls, multithread = n_threads, verbose = FALSE,
                  MAX_CONSIST = 20L,
                  errorEstimationFunction = dada2:::noqualErrfun)
    })
  res
}
cat("[learnErrors]\n")
t <- proc.time()
errF <- learn_with_fallback(fwd_filt, "fwd")
errR <- learn_with_fallback(rev_filt, "rev")
t_errors <- (proc.time() - t)[3] * 1000
cat(sprintf("  errF + errR  (%.1f ms)\n", t_errors))

# 3. Dereplicate (per sample)
cat("[derepFastq]\n")
t <- proc.time()
derepF <- derepFastq(fwd_filt, verbose = FALSE)
derepR <- derepFastq(rev_filt, verbose = FALSE)
names(derepF) <- samples; names(derepR) <- samples
t_derep <- (proc.time() - t)[3] * 1000
cat(sprintf("  done  (%.1f ms)\n", t_derep))

# 4. DADA (pseudo-pool — standard dada2 cross-sample practice)
cat("[dada pool='pseudo']\n")
t <- proc.time()
dadaF <- dada(derepF, err = errF, pool = "pseudo",
              multithread = n_threads, verbose = FALSE)
dadaR <- dada(derepR, err = errR, pool = "pseudo",
              multithread = n_threads, verbose = FALSE)
t_dada <- (proc.time() - t)[3] * 1000
n_asvF <- sum(sapply(dadaF, function(d) length(d$denoised)))
n_asvR <- sum(sapply(dadaR, function(d) length(d$denoised)))
cat(sprintf("  fwd_asvs(total)=%d  rev_asvs(total)=%d  (%.1f ms)\n",
            n_asvF, n_asvR, t_dada))

# 5. Merge pairs (per sample)
cat("[mergePairs]\n")
t <- proc.time()
merged <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = FALSE)
t_merge <- (proc.time() - t)[3] * 1000
n_merged <- sum(sapply(merged, function(d) sum(d$accept)))
cat(sprintf("  merged_pairs=%d  (%.1f ms)\n", n_merged, t_merge))

# 6. Sequence table + chimera removal
cat("[makeSequenceTable + removeBimeraDenovo]\n")
t <- proc.time()
seqtab       <- makeSequenceTable(merged)
seqtab_clean <- removeBimeraDenovo(seqtab, method = "consensus",
                                   multithread = n_threads, verbose = FALSE)
t_chimera <- (proc.time() - t)[3] * 1000
cat(sprintf("  asvs_in=%d  asvs_out=%d  (%.1f ms)\n",
            ncol(seqtab), ncol(seqtab_clean), t_chimera))

t_total_ms <- (proc.time() - t_total)[3] * 1000
cat(sprintf("\nTotal R dada2 time: %.1f ms\n", t_total_ms))

# Per-sample read counts
sample_stats <- lapply(seq_along(samples), function(i) {
  list(sample     = samples[i],
       reads_in   = as.integer(fstats[i, 1]),
       reads_out  = as.integer(fstats[i, 2]))
})

# ASV table: rows = samples, cols = sequences
result <- list(
  tool       = "R dada2",
  n_threads  = n_threads,
  total_ms   = round(t_total_ms, 1),
  stages = list(
    filter_ms       = round(t_filter,  1),
    learn_errors_ms = round(t_errors,  1),
    derep_ms        = round(t_derep,   1),
    dada_ms         = round(t_dada,    1),
    merge_ms        = round(t_merge,   1),
    chimera_ms      = round(t_chimera, 1)
  ),
  samples = sample_stats,
  n_asvs_before_chimera = ncol(seqtab),
  n_asvs_after_chimera  = ncol(seqtab_clean),
  total_abundance       = sum(seqtab_clean),
  asvs = lapply(seq_len(ncol(seqtab_clean)), function(j)
    list(sequence  = colnames(seqtab_clean)[j],
         abundance = as.integer(sum(seqtab_clean[, j]))))
)
writeLines(toJSON(result, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "r_output.json"))

cat(sprintf("\nWrote %s\n", file.path(out_dir, "r_output.json")))
