#!/usr/bin/env Rscript
# Benchmark: dada2rs (Rust core via R extendr) on 3 paired JHI samples.

.libPaths(c(path.expand("~/R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(dada2rs)
  library(jsonlite)
})

args      <- commandArgs(trailingOnly = TRUE)
n_threads <- if (length(args) >= 1) as.integer(args[1]) else 16L
in_dir    <- if (length(args) >= 2) args[2] else stop("input dir required")
out_dir   <- if (length(args) >= 3) args[3] else "/tmp/bench_jhi_out/dada2rs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

Sys.setenv(RAYON_NUM_THREADS = as.character(n_threads))

stems   <- c("JHI-2025-Q1-A-004", "JHI-2025-Q1-A-009", "JHI-2025-Q1-A-010")
fwd_in  <- file.path(in_dir, paste0("raw.", stems, ".1.fq.gz"))
rev_in  <- file.path(in_dir, paste0("raw.", stems, ".2.fq.gz"))
fwd_filt <- file.path(out_dir, paste0(stems, "_R1_filt.fastq.gz"))
rev_filt <- file.path(out_dir, paste0(stems, "_R2_filt.fastq.gz"))
names(fwd_filt) <- stems
names(rev_filt) <- stems

cat(sprintf("[dada2rs] RAYON_NUM_THREADS=%d  samples=%d\n",
            n_threads, length(stems)))
t_total <- proc.time()

cat("[filterAndTrim]\n")
t <- proc.time()
fstats <- filterAndTrim(
  fwd_in, fwd_filt, rev_in, rev_filt,
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

cat("[learnErrors]\n")
t <- proc.time()
errF <- learnErrors(fwd_filt, multithread = n_threads)
errR <- learnErrors(rev_filt, multithread = n_threads)
t_errors <- (proc.time() - t)[3] * 1000
cat(sprintf("  done  (%.1f ms)\n", t_errors))

cat("[derepFastq]\n")
t <- proc.time()
derepF <- derepFastq(fwd_filt)
derepR <- derepFastq(rev_filt)
names(derepF) <- stems; names(derepR) <- stems
t_derep <- (proc.time() - t)[3] * 1000
cat(sprintf("  done  (%.1f ms)\n", t_derep))

cat("[dada pool='pseudo']\n")
t <- proc.time()
dadaF <- dada(derepF, err = errF, pool = "pseudo", multithread = n_threads)
dadaR <- dada(derepR, err = errR, pool = "pseudo", multithread = n_threads)
t_dada <- (proc.time() - t)[3] * 1000
cat(sprintf("  done  (%.1f ms)\n", t_dada))

cat("[mergePairs]\n")
t <- proc.time()
merged <- mapply(function(dF, drF, dR, drR) mergePairs(dF, drF, dR, drR),
                 dadaF, derepF, dadaR, derepR, SIMPLIFY = FALSE)
names(merged) <- stems
t_merge <- (proc.time() - t)[3] * 1000
cat(sprintf("  done  (%.1f ms)\n", t_merge))

cat("[makeSequenceTable + removeBimeraDenovo]\n")
t <- proc.time()
seqtab       <- makeSequenceTable(merged)
seqtab_clean <- removeBimeraDenovo(seqtab, multithread = n_threads)
t_chimera <- (proc.time() - t)[3] * 1000
cat(sprintf("  asvs_in=%d  asvs_out=%d  (%.1f ms)\n",
            ncol(seqtab), ncol(seqtab_clean), t_chimera))

t_total_ms <- (proc.time() - t_total)[3] * 1000
cat(sprintf("\nTotal dada2rs time: %.1f ms\n", t_total_ms))

sample_stats <- lapply(seq_along(stems), function(i) {
  list(sample = stems[i],
       reads_in  = as.integer(fstats[i, 1]),
       reads_out = as.integer(fstats[i, 2]))
})

result <- list(
  tool       = "dada2rs",
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
           file.path(out_dir, "dada2rs_output.json"))
cat(sprintf("\nWrote %s\n", file.path(out_dir, "dada2rs_output.json")))
