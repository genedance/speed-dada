#!/usr/bin/env Rscript
# Benchmark: dada2rs (Rust R-binding) on simulated 16S V3-V4 paired FASTQs.
# Mirrors bench_r.R but loads dada2rs instead of dada2.

.libPaths(c(path.expand("~/R/library"), .libPaths()))
suppressPackageStartupMessages(library(dada2rs))
suppressPackageStartupMessages(library(jsonlite))

R1  <- "/tmp/bench_fastq/R1.fastq"
R2  <- "/tmp/bench_fastq/R2.fastq"
out <- "/tmp/bench_out"
dir.create(out, showWarnings = FALSE)
fR1 <- file.path(out, "rs_filt_R1.fastq")
fR2 <- file.path(out, "rs_filt_R2.fastq")

t_total <- proc.time()

# 1. Filter and trim
cat("[filterAndTrim]\n")
t <- proc.time()
fstats <- filterAndTrim(R1, fR1, R2, fR2,
                         truncLen  = c(230L, 210L),
                         minLen    = 150L,
                         maxEE     = c(3, 5),
                         compress  = FALSE)
t_filter <- (proc.time() - t)[3] * 1000
cat(sprintf("  pairs_in=%d  pairs_out=%d  (%.1f ms)\n",
            fstats[1, 1], fstats[1, 2], t_filter))

# 2. Learn errors
cat("[learnErrors]\n")
t <- proc.time()
errF <- learnErrors(fR1)
errR <- learnErrors(fR2)
t_errors <- (proc.time() - t)[3] * 1000
cat(sprintf("  fwd+rev error models  (%.1f ms)\n", t_errors))

# 3. Dereplicate
cat("[derepFastq]\n")
t <- proc.time()
derepF <- derepFastq(fR1)
derepR <- derepFastq(fR2)
t_derep <- (proc.time() - t)[3] * 1000
cat(sprintf("  fwd_uniq=%d  rev_uniq=%d  (%.1f ms)\n",
            length(derepF$uniques), length(derepR$uniques), t_derep))

# 4. DADA
cat("[dada]\n")
t <- proc.time()
dadaF <- dada(derepF, err = errF)
dadaR <- dada(derepR, err = errR)
t_dada <- (proc.time() - t)[3] * 1000
cat(sprintf("  fwd_asvs=%d  rev_asvs=%d  (%.1f ms)\n",
            length(dadaF$denoised), length(dadaR$denoised), t_dada))

# 5. Merge
cat("[mergePairs]\n")
t <- proc.time()
merged <- mergePairs(dadaF, derepF, dadaR, derepR)
t_merge <- (proc.time() - t)[3] * 1000
cat(sprintf("  pairs=%d  accepted=%d  (%.1f ms)\n",
            nrow(merged), sum(merged$accept), t_merge))

# 6. Sequence table + chimera removal
cat("[makeSequenceTable + removeBimeraDenovo]\n")
t <- proc.time()
seqtab       <- makeSequenceTable(list(s1 = merged))
seqtab_clean <- removeBimeraDenovo(seqtab)
t_chimera <- (proc.time() - t)[3] * 1000
asvs   <- colnames(seqtab_clean)
abunds <- as.integer(seqtab_clean[1, ])
cat(sprintf("  asvs_out=%d  (%.1f ms)\n", length(asvs), t_chimera))

t_total_ms <- (proc.time() - t_total)[3] * 1000
cat(sprintf("\nTotal dada2rs time: %.1f ms\n", t_total_ms))

result <- list(
  tool       = "dada2rs",
  total_ms   = round(t_total_ms, 1),
  stages     = list(
    filter_ms       = round(t_filter,  1),
    learn_errors_ms = round(t_errors,  1),
    derep_ms        = round(t_derep,   1),
    dada_ms         = round(t_dada,    1),
    merge_ms        = round(t_merge,   1),
    chimera_ms      = round(t_chimera, 1)
  ),
  asvs = lapply(seq_along(asvs), function(i)
    list(sequence = asvs[i], abundance = abunds[i]))
)
writeLines(toJSON(result, auto_unbox = TRUE, pretty = TRUE),
           file.path(out, "dada2rs_output.json"))

cat("\nTop ASVs (dada2rs):\n")
ord <- order(-abunds)
for (i in head(ord, 5))
  cat(sprintf("  %s...  abundance=%d\n", substr(asvs[i], 1, 40), abunds[i]))
