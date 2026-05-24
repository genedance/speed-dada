#!/usr/bin/env Rscript
# Benchmark R dada2 pipeline on simulated paired-end FASTQs.
# Pipeline: filterAndTrim → learnErrors → derepFastq → dada → mergePairs → removeBimeraDenovo
# Taxonomy is skipped.

.libPaths(c(path.expand("~/R/library"), .libPaths()))
suppressPackageStartupMessages(library(dada2))

R1     <- "/tmp/bench_fastq/R1.fastq"
R2     <- "/tmp/bench_fastq/R2.fastq"
out    <- "/tmp/bench_out"
fR1    <- file.path(out, "r_filt_R1.fastq")
fR2    <- file.path(out, "r_filt_R2.fastq")
result_json <- file.path(out, "r_output.json")

t_total <- proc.time()

# 1. Filter and trim (paired, lock-step)
cat("[filterAndTrim]\n")
t <- proc.time()
filter_stats <- filterAndTrim(R1, fR1, R2, fR2,
                               truncLen=c(110, 110), minLen=50,
                               maxEE=c(3, 3), compress=FALSE,
                               multithread=TRUE)
cat(sprintf("  pairs_in=%d  pairs_out=%d  (%.1f ms)\n",
            filter_stats[1,1], filter_stats[1,2],
            (proc.time()-t)[3]*1000))

# 2. Learn errors
cat("[learnErrors]\n")
t <- proc.time()
errF <- learnErrors(fR1, multithread=TRUE)
errR <- learnErrors(fR2, multithread=TRUE)
cat(sprintf("  fwd+rev error models  (%.1f ms)\n", (proc.time()-t)[3]*1000))

# 3. Derep
cat("[derepFastq]\n")
t <- proc.time()
derepF <- derepFastq(fR1)
derepR <- derepFastq(fR2)
cat(sprintf("  fwd_uniq=%d  rev_uniq=%d  (%.1f ms)\n",
            length(derepF$uniques), length(derepR$uniques),
            (proc.time()-t)[3]*1000))

# 4. DADA
cat("[dada]\n")
t <- proc.time()
dadaF <- dada(derepF, err=errF, multithread=TRUE)
dadaR <- dada(derepR, err=errR, multithread=TRUE)
cat(sprintf("  fwd_asvs=%d  rev_asvs=%d  (%.1f ms)\n",
            nrow(dadaF$denoised), nrow(dadaR$denoised),
            (proc.time()-t)[3]*1000))

# 5. Merge
cat("[mergePairs]\n")
t <- proc.time()
merged <- mergePairs(dadaF, derepF, dadaR, derepR)
cat(sprintf("  merged_reads=%d  (%.1f ms)\n",
            sum(merged$accept), (proc.time()-t)[3]*1000))

# 6. Sequence table + remove bimeras
cat("[makeSequenceTable + removeBimeraDenovo]\n")
t <- proc.time()
seqtab <- makeSequenceTable(list(s1=merged))
seqtab_clean <- removeBimeraDenovo(seqtab, multithread=TRUE)
asvs   <- colnames(seqtab_clean)
abunds <- as.integer(seqtab_clean[1,])
cat(sprintf("  asvs_out=%d  (%.1f ms)\n",
            length(asvs), (proc.time()-t)[3]*1000))

total_ms <- (proc.time()-t_total)[3]*1000
cat(sprintf("\nTotal R dada2 pipeline time: %.1f ms\n", total_ms))

# Save results
results <- lapply(seq_along(asvs), function(i)
  list(sequence=asvs[i], abundance=abunds[i]))
library(jsonlite)
writeLines(toJSON(results, auto_unbox=TRUE, pretty=TRUE), result_json)

cat("\nTop ASVs (R):\n")
ord <- order(-abunds)
for (i in head(ord, 5)) {
  cat(sprintf("  %s...  abundance=%d\n", substr(asvs[i], 1, 40), abunds[i]))
}
