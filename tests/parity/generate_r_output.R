#!/usr/bin/env Rscript
# Generate reference R dada2 output for the parity test.
# Usage: Rscript generate_r_output.R
# Requires: dada2 R package (BioConductor)

suppressPackageStartupMessages(library(dada2))

set.seed(42)
TRUE_SEQ <- paste(rep(c("A","C","G","T"), 15), collapse="")  # 60 bp

# Generate 1000 reads with 2% errors
mutate_read <- function(seq, err=0.02) {
  bases <- strsplit(seq, "")[[1]]
  for (i in seq_along(bases)) {
    if (runif(1) < err) {
      alt <- setdiff(c("A","C","G","T"), bases[i])
      bases[i] <- sample(alt, 1)
    }
  }
  paste(bases, collapse="")
}

reads  <- sapply(1:1000, function(i) mutate_read(TRUE_SEQ))
quals  <- strrep("I", nchar(TRUE_SEQ))

tmp_fq <- tempfile(fileext=".fastq")
con <- file(tmp_fq, "w")
for (i in seq_along(reads)) {
  writeLines(c(paste0("@read_", i), reads[i], "+", quals), con)
}
close(con)

# Run pipeline
filtered_fq <- tempfile(fileext=".fastq")
out <- filterAndTrim(tmp_fq, filtered_fq, truncLen=nchar(TRUE_SEQ),
                     maxEE=5, minLen=20)
err  <- learnErrors(filtered_fq, multithread=TRUE)
drp  <- derepFastq(filtered_fq)
dada_res <- dada(drp, err=err)

# Format output
seqtab <- makeSequenceTable(list(s1=dada_res))
asvs   <- colnames(seqtab)
abunds <- as.integer(seqtab[1,])

result <- lapply(seq_along(asvs), function(i)
  list(sequence=asvs[i], abundance=abunds[i]))

library(jsonlite)
out_json <- file.path(dirname(sys.frame(1)$ofile), "fixtures", "r_output.json")
write(toJSON(result, auto_unbox=TRUE, pretty=TRUE), out_json)
cat("Written:", out_json, "\n")
