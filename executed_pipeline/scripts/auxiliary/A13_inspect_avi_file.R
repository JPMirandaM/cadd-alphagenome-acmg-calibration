# Inspects the AlphaGenome AVI SNV tabix file: header lines, columns, sequence names, and the
# records for TP53 c.524G>A (chr17:7675088 C>T) and for the same site in CADD v1.7.
suppressPackageStartupMessages(library(Rsamtools))
f <- "C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz"
tf <- TabixFile(f)
cat("== header ==\n"); print(headerTabix(tf))
seqs <- seqnamesTabix(tf); cat("== sequences (", length(seqs), ") ==\n"); print(head(seqs, 30))
chr17 <- if ("chr17" %in% seqs) "chr17" else "17"
x <- scanTabix(tf, param = GenomicRanges::GRanges(chr17, IRanges::IRanges(7675088, 7675088)))[[1]]
cat("== AVI records at 17:7675088 ==\n"); writeLines(x)
y <- scanTabix(TabixFile("C:/cadd/data/whole_genome_SNVs.tsv.gz"),
               param = GenomicRanges::GRanges("17", IRanges::IRanges(7675088, 7675088)))[[1]]
cat("== CADD v1.7 records at 17:7675088 ==\n"); writeLines(y)
