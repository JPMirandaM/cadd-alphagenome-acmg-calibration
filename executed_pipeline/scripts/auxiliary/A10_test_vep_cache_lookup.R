# Checks the column layout of the Ensembl VEP 116 indexed cache (all_vars.gz) and compares the
# gnomAD v4.1 frequencies it stores for TP53 c.524G>A (p.R175H; chr17:7675088 C>T) with those
# obtained directly from the gnomAD v4.1 exomes VCF (A2: AF_exomes = 4.10437e-06).
suppressPackageStartupMessages(library(Rsamtools))
info <- readLines("C:/cadd/vep/homo_sapiens/116_GRCh38/info.txt")
cols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", info, value = TRUE)), ",")[[1]]
f <- "C:/cadd/vep/homo_sapiens/116_GRCh38/17/all_vars.gz"
x <- scanTabix(TabixFile(f, index = paste0(f, ".csi")), param = GenomicRanges::GRanges("17", IRanges::IRanges(7675088, 7675088)))[[1]]
cat("records:", length(x), "| fields per record:", length(strsplit(x[1], "\t")[[1]]), "| variation_cols:", length(cols), "\n")
for (l in x) {
  v <- strsplit(l, "\t")[[1]]
  names(v) <- cols[seq_along(v)]
  print(v[c("chr", "variation_name", "start", "allele_string", "clin_sig", "gnomADe", "gnomADg", "gnomADe_NFE")])
}
