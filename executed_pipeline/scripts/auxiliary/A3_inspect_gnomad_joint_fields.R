suppressPackageStartupMessages(library(Rsamtools))
setwd(tempdir())
url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr17.vcf.bgz"
x <- scanTabix(TabixFile(url), param = GenomicRanges::GRanges("chr17", IRanges::IRanges(7675088, 7675088)))[[1]]
cat("bytes por linea:", nchar(x), "\n")
f <- strsplit(x[1], "\t")[[1]]
cat(f[1:7], "\n")
keys <- strsplit(f[8], ";")[[1]]
cat("n campos INFO:", length(keys), "\n")
print(grep("^(AF|AC|AN|nhomalt|grpmax|fafmax|faf95)[_=]", keys, value = TRUE)[1:60])
