suppressPackageStartupMessages(library(Rsamtools))
url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr17.vcf.bgz"
setwd(tempdir())   # htslib guarda el .tbi remoto en el directorio de trabajo
t0 <- Sys.time()
tf <- TabixFile(url)
x <- tryCatch(scanTabix(tf, param = GenomicRanges::GRanges("chr17", IRanges::IRanges(7675088, 7675088)))[[1]],
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("segundos:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n")
for (l in x) {
  f <- strsplit(l, "\t")[[1]]
  info <- f[8]
  af <- regmatches(info, regexpr("(^|;)AF=[^;]+", info))
  cat(f[1:5], af, "\n")
}
