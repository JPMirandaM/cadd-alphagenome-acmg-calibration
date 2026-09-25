suppressPackageStartupMessages({ library(data.table); library(Rsamtools); library(GenomicRanges) })
gtf <- fread(cmd = NULL, file = "C:/cadd/annot/gencode.v50.basic.annotation.gtf.gz", sep = "\t",
             header = FALSE, skip = "chr1\t", select = c(1, 3, 4, 5),
             col.names = c("chr", "type", "start", "end"))
cds <- gtf[type == "CDS"]
cds_all <- reduce(GRanges(cds$chr, IRanges(cds$start - 10, cds$end + 10)))
cat("CDS+-10 genoma: ", round(sum(width(cds_all)) / 1e6, 1), "Mb en", length(cds_all), "regiones\n")
r21 <- cds_all[seqnames(cds_all) == "chr21"]
cat("chr21:", round(sum(width(r21)) / 1e6, 2), "Mb en", length(r21), "regiones\n")

setwd(tempdir())
url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr21.vcf.bgz"
tf <- open(TabixFile(url))
sub <- r21[1:300]
t0 <- Sys.time()
x <- scanTabix(tf, param = sub)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
nl <- sum(lengths(x)); mb <- sum(vapply(x, function(v) sum(nchar(v, "bytes")), 0)) / 1e6
cat(sprintf("300 regiones (%.3f Mb): %.1f s, %d lineas, %.1f MB\n", sum(width(sub)) / 1e6, secs, nl, mb))
cat(sprintf("Extrapolado a todo el CDS (%.1f Mb), 1 proceso: %.1f h\n",
            sum(width(cds_all)) / 1e6, secs * sum(width(cds_all)) / sum(width(sub)) / 3600))
