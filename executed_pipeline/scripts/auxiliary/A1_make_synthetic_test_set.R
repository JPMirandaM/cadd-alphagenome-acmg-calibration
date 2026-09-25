suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
set.seed(1)
d <- rbindlist(lapply(c("1", "2", "X"), function(ch) {
  pos <- rep(seq(10001L, 4500000L, by = 3L), each = 3)   # cruza varias ventanas de 2 Mb
  data.table(Chrom = ch, Pos = pos, Ref = "A", Alt = c("C", "G", "T"),
             RawScore = round(rnorm(length(pos)), 6))
}))
d[, PHRED := round(pmax(0, RawScore * 12 + 5), 3)]
# empates forzados: mismas PHRED y RawScore en posiciones distintas
d[c(100, 5000, 900000), `:=`(RawScore = 9.5, PHRED = 99)]
dir <- "<SCRATCH>/test"
unlink(dir, recursive = TRUE); dir.create(dir)
tsv <- file.path(dir, "t.tsv")
writeLines(c("## CADD GRCh38-v1.7 test", "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED"), tsv)
fwrite(d, tsv, sep = "\t", append = TRUE, col.names = FALSE)
gz <- bgzip(tsv, overwrite = TRUE); indexTabix(gz, seq = 1, start = 2, end = 2, comment = "#")
# referencia: orden directo en memoria
chr_levels <- c(as.character(1:22), "X", "Y", "MT")
d[, o := match(Chrom, chr_levels)]
setorder(d, -PHRED, -RawScore, o, Pos, Ref, Alt)
fwrite(d[1:20000, .(Rank = .I, Chrom, Pos, Ref, Alt, RawScore, PHRED)], file.path(dir, "expected.tsv"), sep = "\t")
cat("filas:", nrow(d), "\n")
