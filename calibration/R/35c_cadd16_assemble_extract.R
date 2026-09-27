# Aim 2b (A14), step c: assemble the downloaded CADD v1.6 byte ranges into one valid BGZF file
# (complete blocks only, in file order, plus the BGZF end-of-file block), tabix-index it, and extract
# CADD v1.6 scores for the 2026 main reference sets (exact position and alleles).
suppressMessages({ library(data.table); library(Rsamtools); library(GenomicRanges) })
dir16 <- "C:/cadd/data/cadd16"
rng <- fread(file.path(dir16, "byte_ranges.tsv"), colClasses = c("integer", "numeric", "numeric"))
outf <- file.path(dir16, "cadd16_refset_subset.tsv.gz")
con <- file(outf, "wb")
n_blocks <- 0; n_trunc <- 0; prev_end <- -1
for (i in seq_len(nrow(rng))) {
  f <- file.path(dir16, "parts", sprintf("part_%05d.bin", rng$id[i]))
  stopifnot(file.exists(f), file.size(f) == rng$e[i] - rng$s[i] + 1)
  x <- readBin(f, "raw", file.size(f)); p <- 1L; n <- length(x)
  if (rng$s[i] < prev_end) stop("overlapping ranges")        # ranges are disjoint and sorted
  keep_to <- 0L
  while (p + 17L <= n) {
    if (!(x[p] == as.raw(0x1f) && x[p + 1L] == as.raw(0x8b) && x[p + 3L] == as.raw(0x04))) stop(sprintf("not a BGZF block start: part %d offset %d", rng$id[i], p))
    bsize <- as.integer(x[p + 16L]) + 256L * as.integer(x[p + 17L]) + 1L
    if (p + bsize - 1L > n) { n_trunc <- n_trunc + 1; break }
    keep_to <- p + bsize - 1L; n_blocks <- n_blocks + 1; p <- p + bsize
  }
  if (keep_to > 0) writeBin(x[1:keep_to], con)
  prev_end <- rng$e[i]
}
eof <- as.raw(c(0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0))
writeBin(eof, con); close(con)
cat("BGZF blocks kept:", n_blocks, "| truncated trailing blocks dropped:", n_trunc, "| file GB:", round(file.size(outf) / 1e9, 2), "\n")
idx <- indexTabix(outf, seq = 1L, start = 2L, end = 2L, comment = "#")
cat("Index:", idx, "\n")

k <- c("chrom", "pos", "ref", "alt")
keys <- unique(rbind(readRDS("C:/cadd/data/refsets/clinvar_ref.rds")[, ..k], readRDS("C:/cadd/data/refsets/gnomad_ref.rds")[, ..k]))
tf <- TabixFile(outf); sn <- seqnamesTabix(tf)
res <- rbindlist(lapply(intersect(unique(keys$chrom), sn), function(ch) {
  p <- unique(keys[chrom == ch, pos]); gr <- reduce(GRanges(ch, IRanges(p, p)), min.gapwidth = 1L)
  l <- unlist(scanTabix(tf, param = gr), use.names = FALSE)
  x <- tstrsplit(l, "\t", fixed = TRUE)
  data.table(chrom = x[[1]], pos = as.integer(x[[2]]), ref = x[[3]], alt = x[[4]], CADD16_raw = as.numeric(x[[5]]), CADD16_PHRED = as.numeric(x[[6]]))
}))
res <- unique(res, by = k)
sc <- merge(keys, res, by = k, all.x = TRUE)
cat("Reference SNVs:", nrow(keys), "| with a CADD v1.6 score:", sum(!is.na(sc$CADD16_PHRED)), "\n")
# positions present in the subset but with no matching allele (reference mismatch check)
pos_hit <- unique(res[, .(chrom, pos)])
miss <- sc[is.na(CADD16_PHRED)][pos_hit, on = .(chrom, pos), nomatch = NULL]
cat("Unscored SNVs at positions present in the subset (allele mismatch):", nrow(miss), "\n")
saveRDS(sc, "C:/cadd/data/refsets/cadd16_scores_refsets.rds")
