# Aim 2b, step a: plan a partial download of the CADD v1.6 GRCh38 whole-genome SNV file (86.6 GB,
# not stored locally) restricted to the positions of the 2026 main reference sets. The tabix
# linear index (16 kb windows) gives, for each merged region, the compressed-file offset of the
# first BGZF block that can contain it; the byte range of a region runs from that offset to the
# offset of the first window after the region (plus one maximal BGZF block, 64 KiB).
suppressMessages(library(data.table))
tbi_file <- "C:/cadd/data/cadd16/whole_genome_SNVs.tsv.gz.tbi"
file_size <- 86592987071
gap <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(gap)) gap <- 20000L

read_tbi <- function(f) {
  con <- gzfile(f, "rb"); on.exit(close(con))
  magic <- readBin(con, "raw", 4); stopifnot(rawToChar(magic[1:3]) == "TBI")
  hdr <- readBin(con, "integer", 8, size = 4, endian = "little")
  n_ref <- hdr[1]; l_nm <- hdr[8]
  nmr <- readBin(con, "raw", l_nm); z <- which(nmr == as.raw(0))
  nm <- vapply(seq_along(z), function(j) rawToChar(nmr[(c(0, z)[j] + 1):(z[j] - 1)]), "")
  lin <- vector("list", n_ref)
  for (r in seq_len(n_ref)) {
    n_bin <- readBin(con, "integer", 1, size = 4, endian = "little")
    for (b in seq_len(n_bin)) {
      readBin(con, "integer", 1, size = 4, endian = "little")         # bin
      n_chunk <- readBin(con, "integer", 1, size = 4, endian = "little")
      if (n_chunk) readBin(con, "raw", 16 * n_chunk)                     # chunk begin/end (unused)
    }
    n_intv <- readBin(con, "integer", 1, size = 4, endian = "little")
    raw <- readBin(con, "raw", 8 * n_intv)
    m <- matrix(as.integer(raw), nrow = 8)                                # little-endian uint64 -> double
    lin[[r]] <- colSums(m * 256^(0:7))
  }
  setNames(lin, nm)
}
lin <- read_tbi(tbi_file)
cat("Sequences in index:", paste(names(lin), collapse = ","), "\n")

k <- c("chrom", "pos")
pos <- unique(rbind(readRDS("C:/cadd/data/refsets/clinvar_ref.rds")[, ..k], readRDS("C:/cadd/data/refsets/gnomad_ref.rds")[, ..k]))
setorder(pos, chrom, pos)
reg <- pos[, { p <- sort(pos); br <- c(TRUE, diff(p) > gap); g <- cumsum(br)
               .(start = p[br], end = vapply(split(p, g), max, 0L)) }, by = chrom]
cat("Positions:", nrow(pos), "| regions (gap", gap, "):", nrow(reg), "\n")
reg[, `:=`(off_start = NA_real_, off_end = NA_real_)]
for (ch in unique(reg$chrom)) {
  L <- lin[[ch]]; if (is.null(L)) next
  i <- which(reg$chrom == ch)
  ws <- pmin(reg$start[i] %/% 16384L + 1L, length(L)); we <- reg$end[i] %/% 16384L + 2L
  # a zero entry means no record starts in that window; use the previous non-zero entry
  Lf <- L; for (j in seq_along(Lf)[-1]) if (Lf[j] == 0) Lf[j] <- Lf[j - 1]
  reg$off_start[i] <- floor(Lf[ws] / 65536)
  reg$off_end[i] <- ifelse(we <= length(Lf), floor(Lf[pmin(we, length(Lf))] / 65536) + 65536, file_size)
}
reg <- reg[!is.na(off_start)]
# merge overlapping byte ranges
setorder(reg, off_start)
br <- reg[, .(s = off_start, e = pmin(off_end, file_size - 1))]
br[, grp := cumsum(c(TRUE, s[-1] > cummax(e)[-.N]))]
br <- br[, .(s = min(s), e = max(e)), by = grp]
cat(sprintf("Byte ranges: %d | total %.2f GB (%.1f%% of the file)\n", nrow(br), sum(br$e - br$s + 1) / 1e9,
            100 * sum(br$e - br$s + 1) / file_size))
saveRDS(list(regions = reg, ranges = br, gap = gap), sprintf("C:/cadd/data/cadd16/byte_ranges_gap%d.rds", gap))
