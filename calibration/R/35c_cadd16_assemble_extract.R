# Aim 2b (A14), step c: extract CADD v1.6 scores for the 2026 main reference sets from the
# downloaded byte ranges of the official whole-genome file.
# Each range starts at a BGZF block boundary (from the tabix linear index) and ends one maximal BGZF
# block (64 KiB) after the block holding the first record beyond the region. For each range, the
# complete BGZF blocks are kept (a truncated trailing block is dropped) and decompressed; because
# BGZF blocks do not respect line boundaries, the first line (may start mid-record, and precedes the
# region) and the last line (may be cut, and follows the region) of each range are discarded. The
# remaining records are matched to the reference SNVs by chromosome, position and alleles.
# (A first version concatenated the blocks of all ranges into one file for tabix indexing; this
# failed because records were cut at range junctions; log 35c_..._run1_error.log.)
suppressMessages({ library(data.table); library(parallel) })
dir16 <- "C:/cadd/data/cadd16"
rng <- fread(file.path(dir16, "byte_ranges.tsv"), colClasses = c("integer", "numeric", "numeric"))
k <- c("chrom", "pos", "ref", "alt")
keys <- unique(rbind(readRDS("C:/cadd/data/refsets/clinvar_ref.rds")[, ..k], readRDS("C:/cadd/data/refsets/gnomad_ref.rds")[, ..k]))
tmpd <- file.path(dir16, "tmp_extract"); dir.create(tmpd, showWarnings = FALSE)

extract_part <- function(i) {
  outf <- file.path(tmpd, sprintf("part_%05d.rds", rng$id[i]))
  if (file.exists(outf)) return(readRDS(outf)$stats)
  f <- file.path(dir16, "parts", sprintf("part_%05d.bin", rng$id[i]))
  stopifnot(file.size(f) == rng$e[i] - rng$s[i] + 1)
  x <- readBin(f, "raw", file.size(f)); p <- 1L; n <- length(x); keep_to <- 0L; nb <- 0L
  while (p + 17L <= n) {
    if (!(x[p] == as.raw(0x1f) && x[p + 1L] == as.raw(0x8b) && x[p + 3L] == as.raw(0x04))) stop(sprintf("not a BGZF block start: part %d offset %d", rng$id[i], p))
    bsize <- as.integer(x[p + 16L]) + 256L * as.integer(x[p + 17L]) + 1L
    if (p + bsize - 1L > n) break
    keep_to <- p + bsize - 1L; nb <- nb + 1L; p <- p + bsize
  }
  tf <- tempfile(fileext = ".gz"); writeBin(x[1:keep_to], tf)
  con <- gzfile(tf, "r"); l <- readLines(con, warn = FALSE); close(con); unlink(tf)
  l <- l[-c(1L, length(l))]
  d <- fread(text = l, sep = "\t", header = FALSE, col.names = c("chrom", "pos", "ref", "alt", "CADD16_raw", "CADD16_PHRED"),
             colClasses = list(character = 1L), nThread = 1L)
  hit <- d[keys, on = k, nomatch = NULL]
  stats <- data.table(id = rng$id[i], blocks = nb, lines = length(l), first = paste(d$chrom[1], d$pos[1]), last = paste(d$chrom[nrow(d)], d$pos[nrow(d)]), hits = nrow(hit))
  saveRDS(list(hit = hit, stats = stats), outf)
  stats
}

cl <- makeCluster(4L)
clusterEvalQ(cl, suppressMessages(library(data.table)))
clusterExport(cl, c("dir16", "rng", "k", "keys", "tmpd", "extract_part"))
t0 <- Sys.time()
st <- rbindlist(parLapplyLB(cl, seq_len(nrow(rng)), extract_part)); stopCluster(cl)
cat(sprintf("Ranges processed: %d | BGZF blocks: %d | records read: %.0f | reference hits: %d | %.1f min\n", nrow(st), sum(st$blocks),
            sum(as.numeric(st$lines)), sum(st$hits), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
fwrite(st, file.path(dir16, "extract_stats.tsv"), sep = "\t")
res <- unique(rbindlist(lapply(list.files(tmpd, "^part_.*\\.rds$", full.names = TRUE), function(f) readRDS(f)$hit)), by = k)
sc <- merge(keys, res, by = k, all.x = TRUE)
cat("Reference SNVs:", nrow(keys), "| with a CADD v1.6 score:", sum(!is.na(sc$CADD16_PHRED)), "\n")
# consistency check: positions covered by a range but lacking a v1.6 record for the allele
cat("Unscored reference SNVs by chromosome:\n"); print(sc[is.na(CADD16_PHRED), .N, by = chrom][order(-N)][1:10])
saveRDS(sc, "C:/cadd/data/refsets/cadd16_scores_refsets.rds")
