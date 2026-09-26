# Amendment A6: AlphaMissense, REVEL and BayesDel scores for reference variants of the 2019 and 2022
# release sets that are not yet in comparator_scores.rds, with exactly the rules of script 20
# (maximum over transcripts; BayesDel lifted to hg19 with reference-base check).
# Output: C:/cadd/data/refsets/comparator_scores_all_releases.rds (existing + new keys)
suppressPackageStartupMessages({ library(data.table); library(rtracklayer) })
rd <- "C:/cadd/data/refsets"; cd <- "C:/cadd/comparators"; rs <- "C:/cadd/calibration/results/v17_avi/release_series"
k <- c("chrom", "pos", "ref", "alt")
old <- readRDS(file.path(rd, "comparator_scores.rds"))
keys <- unique(rbindlist(lapply(c("r2019", "r2022"), function(r) { s <- readRDS(file.path(rs, sprintf("refsets_%s.rds", r)))
  rbind(s$clin[, ..k], s$gno[, ..k]) })))
keys <- keys[!old, on = k]
cat("new SNVs to score:", nrow(keys), "\n")

am <- fread(file.path(cd, "AlphaMissense_hg38.tsv.gz"), skip = "#CHROM", select = c(1, 2, 3, 4, 9),
            col.names = c("chrom", "pos", "ref", "alt", "AM"), colClasses = list(character = 1))
am[, chrom := sub("^chr", "", chrom)]; am <- am[keys, on = k, nomatch = NULL][, .(AM = max(AM)), by = k]; invisible(gc())
rv <- fread(cmd = sprintf('C:\\Windows\\System32\\tar.exe -xOf "%s" revel_with_transcript_ids',
                          normalizePath(file.path(cd, "revel-v1.3_all_chromosomes.zip"), winslash = "\\")),
            select = c("chr", "grch38_pos", "ref", "alt", "REVEL"), colClasses = list(character = c("chr", "grch38_pos")))
rv <- rv[grch38_pos != "."][, .(chrom = chr, pos = as.integer(grch38_pos), ref, alt, REVEL)]
rv <- rv[keys, on = k, nomatch = NULL][, .(REVEL = max(REVEL)), by = k]; invisible(gc())
chain <- import.chain(file.path(cd, "bayesdel", "hg38ToHg19.over.chain"))
lo <- liftOver(GRanges(paste0("chr", keys$chrom), IRanges(keys$pos, width = 1), strand = "+"), chain)
ok <- which(lengths(lo) == 1); l1 <- unlist(lo[ok]); comp <- c(A = "T", C = "G", G = "C", T = "A")
bd <- data.table(keys[ok], h19_chr = as.character(seqnames(l1)), h19_pos = start(l1), h19_strand = as.character(strand(l1)))
bd[, `:=`(ref19 = ifelse(h19_strand == "-", comp[ref], ref), alt19 = ifelse(h19_strand == "-", comp[alt], alt))]
tb <- TwoBitFile(file.path(cd, "bayesdel", "hg19.2bit")); bd <- bd[h19_chr %in% seqnames(seqinfo(tb))]
bd[, base19 := as.character(getSeq(tb, GRanges(h19_chr, IRanges(h19_pos, width = 1))))]
cat("liftOver: unmapped", sum(lengths(lo) != 1), "| reference mismatch", bd[base19 != ref19, .N], "\n")
bd <- bd[base19 == ref19][, BayesDel := NA_real_]
for (b in c("A", "C", "G", "T")) {
  idx <- which(bd$alt19 == b); if (!length(idx)) next
  q <- GRanges(bd$h19_chr[idx], IRanges(bd$h19_pos[idx], width = 1))
  x <- import(BigWigFile(file.path(cd, "bayesdel", sprintf("NoAFBayesDelMut%s.bw", b))), which = reduce(q))
  ov <- findOverlaps(q, x); bd$BayesDel[idx[queryHits(ov)]] <- score(x)[subjectHits(ov)]
}
new <- Reduce(function(a, b) merge(a, b, by = k, all.x = TRUE), list(keys, am, rv, bd[, c(k, "BayesDel"), with = FALSE]))
all <- rbind(old, new)
stopifnot(!anyDuplicated(all, by = k))
saveRDS(all, file.path(rd, "comparator_scores_all_releases.rds"))
cat("comparator table:", nrow(all), "SNVs (", nrow(new), "new )\n")
