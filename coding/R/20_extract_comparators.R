# Aim 6: scores of the already-recommended predictors for every variant in the 2026 and 2024
# reference sets (ClinVar and gnomAD sets):
#   AlphaMissense (hg38, Zenodo 10813168)          -> AM
#   REVEL v1.3 (grch38_pos column of the official file)  -> REVEL
#   BayesDel noAF (UCSC hg19 bigWig per alternative allele) -> BayesDel, after liftOver
#     hg38 -> hg19 (UCSC chain) and a check of the reference base against hg19.2bit
#     (minus-strand mappings: alleles complemented before the check and the lookup).
# Where a predictor has several values for the same SNV (several transcripts), the maximum is
# used and the number of such SNVs is reported.
# Output: C:/cadd/data/refsets/comparator_scores.rds and coverage summary.
suppressPackageStartupMessages({ library(data.table); library(rtracklayer) })
rd <- "C:/cadd/data/refsets"; cd <- "C:/cadd/comparators"
k <- c("chrom", "pos", "ref", "alt")
s24 <- readRDS("C:/cadd/calibration/results/v17_avi/temporal_holdout/refsets_2024.rds")
keys <- unique(rbindlist(list(readRDS(file.path(rd, "clinvar_ref.rds"))[, ..k], readRDS(file.path(rd, "gnomad_ref.rds"))[, ..k],
                              s24$clin[, ..k], s24$gno[, ..k])))
rm(s24); invisible(gc())
cat("unique SNVs to score:", nrow(keys), "\n")

## AlphaMissense ---------------------------------------------------------------------------------
# Intermediate results are cached so that a failure in a later step does not repeat these reads
cache_am <- file.path(rd, "tmp_am.rds"); cache_rv <- file.path(rd, "tmp_revel.rds")
if (file.exists(cache_am)) { am <- readRDS(cache_am); cat("AlphaMissense: loaded from cache\n") } else {
am <- fread(file.path(cd, "AlphaMissense_hg38.tsv.gz"), skip = "#CHROM", select = c(1, 2, 3, 4, 9),
            col.names = c("chrom", "pos", "ref", "alt", "AM"), colClasses = list(character = 1))
am[, chrom := sub("^chr", "", chrom)]
am <- am[keys, on = k, nomatch = NULL]
cat("AlphaMissense: SNVs with >1 value:", am[, .N, by = k][N > 1, .N], "\n")
am <- am[, .(AM = max(AM)), by = k]; saveRDS(am, cache_am); invisible(gc()) }

## REVEL (streamed from the zip; GRCh38 positions) ----------------------------------------------
# Windows bsdtar by full path (a GNU tar found earlier in PATH reads "C:" as a remote host)
if (file.exists(cache_rv)) { rv <- readRDS(cache_rv); cat("REVEL: loaded from cache\n") } else {
rv <- fread(cmd = sprintf('C:\\Windows\\System32\\tar.exe -xOf "%s" revel_with_transcript_ids',
                          normalizePath(file.path(cd, "revel-v1.3_all_chromosomes.zip"), winslash = "\\")),
            select = c("chr", "grch38_pos", "ref", "alt", "REVEL"), colClasses = list(character = c("chr", "grch38_pos")))
rv <- rv[grch38_pos != "."][, .(chrom = chr, pos = as.integer(grch38_pos), ref, alt, REVEL)]
rv <- rv[keys, on = k, nomatch = NULL]
cat("REVEL: SNVs with >1 value:", rv[, .N, by = k][N > 1, .N], "\n")
rv <- rv[, .(REVEL = max(REVEL)), by = k]; saveRDS(rv, cache_rv); invisible(gc()) }

## BayesDel noAF via liftOver to hg19 -----------------------------------------------------------
chain_file <- file.path(cd, "bayesdel", "hg38ToHg19.over.chain")         # import.chain needs it uncompressed
if (!file.exists(chain_file)) R.utils::gunzip(paste0(chain_file, ".gz"), destname = chain_file, remove = FALSE)
chain <- import.chain(chain_file)
gr <- GRanges(paste0("chr", keys$chrom), IRanges(keys$pos, width = 1), strand = "+")
lo <- liftOver(gr, chain)
n_map <- lengths(lo)
cat("liftOver hg38 -> hg19: 1 hit", sum(n_map == 1), "| 0 hits", sum(n_map == 0), "| >1 hits", sum(n_map > 1), "\n")
ok <- which(n_map == 1)
l1 <- unlist(lo[ok])
comp <- c(A = "T", C = "G", G = "C", T = "A")
bd <- data.table(keys[ok], h19_chr = as.character(seqnames(l1)), h19_pos = start(l1), h19_strand = as.character(strand(l1)))
bd[, `:=`(ref19 = ifelse(h19_strand == "-", comp[ref], ref), alt19 = ifelse(h19_strand == "-", comp[alt], alt))]
tb <- TwoBitFile(file.path(cd, "bayesdel", "hg19.2bit"))
bd <- bd[h19_chr %in% seqnames(seqinfo(tb))]
bd[, base19 := as.character(getSeq(tb, GRanges(h19_chr, IRanges(h19_pos, width = 1))))]
cat("minus-strand mappings:", bd[h19_strand == "-", .N], "| reference base mismatch after liftOver (excluded):", bd[base19 != ref19, .N], "\n")
bd <- bd[base19 == ref19]
# Only intervals present in the bigWig are imported (import with 'which'), so positions without a
# record stay NA; they are not confused with a genuine score of 0.
bd[, BayesDel := NA_real_]
for (b in c("A", "C", "G", "T")) {
  idx <- which(bd$alt19 == b)
  if (!length(idx)) next
  q <- GRanges(bd$h19_chr[idx], IRanges(bd$h19_pos[idx], width = 1))
  x <- import(BigWigFile(file.path(cd, "bayesdel", sprintf("NoAFBayesDelMut%s.bw", b))), which = reduce(q))
  ov <- findOverlaps(q, x)
  bd$BayesDel[idx[queryHits(ov)]] <- score(x)[subjectHits(ov)]
  rm(x, q, ov); invisible(gc())
}
cat("BayesDel: lifted SNVs with a score:", bd[!is.na(BayesDel), .N], "of", nrow(bd), "| exact zeros:", bd[BayesDel == 0, .N], "\n")
bdk <- bd[, c(k, "BayesDel"), with = FALSE]

## Combine and report coverage --------------------------------------------------------------------
sc <- Reduce(function(a, b) merge(a, b, by = k, all.x = TRUE), list(keys, am, rv, bdk))
saveRDS(sc, file.path(rd, "comparator_scores.rds"))
clin <- readRDS(file.path(rd, "clinvar_ref.rds"))[, c(k, "y"), with = FALSE]
cov_tab <- merge(clin, sc, by = k)[, .(n = .N, AM = sum(!is.na(AM)), REVEL = sum(!is.na(REVEL)),
                                       BayesDel = sum(!is.na(BayesDel)),
                                       all_three = sum(!is.na(AM) & !is.na(REVEL) & !is.na(BayesDel))), by = y]
cat("Coverage in the 2026 ClinVar main set (y = 1 P/LP, 0 B/LB):\n"); print(cov_tab)
fwrite(cov_tab, file.path(rd, "comparator_coverage_clinvar2026.tsv"), sep = "\t")
unlink(c(cache_am, cache_rv))
