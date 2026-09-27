# Aim 3 (A13), step c: evaluation of the calibrated thresholds against MAVE functional classes.
# Primary: CADD v1.7 and AVI (main-analysis median thresholds). Secondary: REVEL, BayesDel,
# AlphaMissense (thresholds from the intersection calibration, Aim 6), on SNVs scored by all five.
# Metrics per evidence interval: proportion abnormal and LR = P(interval | abnormal) /
# P(interval | normal) with 95% bootstrap CI (B = 2,000, resampling SNVs); AUROC pooled and per gene
# (>= 30 abnormal and >= 30 normal); paired AVI - CADD AUROC difference.
# Sensitivity: SNVs in the ClinVar calibration reference set excluded.
suppressPackageStartupMessages({ library(data.table); library(rtracklayer) })
set.seed(20260927)
mdir <- "C:/cadd/data/mave"; cdir <- "C:/cadd/data/coding"; cd <- "C:/cadd/comparators"
res_dir <- "C:/cadd/calibration/results/v17_avi/mave"; B <- 2000L
k <- c("chrom", "pos", "ref", "alt")

lab <- readRDS(file.path(mdir, "mave_labels_mane_missense.rds"))
cat("Labelled MANE missense SNVs:", nrow(lab), "| abnormal:", sum(lab$abnormal), "| normal:", sum(lab$abnormal == 0),
    "| genes:", uniqueN(lab$gene), "\n")

## Scores -------------------------------------------------------------------------------------------
chs <- unique(lab$chrom)
sc <- rbindlist(lapply(chs, function(ch) {
  a <- readRDS(file.path(cdir, sprintf("coding_snv_%s.rds", ch)))[, .(chrom, pos, ref, alt, CADD = PHRED)]
  v <- readRDS(file.path(cdir, sprintf("avi_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, AVI = AVI_PHRED)]
  x <- merge(lab[chrom == ch, ..k], unique(a, by = k), by = k, all.x = TRUE)
  merge(x, unique(v, by = k), by = k, all.x = TRUE)
}))
cache <- file.path(mdir, "comparators_mave.rds")
if (file.exists(cache)) comp <- readRDS(cache) else {
  keys <- unique(lab[, ..k])
  am <- fread(file.path(cd, "AlphaMissense_hg38.tsv.gz"), skip = "#CHROM", select = c(1, 2, 3, 4, 9),
              col.names = c("chrom", "pos", "ref", "alt", "AM"), colClasses = list(character = 1))
  am[, chrom := sub("^chr", "", chrom)]; am <- am[keys, on = k, nomatch = NULL][, .(AM = max(AM)), by = k]; invisible(gc())
  rv <- fread(cmd = sprintf('C:\\Windows\\System32\\tar.exe -xOf "%s" revel_with_transcript_ids',
                            normalizePath(file.path(cd, "revel-v1.3_all_chromosomes.zip"), winslash = "\\")),
              select = c("chr", "grch38_pos", "ref", "alt", "REVEL"), colClasses = list(character = c("chr", "grch38_pos")))
  rv <- rv[grch38_pos != "."][, .(chrom = chr, pos = as.integer(grch38_pos), ref, alt, REVEL)]
  rv <- rv[keys, on = k, nomatch = NULL][, .(REVEL = max(REVEL)), by = k]; invisible(gc())
  chain_file <- file.path(cd, "bayesdel", "hg38ToHg19.over.chain"); chain <- import.chain(chain_file)
  lo <- liftOver(GRanges(paste0("chr", keys$chrom), IRanges(keys$pos, width = 1), strand = "+"), chain)
  ok <- which(lengths(lo) == 1); l1 <- unlist(lo[ok]); comp_b <- c(A = "T", C = "G", G = "C", T = "A")
  bd <- data.table(keys[ok], h19_chr = as.character(seqnames(l1)), h19_pos = start(l1), h19_strand = as.character(strand(l1)))
  bd[, `:=`(ref19 = ifelse(h19_strand == "-", comp_b[ref], ref), alt19 = ifelse(h19_strand == "-", comp_b[alt], alt))]
  tb <- TwoBitFile(file.path(cd, "bayesdel", "hg19.2bit")); bd <- bd[h19_chr %in% seqnames(seqinfo(tb))]
  bd[, base19 := as.character(getSeq(tb, GRanges(h19_chr, IRanges(h19_pos, width = 1))))]; bd <- bd[base19 == ref19]
  bd[, BayesDel := NA_real_]
  for (b in c("A", "C", "G", "T")) { idx <- which(bd$alt19 == b); if (!length(idx)) next
    q <- GRanges(bd$h19_chr[idx], IRanges(bd$h19_pos[idx], width = 1))
    xx <- import(BigWigFile(file.path(cd, "bayesdel", sprintf("NoAFBayesDelMut%s.bw", b))), which = reduce(q))
    ov <- findOverlaps(q, xx); bd$BayesDel[idx[queryHits(ov)]] <- score(xx)[subjectHits(ov)] }
  comp <- Reduce(function(a, b) merge(a, b, by = k, all.x = TRUE), list(keys, am, rv, bd[, c(k, "BayesDel"), with = FALSE]))
  saveRDS(comp, cache)
}
d <- merge(merge(lab, sc, by = k), comp, by = k, all.x = TRUE)
cat("Scored: CADD", sum(!is.na(d$CADD)), "| AVI", sum(!is.na(d$AVI)), "| AM", sum(!is.na(d$AM)),
    "| REVEL", sum(!is.na(d$REVEL)), "| BayesDel", sum(!is.na(d$BayesDel)), "\n")

## Thresholds ---------------------------------------------------------------------------------------
num1 <- function(x) suppressWarnings(as.numeric(sub(" .*", "", x)))
tm <- fread("C:/cadd/calibration/results/v17_avi/summary/table_main_thresholds.tsv")
tc <- fread("C:/cadd/calibration/results/v17_avi/comparators/summary/table_thresholds.tsv")
th <- rbind(tm[tool %in% c("CADD", "AVI"), .(tool, level, PP3 = num1(PP3), BP4 = num1(BP4))],
            tc[tool %in% c("REVEL", "BayesDel", "AlphaMissense"), .(tool, level, PP3 = num1(PP3), BP4 = num1(BP4))])
th[tool == "AlphaMissense", tool := "AM"]
fwrite(th, file.path(res_dir, "thresholds_used.tsv"), sep = "\t")
lv <- c("Supporting", "Moderate", "Strong")
interval_of <- function(x, tool) {
  t <- th[tool == ..tool]; out <- rep("Indeterminate", length(x)); out[is.na(x)] <- NA
  for (l in lv) { p <- t[level == l, PP3]; b <- t[level == l, BP4]
    if (length(p) && !is.na(p)) out[!is.na(x) & x >= p] <- paste0("PP3_", l)
    if (length(b) && !is.na(b)) out[!is.na(x) & x <= b] <- paste0("BP4_", l) }
  out
}
ord <- c("PP3_Strong", "PP3_Moderate", "PP3_Supporting", "Indeterminate", "BP4_Supporting", "BP4_Moderate", "BP4_Strong")

auc <- function(s, y) { r <- frank(s, ties.method = "average"); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(as.numeric(r[y == 1])) - n1 * (n1 + 1) / 2) / (as.numeric(n1) * n0) }

lr_table <- function(dd, tool) {
  iv <- interval_of(dd[[tool]], tool); y <- dd$abnormal; keep <- !is.na(iv); iv <- iv[keep]; y <- y[keep]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  base <- rbindlist(lapply(ord, function(o) data.table(interval = o, n_abnormal = sum(iv == o & y == 1), n_normal = sum(iv == o & y == 0))))
  base[, `:=`(prop_abnormal = n_abnormal / pmax(n_abnormal + n_normal, 1), LR = (n_abnormal / n1) / (n_normal / n0))]
  bt <- replicate(B, { i <- sample.int(length(y), replace = TRUE); yy <- y[i]; ii <- iv[i]; m1 <- sum(yy == 1); m0 <- sum(yy == 0)
    vapply(ord, function(o) (sum(ii == o & yy == 1) / m1) / (sum(ii == o & yy == 0) / m0), 0) })
  base[, `:=`(LR_lo = apply(bt, 1, quantile, 0.025, na.rm = TRUE), LR_hi = apply(bt, 1, quantile, 0.975, na.rm = TRUE))]
  base[, `:=`(tool = tool, n_snvs = length(y))][]
}

run_set <- function(dd, label, tools) {
  lr <- rbindlist(lapply(tools, function(t) lr_table(dd[!is.na(get(t))], t)))[, set := label]
  au <- rbindlist(lapply(tools, function(t) { x <- dd[!is.na(get(t))]
    bt <- replicate(B, { i <- sample.int(nrow(x), replace = TRUE); auc(x[[t]][i], x$abnormal[i]) })
    data.table(tool = t, n = nrow(x), AUROC = auc(x[[t]], x$abnormal), lo = quantile(bt, 0.025), hi = quantile(bt, 0.975)) }))[, set := label]
  list(lr = lr, auc = au)
}

primary <- d[!is.na(CADD) & !is.na(AVI)]
inter <- d[!is.na(CADD) & !is.na(AVI) & !is.na(AM) & !is.na(REVEL) & !is.na(BayesDel)]
clin <- readRDS("C:/cadd/data/refsets/clinvar_ref.rds")[, ..k]
indep <- primary[!clin, on = k]
cat("Primary set:", nrow(primary), "| five-tool intersection:", nrow(inter), "| primary without ClinVar reference SNVs:", nrow(indep), "\n")

r1 <- run_set(primary, "primary (CADD & AVI scored)", c("CADD", "AVI"))
r2 <- run_set(inter, "five-tool intersection", c("CADD", "AVI", "REVEL", "BayesDel", "AM"))
r3 <- run_set(indep, "sensitivity: ClinVar reference SNVs excluded", c("CADD", "AVI"))
lr_all <- rbind(r1$lr, r2$lr, r3$lr); auc_all <- rbind(r1$auc, r2$auc, r3$auc)
fwrite(lr_all, file.path(res_dir, "interval_LR.tsv"), sep = "\t"); fwrite(auc_all, file.path(res_dir, "auroc_pooled.tsv"), sep = "\t")

# paired AVI - CADD AUROC difference (primary set)
bt <- replicate(B, { i <- sample.int(nrow(primary), replace = TRUE); auc(primary$AVI[i], primary$abnormal[i]) - auc(primary$CADD[i], primary$abnormal[i]) })
pdiff <- data.table(set = "primary", diff_AVI_minus_CADD = auc(primary$AVI, primary$abnormal) - auc(primary$CADD, primary$abnormal),
                    lo = quantile(bt, 0.025), hi = quantile(bt, 0.975))
fwrite(pdiff, file.path(res_dir, "auroc_paired_difference.tsv"), sep = "\t")

# per gene
pg <- primary[, .(n_abnormal = sum(abnormal), n_normal = sum(abnormal == 0)), by = gene][n_abnormal >= 30 & n_normal >= 30]
per_gene <- rbindlist(lapply(pg$gene, function(g) { x <- primary[gene == g]; xi <- inter[gene == g]
  data.table(gene = g, n = nrow(x), n_abnormal = sum(x$abnormal), AUROC_CADD = auc(x$CADD, x$abnormal), AUROC_AVI = auc(x$AVI, x$abnormal),
             n_inter = nrow(xi),
             AUROC_REVEL = if (sum(xi$abnormal) >= 30 && sum(xi$abnormal == 0) >= 30) auc(xi$REVEL, xi$abnormal) else NA_real_,
             AUROC_BayesDel = if (sum(xi$abnormal) >= 30 && sum(xi$abnormal == 0) >= 30) auc(xi$BayesDel, xi$abnormal) else NA_real_,
             AUROC_AM = if (sum(xi$abnormal) >= 30 && sum(xi$abnormal == 0) >= 30) auc(xi$AM, xi$abnormal) else NA_real_) }))
fwrite(per_gene[order(gene)], file.path(res_dir, "auroc_per_gene.tsv"), sep = "\t")

options(width = 200)
cat("\n== Pooled AUROC ==\n"); print(auc_all)
cat("\n== Paired AVI - CADD ==\n"); print(pdiff)
cat("\n== Interval LRs ==\n"); print(lr_all[, .(set, tool, interval, n_abnormal, n_normal, prop_abnormal = round(prop_abnormal, 3), LR = round(LR, 2), LR_lo = round(LR_lo, 2), LR_hi = round(LR_hi, 2))], nrows = 200)
cat("\n== Per gene ==\n"); print(per_gene[order(gene)], digits = 3)
cat("\nGenes where AVI AUROC > CADD AUROC:", per_gene[AUROC_AVI > AUROC_CADD, .N], "of", nrow(per_gene), "\n")
