# Main article, amendment A17 (Aim 9 extension): functional check, independent of ClinVar, of the main CADD v1.7 and
# AVI thresholds outside missense variants.
# Data: (1) SGE, the 60,382 SNVs of the 12 saturation genome editing sets built for Paper 2 (script 62; abnormal vs
#       normal), by region class; classes with >= 30 abnormal and >= 30 normal SNVs; missense as a control.
#       (2) MPRA, record-level set of R13 (mpra_scored.rds; rule of A10), pooled and per element (descriptive).
# Scores: CADD v1.7 PHRED and AVI Phred (2 decimals) from the local files of the main article.
# Thresholds: main thresholds (median of the five seeds, discounted) in every class; splice-region thresholds (A10,
#       script 30b) in addition in the intronic splice region class.
# Metrics: interval LR with 0.5 added to zero counts, 95% CI from RB bootstrap resamples of SNVs (seed 20261101),
#       intervals with >= 10 SNVs; direction check; proportion of abnormal and normal SNVs with BP4; AUROC of CADD and
#       AVI with the paired bootstrap difference. Sensitivity: without SNVs of the ClinVar reference sets (R1, R13).
# Usage: Rscript 65_functional_noncoding_main_thresholds.R [RB = 1000]
# Output: C:/cadd/calibration/results/v17_avi/functional_noncoding/ (or functional_noncoding_test_RB<RB> if RB < 1000)
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
a <- commandArgs(TRUE); RB <- if (length(a)) as.integer(a[1]) else 1000L
out <- if (RB < 1000L) sprintf("C:/cadd/calibration/results/v17_avi/functional_noncoding_test_RB%d", RB) else
  "C:/cadd/calibration/results/v17_avi/functional_noncoding"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt"); SEED <- 20261101L
req <- c(Supporting = 2.406, Moderate = 5.790, Strong = 33.53)
levs <- c("PP3_Strong", "PP3_Moderate", "PP3_Supporting", "Indeterminate", "BP4_Supporting", "BP4_Moderate", "BP4_Strong")

## Score lookup (function of script 30) -------------------------------------------------------------
lookup_chr <- function(file, pref, x, col, round2 = FALSE) {
  p <- sort(unique(x$pos)); grp <- cumsum(c(1L, diff(p) > 1L)); rg <- data.table(s = p, g = grp)[, .(start = min(s), end = max(s)), by = g]
  tf <- open(TabixFile(file)); on.exit(close(tf))
  r <- rbindlist(lapply(split(rg, ceiling(seq_len(nrow(rg)) / 5000)), function(b) {
    l <- unlist(scanTabix(tf, param = GenomicRanges::GRanges(paste0(pref, x$chrom[1]), IRanges::IRanges(b$start, b$end))), use.names = FALSE)
    if (!length(l)) return(NULL); fread(text = l, sep = "\t", header = FALSE, select = c(2, 3, 4, 6), col.names = c("pos", "ref", "alt", "sc"), nThread = 1) }))
  if (!nrow(r)) return(NULL)
  if (round2) r[, sc := round(sc, 2)]; r[, chrom := x$chrom[1]]; setnames(r, "sc", col); r[x[, ..k], on = k, nomatch = NULL] }
fc <- "C:/cadd/data/whole_genome_SNVs.tsv.gz"; fa <- "C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz"
sn_cadd <- seqnamesTabix(TabixFile(fc)); sn_avi <- sub("^chr", "", seqnamesTabix(TabixFile(fa)))

sge <- unique(readRDS("C:/cadd/calibration/results/paper2/functional/sge_labels.rds"), by = k)
u <- unique(sge[, ..k])
s1 <- rbindlist(lapply(split(u[chrom %chin% sn_cadd], by = "chrom"), function(x) lookup_chr(fc, "", x, "CADD_PHRED")))
s2 <- rbindlist(lapply(split(u[chrom %chin% sn_avi], by = "chrom"), function(x) lookup_chr(fa, "chr", x, "AVI_PHRED", TRUE)))
sge <- merge(merge(sge, unique(s1, by = k), by = k, all.x = TRUE), unique(s2, by = k), by = k, all.x = TRUE)
cat(sprintf("SGE: %d SNVs; with CADD %d; with AVI %d; with both %d\n", nrow(sge), sum(!is.na(sge$CADD_PHRED)),
            sum(!is.na(sge$AVI_PHRED)), sum(!is.na(sge$CADD_PHRED) & !is.na(sge$AVI_PHRED))))
sge <- sge[!is.na(CADD_PHRED) & !is.na(AVI_PHRED)]
mp <- readRDS("C:/cadd/calibration/results/v17_avi/noncoding/mpra_scored.rds")[, y := functional]
cat(sprintf("MPRA: %d SNV-element records (%d functional)\n", nrow(mp), sum(mp$y)))

## Thresholds ---------------------------------------------------------------------------------------
thr_main <- fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main",
              .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)]
thr_sr <- fread("C:/cadd/calibration/results/v17_avi/noncoding/splice_region_intronic/thresholds.tsv")[,
              .(tool, level, PP3 = PP3_discounted, BP4 = BP4_discounted)]
cat("\nMain thresholds:\n"); print(thr_main); cat("Splice-region thresholds:\n"); print(thr_sr)
ivl <- function(s, t) { lv <- rep("Indeterminate", length(s))
  for (L in c("Supporting", "Moderate", "Strong")) { v <- t[level == L, PP3]; if (length(v) && !is.na(v)) lv[s >= v] <- paste0("PP3_", L) }
  for (L in c("Supporting", "Moderate", "Strong")) { v <- t[level == L, BP4]; if (length(v) && !is.na(v)) lv[s <= v] <- paste0("BP4_", L) }; lv }

## Metrics ------------------------------------------------------------------------------------------
lrv <- function(a, n, nA, nN) { a <- ifelse(a == 0, 0.5, a); n <- ifelse(n == 0, 0.5, n); (a / nA) / (n / nN) }
auroc <- function(s, y) { r <- rank(s); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0)); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
eval_group <- function(d, group, tset, t, tl) {
  col <- if (tl == "CADD") "CADD_PHRED" else "AVI_PHRED"
  iv <- match(ivl(d[[col]], t[tool == tl]), levs); y <- d$y; nA <- sum(y == 1); nN <- sum(y == 0)
  A <- tabulate(iv[y == 1], 7); N <- tabulate(iv[y == 0], 7); lr <- lrv(A, N, nA, nN)
  set.seed(SEED); i1 <- which(y == 1); i0 <- which(y == 0)
  bs <- vapply(seq_len(RB), function(b) { j1 <- iv[i1[sample.int(nA, nA, TRUE)]]; j0 <- iv[i0[sample.int(nN, nN, TRUE)]]
    lrv(tabulate(j1, 7), tabulate(j0, 7), nA, nN) }, numeric(7))
  r <- data.table(group = group, thresholds = tset, tool = tl, iv = levs, a = A, n = N, nA = nA, nN = nN, LR = lr,
                  lo = apply(bs, 1, quantile, 0.025), hi = apply(bs, 1, quantile, 0.975))
  r[, `:=`(dir = sub("_.*", "", iv), level = sub(".*_", "", iv))]
  r[, evaluable := iv != "Indeterminate" & a + n >= 10]
  r[, direction_ok := fifelse(dir == "PP3", LR > 1, fifelse(dir == "BP4", LR < 1, NA))]
  r[, meets_clingen := fifelse(dir == "PP3", LR >= req[level], fifelse(dir == "BP4", LR <= 1 / req[level], NA))]
  r[iv == "Indeterminate", c("direction_ok", "meets_clingen") := NA]
  r[, c("dir", "level") := NULL][] }
bp4_rates <- function(d, group, tset, t) rbindlist(lapply(c("CADD", "AVI"), function(tl) {
  col <- if (tl == "CADD") "CADD_PHRED" else "AVI_PHRED"; iv <- ivl(d[[col]], t[tool == tl])
  data.table(group = group, thresholds = tset, tool = tl, nA = sum(d$y == 1), nN = sum(d$y == 0),
             abnormal_BP4 = mean(grepl("^BP4", iv[d$y == 1])), normal_BP4 = mean(grepl("^BP4", iv[d$y == 0])),
             abnormal_PP3 = mean(grepl("^PP3", iv[d$y == 1])), normal_PP3 = mean(grepl("^PP3", iv[d$y == 0]))) }))
auc_group <- function(d, group) {
  a1 <- auroc(d$CADD_PHRED, d$y); a2 <- auroc(d$AVI_PHRED, d$y); set.seed(SEED)
  bs <- replicate(RB, { i <- sample.int(nrow(d), replace = TRUE); auroc(d$AVI_PHRED[i], d$y[i]) - auroc(d$CADD_PHRED[i], d$y[i]) })
  data.table(group = group, n = nrow(d), nA = sum(d$y == 1), CADD_AUROC = a1, AVI_AUROC = a2, diff_AVI_CADD = a2 - a1,
             lo = quantile(bs, 0.025), hi = quantile(bs, 0.975)) }

run_set <- function(tag, sg, m) {
  cls <- sg[, .(nA = sum(y == 1), nN = sum(y == 0)), by = class][nA >= 30 & nN >= 30, class]
  cat(sprintf("\n== %s: SGE %d SNVs; evaluable classes: %s ==\n", tag, nrow(sg), paste(cls, collapse = ", ")))
  jobs <- c(lapply(cls, function(cl) list(d = sg[class == cl], g = cl, ts = "main", t = thr_main)),
            if ("splice_region_intronic" %chin% cls) list(list(d = sg[class == "splice_region_intronic"], g = "splice_region_intronic", ts = "splice_region", t = thr_sr)))
  iv <- rbindlist(lapply(jobs, function(j) rbind(eval_group(j$d, j$g, j$ts, j$t, "CADD"), eval_group(j$d, j$g, j$ts, j$t, "AVI"))))
  br <- rbindlist(lapply(jobs, function(j) bp4_rates(j$d, j$g, j$ts, j$t)))
  au <- rbindlist(lapply(cls, function(cl) auc_group(sg[class == cl], cl)))
  fwrite(iv, file.path(out, sprintf("sge_intervals%s.tsv", tag)), sep = "\t"); fwrite(br, file.path(out, sprintf("sge_bp4_pp3_rates%s.tsv", tag)), sep = "\t")
  fwrite(au, file.path(out, sprintf("sge_auroc%s.tsv", tag)), sep = "\t")
  print(iv[evaluable == TRUE, .(group, thresholds, tool, iv, a, n, LR = signif(LR, 3), lo = signif(lo, 3), hi = signif(hi, 3), direction_ok, meets_clingen)])
  print(br); print(au)
  # MPRA, pooled over elements (bootstrap) and per element (point estimates)
  miv <- rbind(eval_group(m, "MPRA_all", "main", thr_main, "CADD"), eval_group(m, "MPRA_all", "main", thr_main, "AVI"))
  mbr <- bp4_rates(m, "MPRA_all", "main", thr_main); mau <- auc_group(m, "MPRA_all")
  mel <- rbindlist(lapply(split(m, by = "Element"), function(x) rbindlist(lapply(c("CADD", "AVI"), function(tl) {
    col <- if (tl == "CADD") "CADD_PHRED" else "AVI_PHRED"; iv <- match(ivl(x[[col]], thr_main[tool == tl]), levs)
    A <- tabulate(iv[x$y == 1], 7); N <- tabulate(iv[x$y == 0], 7)
    data.table(element = x$Element[1], tool = tl, iv = levs, a = A, n = N, nA = sum(x$y == 1), nN = sum(x$y == 0),
               LR = lrv(A, N, sum(x$y == 1), sum(x$y == 0))) }))))
  fwrite(miv, file.path(out, sprintf("mpra_intervals%s.tsv", tag)), sep = "\t"); fwrite(rbind(mbr), file.path(out, sprintf("mpra_bp4_pp3_rates%s.tsv", tag)), sep = "\t")
  fwrite(mau, file.path(out, sprintf("mpra_auroc%s.tsv", tag)), sep = "\t"); fwrite(mel, file.path(out, sprintf("mpra_intervals_by_element%s.tsv", tag)), sep = "\t")
  cat(sprintf("\n-- MPRA %s: %d records (%d functional) --\n", tag, nrow(m), sum(m$y)))
  print(miv[, .(tool, iv, a, n, LR = signif(LR, 3), lo = signif(lo, 3), hi = signif(hi, 3), evaluable, direction_ok, meets_clingen)]); print(mbr); print(mau)
  e <- mel[a + n >= 10 & iv != "Indeterminate" & nA >= 10 & nN >= 10]
  cat("Per element (intervals with >= 10 records; elements with >= 10 functional and non-functional): BP4 intervals with LR < 1\n")
  print(e[grepl("^BP4", iv), .(n_intervals = .N, LR_below_1 = sum(LR < 1)), by = .(tool, iv)])
}

run_set("", sge, mp)
ref <- unique(rbind(readRDS("C:/cadd/data/refsets/clinvar_ref.rds")[, ..k],
                    readRDS("C:/cadd/calibration/results/v17_avi/noncoding/clinvar_noncoding_reference.rds")[, ..k]))
run_set("_sens_no_clinvar_ref", sge[!ref, on = k], mp[!ref, on = k])
cat("\nFinished", format(Sys.time()), "\n")
