# Aim 2c - temporal hold-out validation (plan amendment A2).
#  1. Rebuild the 2026-09-23 reference sets from the all-MANE gnomAD table (13a) and check they are
#     identical to the sets used in the main calibration (consistency of the criteria code).
#  2. Build the 2024-01-27 reference sets with the same criteria.
#  3. Calibrate CADD v1.7 PHRED and AVI Phred on the 2024 sets (B = 10,000, seed 20260925).
#  4. Test set = 2026 main-set variants absent from the 2024 set (label changes excluded, counted).
#  5. Interval-based positive likelihood ratios on the test set (1,000 bootstrap replicates),
#     compared with the LR required for each evidence strength.
# Output: C:/cadd/calibration/results/v17_avi/temporal_holdout/
# Usage: Rscript 18_temporal_holdout.R [workers]
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/pejaver_port_fast.R")
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 8L
cdir <- "C:/cadd/data/coding"; out <- "C:/cadd/calibration/results/v17_avi/temporal_holdout"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat(sprintf(...), "\n")
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")

## ClinVar reader (same rules as 13_build_reference_sets.R) --------------------------------
read_clinvar <- function(f) {
  cv <- fread(f, sep = "\t", skip = "#CHROM", select = c(1, 2, 3, 4, 5, 8),
              col.names = c("chrom", "pos", "clinvar_id", "ref", "alt", "info"), colClasses = list(character = c(1, 3)))
  getk <- function(x, k) { v <- sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", k, "="), x)] <- NA; v }
  cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
  cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                      grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                      grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
  dg <- cv[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo),
           unique(unlist(lapply(strsplit(geneinfo, "|", fixed = TRUE), function(z) sub(":.*", "", z))))]
  list(dg = dg, snv = cv[nchar(ref) == 1 & nchar(alt) == 1, .(chrom, pos, ref, alt, clinvar_id, clnsig, stars)], n = nrow(cv))
}
cv26 <- read_clinvar("C:/cadd/annot/clinvar.vcf.gz")
cv24 <- read_clinvar("C:/cadd/annot/clinvar_archive/clinvar_20240127.vcf.gz")
log("ClinVar 2026-09-23: %d records, %d disease genes | 2024-01-27: %d records, %d disease genes",
    cv26$n, length(cv26$dg), cv24$n, length(cv24$dg))
dgu <- union(cv26$dg, cv24$dg)

## MANE missense SNVs in disease genes of either release + exon-edge flag (as in 13) --------------
pos <- readRDS(file.path(cdir, "mane_cds_positions.rds"))
setorder(pos, transcript_id, cds_pos)
pos[, jump := c(abs(diff(pos)) != 1L, FALSE), by = transcript_id][, idx := seq_len(.N), by = transcript_id]
jn <- pos[jump == TRUE, .(transcript_id, k = idx)]
flag <- rbindlist(lapply(c(-2:0, 1:3), function(o) jn[, .(transcript_id, idx = k + o)]))
pos[, near := FALSE][flag, on = .(transcript_id, idx), near := TRUE]
edge <- unique(pos[near == TRUE, .(chrom, pos)])[, exon_edge := TRUE]
rm(pos, jn, flag); invisible(gc())
# Memory-lean (v2): disease-gene membership for each release is flagged per chromosome from the
# per-transcript gene names (no string splitting over the whole table); only needed columns kept.
mis <- rbindlist(lapply(list.files(cdir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f) {
  d <- readRDS(f)[, .(chrom, pos, ref, alt, gene_name, consequence)]
  bad <- unique(d[consequence %chin% c("stop_gained", "stop_lost", "start_lost", "incomplete_terminal_codon_variant"), .(chrom, pos, alt)])
  d <- d[consequence == "missense"][!bad, on = .(chrom, pos, alt)][, consequence := NULL]
  d <- d[d[gene_name %chin% dgu, unique(.SD), .SDcols = c("chrom", "pos", "alt")], on = .(chrom, pos, alt)]
  r <- d[, .(gene = paste(sort(unique(gene_name)), collapse = "&"),
             dg26 = any(gene_name %chin% cv26$dg), dg24 = any(gene_name %chin% cv24$dg)), by = .(chrom, pos, ref, alt)]
  rm(d); invisible(gc()); r
}))
mis[edge, on = .(chrom, pos), exon_edge := i.exon_edge][is.na(exon_edge), exon_edge := FALSE]
af <- readRDS(file.path(cdir, "gnomad_af_mane_coding.rds"))
af <- af[, .(chrom, pos, ref, alt, AF_max)]; invisible(gc())
mis[af, on = .(chrom, pos, ref, alt), AF_max := i.AF_max][is.na(AF_max), AF_max := 0]
rm(af, edge); invisible(gc())
log("missense SNVs in disease genes of either release: %s", format(nrow(mis), big.mark = ","))

attach_scores <- function(d) rbindlist(lapply(split(d, d$chrom), function(x) {
  ch <- x$chrom[1]
  s1 <- readRDS(file.path(cdir, sprintf("cadd_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, CADD_PHRED = PHRED)]
  s2 <- readRDS(file.path(cdir, sprintf("avi_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, AVI_PHRED = round(AVI_PHRED, 2))]
  merge(merge(x, s1, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE), s2, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE)
}))
build_sets <- function(cv, flag) {
  m <- mis[get(flag) == TRUE]
  cvm <- merge(m, cv$snv, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE)
  clin <- cvm[!is.na(clnsig) & (clnsig %chin% plp | clnsig %chin% blb) & stars >= 1 & AF_max < 0.01][, y := as.integer(clnsig %chin% plp)]
  gno <- cvm[is.na(clinvar_id) & AF_max > 0 & AF_max < 0.01][!cv$snv, on = .(chrom, pos, ref, alt)]
  list(clin = attach_scores(clin), gno = attach_scores(gno))
}

## 1. Consistency check against the main-analysis sets -------------------------------------------
s26 <- build_sets(cv26, "dg26")
old_c <- readRDS("C:/cadd/data/refsets/clinvar_ref.rds"); old_g <- readRDS("C:/cadd/data/refsets/gnomad_ref.rds")
k <- c("chrom", "pos", "ref", "alt")
cmp <- function(a, b, cols) { a <- a[order(chrom, pos, alt), ..cols]; b <- b[order(chrom, pos, alt), ..cols]
  isTRUE(all.equal(a, b, check.attributes = FALSE)) }
log("Consistency 2026: ClinVar set identical (keys, y, stars, AF_max, scores): %s | gnomAD set identical (keys, AF_max, scores): %s",
    cmp(s26$clin, old_c, c(k, "y", "stars", "AF_max", "CADD_PHRED", "AVI_PHRED")),
    cmp(s26$gno, old_g, c(k, "AF_max", "CADD_PHRED", "AVI_PHRED")))
log("  sizes rebuilt: %d / %d ClinVar, %d / %d gnomAD", nrow(s26$clin), nrow(old_c), nrow(s26$gno), nrow(old_g))

## 2. 2024 sets and test set ---------------------------------------------------------------------
rm(s26); invisible(gc())
s24 <- build_sets(cv24, "dg24")
log("2024-01-27 sets: %d P/LP, %d B/LB (%d genes); gnomAD %d", s24$clin[y == 1, .N], s24$clin[y == 0, .N],
    uniqueN(s24$clin$gene), nrow(s24$gno))
both <- merge(old_c[, c(k, "y"), with = FALSE], s24$clin[, c(k, "y"), with = FALSE], by = k, suffixes = c("_26", "_24"))
changed <- both[y_26 != y_24]
test <- old_c[!s24$clin, on = k]
log("Test set (in 2026 main set, absent from 2024 set): %d P/LP, %d B/LB | label changed between releases (excluded): %d (P->B %d, B->P %d)",
    test[y == 1, .N], test[y == 0, .N], nrow(changed), changed[y_24 == 1 & y_26 == 0, .N], changed[y_24 == 0 & y_26 == 1, .N])
saveRDS(s24, file.path(out, "refsets_2024.rds")); saveRDS(test, file.path(out, "test_set.rds"))

## 3-5. Calibrate on 2024, evaluate on the test set ----------------------------------------------
req <- posterior_thresholds()$LR; names(req) <- c("VeryStrong", "Strong", "Moderate", "Supporting")
tools <- c(CADD = "CADD_PHRED", AVI = "AVI_PHRED")
res_all <- list()
for (tn in names(tools)) {
  x <- s24$clin[[tools[[tn]]]]; y <- s24$clin$y; g <- s24$gno[[tools[[tn]]]]
  t0 <- Sys.time()
  r <- calibrate_fast(x, y, g, B = 10000, workers = workers, seed = 20260925)
  log("[%s] calibrated on 2024: %.1f min", tn, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  th <- as.data.table(r$thresholds); th[, tool := tn]
  fwrite(th, file.path(out, sprintf("thresholds_2024_%s.tsv", tn)), sep = "\t")
  fwrite(as.data.table(r$curves), file.path(out, sprintf("posterior_curves_2024_%s.tsv", tn)), sep = "\t")
  fwrite(as.data.table(r$levels), file.path(out, sprintf("posterior_levels_2024_%s.tsv", tn)), sep = "\t")
  P <- setNames(th$PP3_discounted, th$level); Bn <- setNames(th$BP4_discounted, th$level)
  s <- test[[tools[[tn]]]]; yt <- test$y
  interval <- function(s) {
    lv <- rep("Indeterminate", length(s))
    for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) {
      if (!is.na(P[L])) lv[s >= P[L]] <- paste0("PP3_", L)
    }
    for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) {
      if (!is.na(Bn[L])) lv[s <= Bn[L]] <- paste0("BP4_", L)
    }
    lv
  }
  iv <- interval(s)
  lr <- function(iv, yt) { levs <- sort(unique(iv))
    vapply(levs, function(L) (mean(iv[yt == 1] == L)) / (mean(iv[yt == 0] == L)), 0) }
  est <- lr(iv, yt)
  set.seed(20260925)
  bs <- replicate(1000, { i <- sample.int(length(yt), replace = TRUE); v <- lr(iv[i], yt[i]); v[names(est)] })
  tab <- data.table(tool = tn, interval = names(est), n_PLP = as.integer(table(factor(iv[yt == 1], levels = names(est)))),
                    n_BLB = as.integer(table(factor(iv[yt == 0], levels = names(est)))), LR = est,
                    lower95 = apply(bs, 1, quantile, 0.025, na.rm = TRUE), upper95 = apply(bs, 1, quantile, 0.975, na.rm = TRUE))
  tab[, level := sub("^(PP3|BP4)_", "", interval)]
  tab[grepl("^PP3", interval), `:=`(LR_required = req[level], criterion = "LR >= required")]
  tab[grepl("^BP4", interval), `:=`(LR_required = 1 / req[level], criterion = "LR <= 1/required")]
  tab[, meets_point := fifelse(grepl("^PP3", interval), LR >= LR_required, fifelse(grepl("^BP4", interval), LR <= LR_required, NA))]
  tab[, meets_bound := fifelse(grepl("^PP3", interval), lower95 >= LR_required, fifelse(grepl("^BP4", interval), upper95 <= LR_required, NA))]
  res_all[[tn]] <- tab
}
res <- rbindlist(res_all)
print(res[, .(tool, interval, n_PLP, n_BLB, LR = signif(LR, 3), lower95 = signif(lower95, 3), upper95 = signif(upper95, 3),
              LR_required = signif(LR_required, 3), meets_point, meets_bound)])
fwrite(res, file.path(out, "interval_LR_test_set.tsv"), sep = "\t")

main <- fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main" & seed == 20260925]
drift <- rbindlist(lapply(names(tools), function(tn) merge(
  fread(file.path(out, sprintf("thresholds_2024_%s.tsv", tn)))[, .(tool, level, PP3_2024 = PP3_discounted, BP4_2024 = BP4_discounted)],
  main[tool == tn, .(tool, level, PP3_2026 = PP3_discounted, BP4_2026 = BP4_discounted)], by = c("tool", "level"))))
cat("\nThresholds calibrated on ClinVar 2024-01-27 vs 2026-09-23 (seed 20260925):\n"); print(drift)
fwrite(drift, file.path(out, "thresholds_2024_vs_2026.tsv"), sep = "\t")
