# Aim 8 (amendment A4): population frequencies of MANE missense SNVs that reach PP3 levels with
# CADD v1.7 or AVI (main-analysis median thresholds), across gnomAD v4.1 genetic ancestry groups.
# Output: C:/cadd/calibration/results/v17_avi/ancestry/
suppressPackageStartupMessages(library(data.table))
cdir <- "C:/cadd/data/coding"; out <- "C:/cadd/calibration/results/v17_avi/ancestry"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
groups <- c("AFR", "AMR", "ASJ", "EAS", "FIN", "MID", "NFE", "SAS", "AMI")
thr <- fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main"]
thr <- thr[, .(PP3 = median(PP3_discounted)), by = .(tool, level)]
thrs <- list(CADD = c(Supporting = thr[tool == "CADD" & level == "Supporting", PP3], Moderate = thr[tool == "CADD" & level == "Moderate", PP3]),
             AVI  = c(Supporting = thr[tool == "AVI" & level == "Supporting", PP3],  Moderate = thr[tool == "AVI" & level == "Moderate", PP3]))
print(thrs)
lvl <- function(s, t) fifelse(s >= t["Moderate"], "Moderate", fifelse(s >= t["Supporting"], "Supporting", "none"))

af_all <- readRDS(file.path(cdir, "gnomad_af_mane_coding.rds"))
n_obs_missense <- 0L
pp3 <- rbindlist(lapply(list.files(cdir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f) {
  d <- readRDS(f)[, .(chrom, pos, ref, alt, gene_name, consequence)]
  ch <- d$chrom[1]
  bad <- unique(d[consequence %chin% c("stop_gained", "stop_lost", "start_lost", "incomplete_terminal_codon_variant"), .(chrom, pos, alt)])
  d <- d[consequence == "missense"][!bad, on = .(chrom, pos, alt)]
  v <- d[, .(gene = paste(sort(unique(gene_name)), collapse = "&")), by = .(chrom, pos, ref, alt)]
  rm(d)
  a <- af_all[chrom == ch]
  v <- merge(v, a, by = c("chrom", "pos", "ref", "alt"))
  v <- v[!is.na(AF_max) & AF_max > 0]
  n_obs_missense <<- n_obs_missense + nrow(v)
  for (g in groups) {
    e <- paste0("gnomADe_", g); gg <- paste0("gnomADg_", g)
    ev <- if (e %in% names(v)) v[[e]] else NA_real_; gv <- if (gg %in% names(v)) v[[gg]] else NA_real_
    set(v, j = paste0("AF_", g), value = pmax(ev, gv, na.rm = TRUE))
  }
  s1 <- readRDS(file.path(cdir, sprintf("cadd_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, CADD = PHRED)]
  s2 <- readRDS(file.path(cdir, sprintf("avi_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, AVI = round(AVI_PHRED, 2))]
  v <- merge(merge(v, s1, by = c("chrom", "pos", "ref", "alt")), s2, by = c("chrom", "pos", "ref", "alt"))
  v[, `:=`(PP3_CADD = lvl(CADD, thrs$CADD), PP3_AVI = lvl(AVI, thrs$AVI))]
  v <- v[PP3_CADD != "none" | PP3_AVI != "none", c("chrom", "pos", "ref", "alt", "gene", "AF_max", paste0("AF_", groups),
                                                     "CADD", "AVI", "PP3_CADD", "PP3_AVI"), with = FALSE]
  cat(ch, nrow(v), "\n"); v
}))
rm(af_all); invisible(gc())
afc <- paste0("AF_", groups)
pp3[, grp_AFmax := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = afc]
pp3[, grpmax := groups[max.col(replace(as.matrix(.SD), is.na(as.matrix(.SD)), -1), ties.method = "first")], .SDcols = afc]
pp3[, n_groups_gt_0.01 := rowSums(as.matrix(.SD) > 0.01, na.rm = TRUE), .SDcols = afc]
saveRDS(pp3, file.path(out, "pp3_level_missense_observed_in_gnomad.rds"))
cat("MANE missense SNVs observed in gnomAD v4.1:", n_obs_missense, "\n")

wilson <- function(x, n) { if (n == 0) return(c(NA, NA)); z <- 1.96; p <- x / n
  c((p + z^2/(2*n) - z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n), (p + z^2/(2*n) + z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)) }

res1 <- list(); res2 <- list(); res3 <- list()
for (tool in c("CADD", "AVI")) for (L in c("Supporting+", "Moderate")) {
  sel <- if (L == "Moderate") pp3[[paste0("PP3_", tool)]] == "Moderate" else pp3[[paste0("PP3_", tool)]] != "none"
  x <- pp3[sel]
  n <- nrow(x)
  res1[[paste(tool, L)]] <- data.table(tool = tool, PP3 = L, n_observed = n,
    n_grpAF_gt_0.01 = x[grp_AFmax > 0.01, .N], pct_gt_0.01 = 100 * x[grp_AFmax > 0.01, .N] / n,
    n_grpAF_gt_0.05 = x[grp_AFmax > 0.05, .N], pct_gt_0.05 = 100 * x[grp_AFmax > 0.05, .N] / n)
  res2[[paste(tool, L)]] <- x[grp_AFmax > 0.01, .N, by = grpmax][, `:=`(tool = tool, PP3 = L, pct = round(100 * N / sum(N), 1))]
  for (g in groups) {
    a <- x[[paste0("AF_", g)]]; obs <- sum(!is.na(a) & a > 0); c1 <- sum(!is.na(a) & a > 0.01); c5 <- sum(!is.na(a) & a > 0.05)
    only <- sum(!is.na(a) & a > 0.01 & x$n_groups_gt_0.01 == 1)
    w1 <- wilson(c1, obs); w5 <- wilson(c5, obs)
    res3[[paste(tool, L, g)]] <- data.table(tool = tool, PP3 = L, group = g, n_observed_in_group = obs,
      n_gt_0.01 = c1, pct_gt_0.01 = 100 * c1 / obs, lo_0.01 = 100 * w1[1], hi_0.01 = 100 * w1[2],
      n_gt_0.05 = c5, pct_gt_0.05 = 100 * c5 / obs, lo_0.05 = 100 * w5[1], hi_0.05 = 100 * w5[2],
      n_group_specific_gt_0.01 = only)
  }
}
t1 <- rbindlist(res1); t2 <- rbindlist(res2); t3 <- rbindlist(res3)
cat("\n== PP3-level missense observed in gnomAD with group AF > 1% / > 5% ==\n"); print(t1, digits = 3)
cat("\n== grpmax among PP3-level variants with group AF > 1% ==\n"); print(t2[order(tool, PP3, -N)])
cat("\n== Per group: % of PP3-level variants observed in the group that exceed 1% / 5% there ==\n")
print(t3[, .(tool, PP3, group, n_observed_in_group, pct_gt_0.01 = round(pct_gt_0.01, 3), CI_0.01 = sprintf("%.3f-%.3f", lo_0.01, hi_0.01),
             pct_gt_0.05 = round(pct_gt_0.05, 3), n_group_specific_gt_0.01)])
fwrite(t1, file.path(out, "table_common_pp3.tsv"), sep = "\t"); fwrite(t2, file.path(out, "table_grpmax.tsv"), sep = "\t")
fwrite(t3, file.path(out, "table_per_group.tsv"), sep = "\t")
