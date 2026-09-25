# Aim 6 (amendment A3): calibration of CADD v1.7, AVI, REVEL, AlphaMissense and BayesDel on the
# 2026 reference variants scored by all five predictors (intersection), 5 seeds each, B = 10,000;
# discrimination on the intersection and paired bootstrap differences versus AVI.
# Output: C:/cadd/calibration/results/v17_avi/comparators/<tool>/seed<seed>/ and .../comparators/summary/
# Usage: Rscript 21_calibrate_comparators.R [workers]
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/pejaver_port_fast.R")
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 10L
root <- "C:/cadd/calibration/results/v17_avi/comparators"; sdir <- file.path(root, "summary")
dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
seeds <- c(20260925L, 101L, 202L, 303L, 404L)
sc <- readRDS("C:/cadd/data/refsets/comparator_scores.rds")
sc[, `:=`(REVEL = round(REVEL, 3), AM = round(AM, 4), BayesDel = round(BayesDel, 2))]
tools <- list(CADD = list(col = "CADD_PHRED", inc = 0.1), AVI = list(col = "AVI_PHRED", inc = 0.1),
              REVEL = list(col = "REVEL", inc = 0.001), AlphaMissense = list(col = "AM", inc = 0.001),
              BayesDel = list(col = "BayesDel", inc = 0.01))
cols <- vapply(tools, `[[`, "", "col")
clin <- merge(readRDS("C:/cadd/data/refsets/clinvar_ref.rds"), sc, by = k)
gno  <- merge(readRDS("C:/cadd/data/refsets/gnomad_ref.rds"), sc, by = k)
clin <- clin[complete.cases(clin[, ..cols])]; gno <- gno[complete.cases(gno[, ..cols])]
cat(sprintf("Intersection: ClinVar %d P/LP, %d B/LB (%d genes); gnomAD %d\n",
            clin[y == 1, .N], clin[y == 0, .N], uniqueN(clin$gene), nrow(gno)))
fwrite(data.table(P_LP = clin[y == 1, .N], B_LB = clin[y == 0, .N], genes = uniqueN(clin$gene), gnomAD = nrow(gno)),
       file.path(sdir, "intersection_sizes.tsv"), sep = "\t")

for (tn in names(tools)) for (s in seeds) {
  out <- file.path(root, tn, paste0("seed", s))
  if (file.exists(file.path(out, "thresholds.tsv"))) { cat("skip", tn, s, "\n"); next }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  r <- calibrate_fast(clin[[tools[[tn]]$col]], clin$y, gno[[tools[[tn]]$col]], B = 10000, workers = workers,
                      seed = s, increment = tools[[tn]]$inc)
  fwrite(r$thresholds, file.path(out, "thresholds.tsv"), sep = "\t")
  fwrite(r$levels, file.path(out, "posterior_levels.tsv"), sep = "\t")
  fwrite(r$curves, file.path(out, "posterior_curves.tsv"), sep = "\t")
  fwrite(r$boot_thresholds, file.path(out, "bootstrap_thresholds.tsv"), sep = "\t")
  writeLines(capture.output(str(r$params)), file.path(out, "parameters.txt"))
  cat(sprintf("[%s] %s seed %d: %.1f min\n", format(Sys.time(), "%H:%M"), tn, s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## Summary of thresholds -------------------------------------------------------------------------
num <- c("PP3_point", "PP3_discounted", "BP4_point", "BP4_discounted", "PP3_boot_NA_frac", "BP4_boot_NA_frac")
all <- rbindlist(lapply(names(tools), function(tn) rbindlist(lapply(seeds, function(s)
  data.table(tool = tn, seed = s, fread(file.path(root, tn, paste0("seed", s), "thresholds.tsv"),
                                         colClasses = c(level = "character", setNames(rep("numeric", length(num)), num))))))))
fwrite(all, file.path(sdir, "thresholds_all_runs.tsv"), sep = "\t")
fmt <- function(v) if (all(is.na(v))) "not reached" else if (any(is.na(v))) sprintf("%s (reached in %d/%d seeds)", format(median(v, na.rm = TRUE)), sum(!is.na(v)), length(v)) else
  if (min(v) == max(v)) format(min(v)) else sprintf("%s [%s-%s]", format(median(v)), format(min(v)), format(max(v)))
all[, level := factor(level, levels = c("Supporting", "Moderate", "Strong", "VeryStrong"))]
tab <- all[, .(PP3 = fmt(PP3_discounted), BP4 = fmt(BP4_discounted)), by = .(tool, level)][order(tool, level)]
cat("\n== Thresholds on the intersection (median [range] over 5 seeds) ==\n"); print(tab)
fwrite(tab, file.path(sdir, "table_thresholds.tsv"), sep = "\t")

## Evidence levels assigned to the reference variants (median thresholds; in-sample) --------------
med <- all[, .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)]
assign_level <- function(s, m) {
  lv <- rep("Indeterminate", length(s))
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { t <- m[level == L, PP3]; if (length(t) && !is.na(t)) lv[s >= t] <- paste0("PP3_", L) }
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { t <- m[level == L, BP4]; if (length(t) && !is.na(t)) lv[s <= t] <- paste0("BP4_", L) }
  lv
}
ev <- rbindlist(lapply(names(tools), function(tn) data.table(tool = tn, class = ifelse(clin$y == 1, "P/LP", "B/LB"),
                                                           evidence = assign_level(clin[[tools[[tn]]$col]], med[tool == tn]))))
evt <- dcast(ev[, .N, by = .(tool, class, evidence)][, pct := round(100 * N / sum(N), 1), by = .(tool, class)],
             tool + evidence ~ class, value.var = "pct", fill = 0)
cat("\n== % of intersection reference variants per evidence level ==\n"); print(evt[order(tool, evidence)])
fwrite(evt, file.path(sdir, "table_evidence_levels.tsv"), sep = "\t")

## Discrimination and paired differences versus AVI -----------------------------------------------
auroc <- function(s, y) { r <- rank(s); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0)); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
auprc <- function(s, y) {
  d <- data.table(s, y)[, .(tp = sum(y == 1), fp = sum(y == 0)), by = s][order(-s)]
  d[, `:=`(TP = cumsum(tp), FP = cumsum(fp))]; d[, `:=`(prec = TP / (TP + FP), rec = TP / sum(tp))]
  sum(diff(c(0, d$rec)) * d$prec)
}
y <- clin$y; S <- lapply(tools, function(t) clin[[t$col]])
est <- rbindlist(lapply(names(S), function(tn) data.table(tool = tn, AUROC = auroc(S[[tn]], y), AUPRC = auprc(S[[tn]], y))))
set.seed(20260925)
bs <- rbindlist(lapply(1:1000, function(b) { i <- sample.int(length(y), replace = TRUE)
  rbindlist(lapply(names(S), function(tn) data.table(b = b, tool = tn, AUROC = auroc(S[[tn]][i], y[i]), AUPRC = auprc(S[[tn]][i], y[i])))) }))
ci <- bs[, .(AUROC_lo = quantile(AUROC, 0.025), AUROC_hi = quantile(AUROC, 0.975),
             AUPRC_lo = quantile(AUPRC, 0.025), AUPRC_hi = quantile(AUPRC, 0.975)), by = tool]
d <- merge(bs, bs[tool == "AVI", .(b, AUROC_avi = AUROC, AUPRC_avi = AUPRC)], by = "b")[tool != "AVI"]
pair <- d[, .(dAUROC_vs_AVI = NA_real_, dAUROC_lo = quantile(AUROC - AUROC_avi, 0.025), dAUROC_hi = quantile(AUROC - AUROC_avi, 0.975),
              n_dAUROC_ge0 = sum(AUROC - AUROC_avi >= 0),
              dAUPRC_lo = quantile(AUPRC - AUPRC_avi, 0.025), dAUPRC_hi = quantile(AUPRC - AUPRC_avi, 0.975),
              n_dAUPRC_ge0 = sum(AUPRC - AUPRC_avi >= 0)), by = tool]
disc <- merge(merge(est, ci, by = "tool"), pair, by = "tool", all.x = TRUE)
disc[, dAUROC_vs_AVI := AUROC - est[tool == "AVI", AUROC]][, dAUPRC_vs_AVI := AUPRC - est[tool == "AVI", AUPRC]]
cat("\n== Discrimination on the intersection (95% bootstrap CI; paired differences vs AVI, 1,000 replicates) ==\n")
print(disc[order(-AUROC)], digits = 4)
fwrite(disc, file.path(sdir, "table_discrimination.tsv"), sep = "\t")
