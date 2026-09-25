# Paired comparison of discrimination between AVI and CADD v1.7 on the main reference set.
# ADDED AFTER analysis plan v1 (2026-09-25): the plan specified AUROC/AUPRC per predictor with
# separate bootstrap intervals; this analysis adds the paired difference. Each of 1,000 bootstrap
# replicates resamples variants once and evaluates both predictors on the same sample.
# Output: calibration/results/v17_avi/summary/table_discrimination_paired.tsv
suppressPackageStartupMessages(library(data.table))
clin <- readRDS("C:/cadd/data/refsets/clinvar_ref.rds")
y <- clin$y; a <- clin$AVI_PHRED; c <- clin$CADD_PHRED
auroc <- function(s, y) { r <- rank(s); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0)); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
auprc <- function(s, y) {
  d <- data.table(s, y)[, .(tp = sum(y == 1), fp = sum(y == 0)), by = s][order(-s)]
  d[, `:=`(TP = cumsum(tp), FP = cumsum(fp))]; d[, `:=`(prec = TP / (TP + FP), rec = TP / sum(tp))]
  sum(diff(c(0, d$rec)) * d$prec)
}
est <- c(AUROC = auroc(a, y) - auroc(c, y), AUPRC = auprc(a, y) - auprc(c, y))
set.seed(20260925)
bs <- replicate(1000, { i <- sample.int(length(y), replace = TRUE)
  c(auroc(a[i], y[i]) - auroc(c[i], y[i]), auprc(a[i], y[i]) - auprc(c[i], y[i])) })
res <- data.table(metric = names(est), difference_AVI_minus_CADD = est,
                  lower95 = apply(bs, 1, quantile, 0.025), upper95 = apply(bs, 1, quantile, 0.975),
                  replicates_with_difference_le_0 = apply(bs, 1, function(v) sum(v <= 0)), B = 1000L, seed = 20260925L)
print(res, digits = 4)
fwrite(res, "C:/cadd/calibration/results/v17_avi/summary/table_discrimination_paired.tsv", sep = "\t")
