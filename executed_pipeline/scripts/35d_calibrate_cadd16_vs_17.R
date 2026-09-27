# Aim 2b (A14), step d: CADD v1.6 vs v1.7 on identical reference data (2026 main reference sets).
# Both versions are calibrated on the SAME variants: reference and gnomAD SNVs with a v1.6 score
# (the v1.7 main-analysis thresholds are reused only if this subset equals the main set).
# Five seeds, B = 10,000, increment 0.1, validated fast port. Also: AUROC of each version with the
# paired bootstrap difference (1,000 replicates) and the cross-tabulation of evidence levels
# (each version at its own median thresholds) on the ClinVar reference variants.
# Output: C:/cadd/calibration/results/v17_avi/cadd16_vs_17/
# Usage: Rscript 35d_calibrate_cadd16_vs_17.R [workers]
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R")
source("C:/cadd/calibration/R/pejaver_port_fast.R")
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 8L
B <- 10000L; seeds <- c(20260925L, 101L, 202L, 303L, 404L)
root <- "C:/cadd/calibration/results/v17_avi/cadd16_vs_17"; dir.create(root, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
v16 <- readRDS("C:/cadd/data/refsets/cadd16_scores_refsets.rds")[, c(k, "CADD16_PHRED"), with = FALSE]
clin <- merge(readRDS("C:/cadd/data/refsets/clinvar_ref.rds"), v16, by = k, all.x = TRUE)
gno <- merge(readRDS("C:/cadd/data/refsets/gnomad_ref.rds"), v16, by = k, all.x = TRUE)
cat(sprintf("ClinVar reference: %d, with v1.6 score %d | gnomAD set: %d, with v1.6 score %d\n",
            nrow(clin), sum(!is.na(clin$CADD16_PHRED)), nrow(gno), sum(!is.na(gno$CADD16_PHRED))))
clin <- clin[!is.na(CADD16_PHRED) & !is.na(CADD_PHRED)]; gno <- gno[!is.na(CADD16_PHRED) & !is.na(CADD_PHRED)]
cat(sprintf("Comparison set: %d P/LP, %d B/LB, %d gnomAD\n", sum(clin$y == 1), sum(clin$y == 0), nrow(gno)))
tools <- c(CADD16 = "CADD16_PHRED", CADD17 = "CADD_PHRED")

for (tn in names(tools)) for (s in seeds) {
  out <- file.path(root, tn, paste0("seed", s))
  if (file.exists(file.path(out, "thresholds.tsv"))) { cat("skip", tn, s, "\n"); next }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  x <- clin[[tools[[tn]]]]; y <- clin$y; g <- gno[[tools[[tn]]]]
  t0 <- Sys.time()
  res <- calibrate_fast(x, y, g, B = B, workers = workers, seed = s)
  fwrite(res$thresholds, file.path(out, "thresholds.tsv"), sep = "\t")
  fwrite(res$curves, file.path(out, "posterior_curves.tsv"), sep = "\t")
  fwrite(res$boot_thresholds, file.path(out, "bootstrap_thresholds.tsv"), sep = "\t")
  writeLines(capture.output(str(res$params)), file.path(out, "parameters.txt"))
  cat(sprintf("[%s] %s seed %d: %.1f min\n", format(Sys.time(), "%H:%M"), tn, s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

all <- rbindlist(lapply(names(tools), function(tn) rbindlist(lapply(seeds, function(s)
  data.table(tool = tn, seed = s, fread(file.path(root, tn, paste0("seed", s), "thresholds.tsv"), colClasses = c(level = "character")))))))
fwrite(all, file.path(root, "thresholds_all_runs.tsv"), sep = "\t")
med <- all[, .(PP3 = median(PP3_discounted), PP3_min = min(PP3_discounted), PP3_max = max(PP3_discounted),
               BP4 = median(BP4_discounted), BP4_min = min(BP4_discounted), BP4_max = max(BP4_discounted)), by = .(tool, level)]
fwrite(med, file.path(root, "thresholds_median.tsv"), sep = "\t")
cat("\n== Thresholds (median [range] over five seeds) ==\n"); print(med)

auroc <- function(s, y) { r <- rank(s); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0)); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
set.seed(20260925)
bt <- replicate(1000, { i <- sample.int(nrow(clin), replace = TRUE)
  c(auroc(clin$CADD16_PHRED[i], clin$y[i]), auroc(clin$CADD_PHRED[i], clin$y[i])) })
disc <- data.table(metric = c("AUROC CADD v1.6", "AUROC CADD v1.7", "difference v1.7 - v1.6"),
                   estimate = c(auroc(clin$CADD16_PHRED, clin$y), auroc(clin$CADD_PHRED, clin$y),
                                auroc(clin$CADD_PHRED, clin$y) - auroc(clin$CADD16_PHRED, clin$y)),
                   lower95 = c(quantile(bt[1, ], 0.025), quantile(bt[2, ], 0.025), quantile(bt[2, ] - bt[1, ], 0.025)),
                   upper95 = c(quantile(bt[1, ], 0.975), quantile(bt[2, ], 0.975), quantile(bt[2, ] - bt[1, ], 0.975)))
fwrite(disc, file.path(root, "discrimination.tsv"), sep = "\t"); cat("\n== Discrimination ==\n"); print(disc)

lvl <- function(s, tn) { t <- med[tool == tn]; out <- rep("Indeterminate", length(s))
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { p <- t[level == L, PP3]; b <- t[level == L, BP4]
    if (length(p) && !is.na(p)) out[s >= p] <- paste0("PP3_", L); if (length(b) && !is.na(b)) out[s <= b] <- paste0("BP4_", L) }
  out }
clin[, `:=`(lv16 = lvl(CADD16_PHRED, "CADD16"), lv17 = lvl(CADD_PHRED, "CADD17"))]
xt <- clin[, .N, by = .(y, lv16, lv17)][order(y, lv16, lv17)]
fwrite(xt, file.path(root, "evidence_crosstab.tsv"), sep = "\t")
chg <- clin[, .(n = .N, changed = sum(lv16 != lv17), pct_changed = round(100 * mean(lv16 != lv17), 1),
                direction_reversed = sum((grepl("^PP3", lv16) & grepl("^BP4", lv17)) | (grepl("^BP4", lv16) & grepl("^PP3", lv17)))),
            by = .(set = fifelse(y == 1, "P/LP", "B/LB"))]
fwrite(chg, file.path(root, "evidence_changes.tsv"), sep = "\t"); cat("\n== Evidence level changes v1.6 -> v1.7 ==\n"); print(chg)
cat("Spearman correlation v1.6 vs v1.7 (reference variants):", round(cor(clin$CADD16_PHRED, clin$CADD_PHRED, method = "spearman"), 3), "\n")
