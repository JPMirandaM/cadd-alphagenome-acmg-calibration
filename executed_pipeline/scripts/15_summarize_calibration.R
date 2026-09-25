# Summarises the CADD v1.7 / AVI calibration (script 14): main-analysis thresholds (median and range
# over 5 seeds), sensitivity analyses S1-S4, discrimination, and the fraction of reference variants
# reaching each evidence level. Writes tables to C:/cadd/calibration/results/v17_avi/summary/.
suppressPackageStartupMessages(library(data.table))
root <- "C:/cadd/calibration/results/v17_avi"; out <- file.path(root, "summary")
dir.create(out, showWarnings = FALSE)
files <- list.files(root, "^thresholds\\.tsv$", recursive = TRUE, full.names = TRUE)   # not bootstrap_thresholds.tsv
num <- c("PP3_point", "PP3_discounted", "BP4_point", "BP4_discounted", "PP3_boot_NA_frac", "BP4_boot_NA_frac")
all <- rbindlist(lapply(files, function(f) {
  p <- strsplit(sub(paste0(root, "/"), "", f), "/")[[1]]
  data.table(analysis = p[1], tool = p[2], seed = sub("seed", "", p[3]),
             fread(f, colClasses = c(level = "character", setNames(rep("numeric", length(num)), num))))
}))
all <- all[analysis != "summary"]
cat("runs:", uniqueN(all[, .(analysis, tool, seed)]), "\n")
lev_order <- c("Supporting", "Moderate", "Strong", "VeryStrong")
all[, level := factor(level, levels = lev_order)]
fwrite(all, file.path(out, "thresholds_all_runs.tsv"), sep = "\t")

fmt <- function(v) if (all(is.na(v))) "not reached" else if (any(is.na(v))) sprintf("%s (reached in %d/%d seeds)", format(median(v, na.rm = TRUE)), sum(!is.na(v)), length(v)) else
  if (min(v) == max(v)) format(min(v)) else sprintf("%s [%s-%s]", format(median(v)), format(min(v)), format(max(v)))
main <- all[analysis == "main", .(PP3 = fmt(PP3_discounted), BP4 = fmt(BP4_discounted),
                                   PP3_fail = sprintf("%.3f", median(PP3_boot_NA_frac)), BP4_fail = sprintf("%.3f", median(BP4_boot_NA_frac))),
            by = .(tool, level)][order(tool, level)]
cat("\n== Main analysis (5 seeds): PP3 score >= threshold; BP4 score <= threshold ==\n"); print(main)
fwrite(main, file.path(out, "table_main_thresholds.tsv"), sep = "\t")

sens <- dcast(all[analysis != "main" | seed == "20260925"], tool + level ~ analysis, value.var = "PP3_discounted")
sensb <- dcast(all[analysis != "main" | seed == "20260925"], tool + level ~ analysis, value.var = "BP4_discounted")
cat("\n== PP3 thresholds, main (seed 20260925) and sensitivity analyses ==\n"); print(sens[order(tool, level)])
cat("\n== BP4 thresholds, main (seed 20260925) and sensitivity analyses ==\n"); print(sensb[order(tool, level)])
fwrite(sens, file.path(out, "table_sensitivity_PP3.tsv"), sep = "\t"); fwrite(sensb, file.path(out, "table_sensitivity_BP4.tsv"), sep = "\t")

disc <- rbindlist(lapply(c("CADD", "AVI"), function(t) data.table(tool = t, fread(file.path(root, "main", t, "seed20260925", "discrimination.tsv")))))
cat("\n== Discrimination (main reference set) ==\n"); print(disc, digits = 4)
fwrite(disc, file.path(out, "table_discrimination.tsv"), sep = "\t")

# Fraction of reference variants meeting each level with the main-analysis median thresholds
clin <- readRDS("C:/cadd/data/refsets/clinvar_ref.rds")
col <- c(CADD = "CADD_PHRED", AVI = "AVI_PHRED")
med <- all[analysis == "main", .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)]
yield <- rbindlist(lapply(names(col), function(t) {
  s <- clin[[col[[t]]]]; m <- med[tool == t]
  thr <- function(crit, lv) m[level == lv, get(crit)]
  lvl_p <- fcase(!is.na(thr("PP3", "VeryStrong")) & s >= thr("PP3", "VeryStrong"), "PP3_VeryStrong",
                 !is.na(thr("PP3", "Strong")) & s >= thr("PP3", "Strong"), "PP3_Strong",
                 s >= thr("PP3", "Moderate"), "PP3_Moderate", s >= thr("PP3", "Supporting"), "PP3_Supporting",
                 !is.na(thr("BP4", "VeryStrong")) & s <= thr("BP4", "VeryStrong"), "BP4_VeryStrong",
                 !is.na(thr("BP4", "Strong")) & s <= thr("BP4", "Strong"), "BP4_Strong",
                 s <= thr("BP4", "Moderate"), "BP4_Moderate", s <= thr("BP4", "Supporting"), "BP4_Supporting",
                 default = "Indeterminate")
  data.table(tool = t, class = ifelse(clin$y == 1, "P/LP", "B/LB"), evidence = lvl_p)
}))
yt <- yield[, .N, by = .(tool, class, evidence)][, pct := round(100 * N / sum(N), 1), by = .(tool, class)]
yt <- dcast(yt, tool + evidence ~ class, value.var = "pct", fill = 0)
ord <- c("PP3_VeryStrong", "PP3_Strong", "PP3_Moderate", "PP3_Supporting", "Indeterminate", "BP4_Supporting", "BP4_Moderate", "BP4_Strong", "BP4_VeryStrong")
yt[, evidence := factor(evidence, levels = ord)]
cat("\n== % of ClinVar reference variants assigned to each evidence level (main median thresholds) ==\n"); print(yt[order(tool, evidence)])
fwrite(yt[order(tool, evidence)], file.path(out, "table_evidence_levels_reference.tsv"), sep = "\t")
