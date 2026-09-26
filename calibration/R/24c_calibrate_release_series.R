# Amendment A6: threshold drift across ClinVar releases and before/after the 2022 ClinGen
# recommendations, for REVEL, BayesDel, AlphaMissense, CADD v1.7 and AVI (AVI = negative control).
#  (a) releases 2019-12-02, 2022-11-19, 2024-01-27, 2026-09-23: variants scored by all five
#  (b) 2026 main set split into "pre" (in the 2022-11-19 reference set) and "post" (absent from it)
# 1 seed (20260925), B = 10,000, same precision/increments as amendment A3.
# Output: C:/cadd/calibration/results/v17_avi/release_series/{runs/, summary tables}
# Usage: Rscript 24c_calibrate_release_series.R [workers]
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/pejaver_port_fast.R")
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 10L
rs <- "C:/cadd/calibration/results/v17_avi/release_series"; runs <- file.path(rs, "runs")
k <- c("chrom", "pos", "ref", "alt")
sc <- readRDS("C:/cadd/data/refsets/comparator_scores_all_releases.rds")
sc[, `:=`(REVEL = round(REVEL, 3), AM = round(AM, 4), BayesDel = round(BayesDel, 2))]
tools <- list(REVEL = list(col = "REVEL", inc = 0.001), BayesDel = list(col = "BayesDel", inc = 0.01),
              AlphaMissense = list(col = "AM", inc = 0.001), CADD = list(col = "CADD_PHRED", inc = 0.1),
              AVI = list(col = "AVI_PHRED", inc = 0.1))
cols <- vapply(tools, `[[`, "", "col")
prep <- function(s) {
  c <- merge(s$clin[, c(k, "y", "CADD_PHRED", "AVI_PHRED"), with = FALSE], sc, by = k)
  g <- merge(s$gno[, c(k, "CADD_PHRED", "AVI_PHRED"), with = FALSE], sc, by = k)
  list(clin = c[complete.cases(c[, ..cols])], gno = g[complete.cases(g[, ..cols])])
}
sets <- list(
  r2019 = prep(readRDS(file.path(rs, "refsets_r2019.rds"))),
  r2022 = prep(readRDS(file.path(rs, "refsets_r2022.rds"))),
  r2024 = prep(readRDS("C:/cadd/calibration/results/v17_avi/temporal_holdout/refsets_2024.rds")),
  r2026 = prep(list(clin = readRDS("C:/cadd/data/refsets/clinvar_ref.rds"), gno = readRDS("C:/cadd/data/refsets/gnomad_ref.rds"))))
pre_keys <- readRDS(file.path(rs, "refsets_r2022.rds"))$clin[, ..k]
sets$post2022 <- list(clin = sets$r2026$clin[!pre_keys, on = k], gno = sets$r2026$gno)
sets$pre2022  <- list(clin = sets$r2026$clin[pre_keys, on = k, nomatch = NULL], gno = sets$r2026$gno)
sizes <- rbindlist(lapply(names(sets), function(n) data.table(set = n, P_LP = sets[[n]]$clin[y == 1, .N],
                   B_LB = sets[[n]]$clin[y == 0, .N], gnomAD = nrow(sets[[n]]$gno))))
print(sizes); fwrite(sizes, file.path(rs, "set_sizes.tsv"), sep = "\t")

for (n in names(sets)) for (tn in names(tools)) {
  out <- file.path(runs, n, tn)
  if (file.exists(file.path(out, "thresholds.tsv"))) { cat("skip", n, tn, "\n"); next }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  r <- calibrate_fast(sets[[n]]$clin[[tools[[tn]]$col]], sets[[n]]$clin$y, sets[[n]]$gno[[tools[[tn]]$col]],
                      B = 10000, workers = workers, seed = 20260925, increment = tools[[tn]]$inc)
  fwrite(r$thresholds, file.path(out, "thresholds.tsv"), sep = "\t")
  fwrite(r$curves, file.path(out, "posterior_curves.tsv"), sep = "\t")
  fwrite(r$boot_thresholds, file.path(out, "bootstrap_thresholds.tsv"), sep = "\t")
  writeLines(capture.output(str(r$params)), file.path(out, "parameters.txt"))
  cat(sprintf("[%s] %s %s: %.1f min\n", format(Sys.time(), "%H:%M"), n, tn, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

num <- c("PP3_point", "PP3_discounted", "BP4_point", "BP4_discounted", "PP3_boot_NA_frac", "BP4_boot_NA_frac")
all <- rbindlist(lapply(names(sets), function(n) rbindlist(lapply(names(tools), function(tn)
  data.table(set = n, tool = tn, fread(file.path(runs, n, tn, "thresholds.tsv"),
             colClasses = c(level = "character", setNames(rep("numeric", length(num)), num))))))))
fwrite(all, file.path(rs, "thresholds_all.tsv"), sep = "\t")
lv <- c("Supporting", "Moderate", "Strong", "VeryStrong")
w <- function(crit) { x <- dcast(all, tool + level ~ set, value.var = crit)[, level := factor(level, levels = lv)]
  setcolorder(x, c("tool", "level", "r2019", "r2022", "r2024", "r2026", "pre2022", "post2022")); x[order(tool, level)] }
cat("\n== PP3 thresholds (score >=) ==\n"); print(w("PP3_discounted"))
cat("\n== BP4 thresholds (score <=) ==\n"); print(w("BP4_discounted"))
fwrite(w("PP3_discounted"), file.path(rs, "table_PP3_by_set.tsv"), sep = "\t")
fwrite(w("BP4_discounted"), file.path(rs, "table_BP4_by_set.tsv"), sep = "\t")
