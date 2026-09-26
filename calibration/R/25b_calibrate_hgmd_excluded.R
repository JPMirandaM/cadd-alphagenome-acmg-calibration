# Amendment A7 (b, c): exclude SNVs at positions of single-nucleotide HGMD Public 2020.4 records from
# the ClinVar and gnomAD reference sets (intersection of the five predictors), for the 2019-12-02
# and 2026-09-23 releases, and recalibrate REVEL, BayesDel, AlphaMissense, CADD v1.7 and AVI
# (1 seed, 20260925, B = 10,000; precision/increments as in amendment A3).
# Output: C:/cadd/calibration/results/v17_avi/hgmd_excluded/
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/pejaver_port_fast.R")
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 10L
out <- "C:/cadd/calibration/results/v17_avi/hgmd_excluded"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
rs <- "C:/cadd/calibration/results/v17_avi/release_series"; k <- c("chrom", "pos", "ref", "alt")
hg <- readRDS("C:/cadd/data/refsets/hgmd_public_2020.4.rds")[start == end, .(chrom, pos = start)]
hg <- unique(hg)
cat("HGMD Public single-nucleotide positions:", nrow(hg), "\n")
sc <- readRDS("C:/cadd/data/refsets/comparator_scores_all_releases.rds")
sc[, `:=`(REVEL = round(REVEL, 3), AM = round(AM, 4), BayesDel = round(BayesDel, 2))]
tools <- list(REVEL = list(col = "REVEL", inc = 0.001), BayesDel = list(col = "BayesDel", inc = 0.01),
              AlphaMissense = list(col = "AM", inc = 0.001), CADD = list(col = "CADD_PHRED", inc = 0.1),
              AVI = list(col = "AVI_PHRED", inc = 0.1))
cols <- vapply(tools, `[[`, "", "col")
prep <- function(s) {
  c <- merge(s$clin[, c(k, "y", "CADD_PHRED", "AVI_PHRED"), with = FALSE], sc, by = k)
  g <- merge(s$gno[, c(k, "CADD_PHRED", "AVI_PHRED"), with = FALSE], sc, by = k)
  c <- c[complete.cases(c[, ..cols])]; g <- g[complete.cases(g[, ..cols])]
  cx <- c[!hg, on = .(chrom, pos)]; gx <- g[!hg, on = .(chrom, pos)]
  cat(sprintf("  ClinVar: %d P/LP, %d B/LB -> after exclusion %d P/LP, %d B/LB | gnomAD %d -> %d\n",
              c[y == 1, .N], c[y == 0, .N], cx[y == 1, .N], cx[y == 0, .N], nrow(g), nrow(gx)))
  list(clin = cx, gno = gx, removed = data.table(P_LP_removed = c[y == 1, .N] - cx[y == 1, .N],
       B_LB_removed = c[y == 0, .N] - cx[y == 0, .N], gnomAD_removed = nrow(g) - nrow(gx)))
}
cat("2019-12-02:\n"); s19 <- prep(readRDS(file.path(rs, "refsets_r2019.rds")))
cat("2026-09-23:\n"); s26 <- prep(list(clin = readRDS("C:/cadd/data/refsets/clinvar_ref.rds"), gno = readRDS("C:/cadd/data/refsets/gnomad_ref.rds")))
sets <- list(r2019_noHGMD = s19, r2026_noHGMD = s26)
fwrite(rbindlist(lapply(names(sets), function(n) data.table(set = n, P_LP = sets[[n]]$clin[y == 1, .N], B_LB = sets[[n]]$clin[y == 0, .N],
       gnomAD = nrow(sets[[n]]$gno), sets[[n]]$removed))), file.path(out, "set_sizes.tsv"), sep = "\t")

for (n in names(sets)) for (tn in names(tools)) {
  o <- file.path(out, "runs", n, tn)
  if (file.exists(file.path(o, "thresholds.tsv"))) { cat("skip", n, tn, "\n"); next }
  dir.create(o, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  r <- calibrate_fast(sets[[n]]$clin[[tools[[tn]]$col]], sets[[n]]$clin$y, sets[[n]]$gno[[tools[[tn]]$col]],
                      B = 10000, workers = workers, seed = 20260925, increment = tools[[tn]]$inc)
  fwrite(r$thresholds, file.path(o, "thresholds.tsv"), sep = "\t")
  fwrite(r$curves, file.path(o, "posterior_curves.tsv"), sep = "\t")
  fwrite(r$boot_thresholds, file.path(o, "bootstrap_thresholds.tsv"), sep = "\t")
  writeLines(capture.output(str(r$params)), file.path(o, "parameters.txt"))
  cat(sprintf("[%s] %s %s: %.1f min\n", format(Sys.time(), "%H:%M"), n, tn, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

num <- c("PP3_point", "PP3_discounted", "BP4_point", "BP4_discounted", "PP3_boot_NA_frac", "BP4_boot_NA_frac")
ex <- rbindlist(lapply(names(sets), function(n) rbindlist(lapply(names(tools), function(tn)
  data.table(set = n, tool = tn, fread(file.path(out, "runs", n, tn, "thresholds.tsv"),
             colClasses = c(level = "character", setNames(rep("numeric", length(num)), num))))))))
base <- fread(file.path(rs, "thresholds_all.tsv"))[set %in% c("r2019", "r2026")]
base[, set := paste0(set, "_all")]
all <- rbind(base[, names(ex), with = FALSE], ex)
pub <- data.table(tool = c(rep("REVEL", 7), rep("BayesDel", 5), rep("AlphaMissense", 7)),
  crit = c("PP3", "PP3", "PP3", "BP4", "BP4", "BP4", "BP4", "PP3", "PP3", "PP3", "BP4", "BP4", "PP3", "PP3", "PP3", "BP4", "BP4", "BP4", "PP3"),
  level = c("Supporting", "Moderate", "Strong", "Supporting", "Moderate", "Strong", "VeryStrong",
            "Supporting", "Moderate", "Strong", "Supporting", "Moderate",
            "Supporting", "Moderate", "Strong", "Supporting", "Moderate", "Strong", "VeryStrong"),
  published = c(0.644, 0.773, 0.932, 0.290, 0.183, 0.016, 0.003, 0.13, 0.27, 0.50, -0.18, -0.36,
                0.792, 0.906, 0.972, 0.169, 0.099, 0.070, 0.990))
long <- rbind(all[, .(set, tool, level, crit = "PP3", value = PP3_discounted)], all[, .(set, tool, level, crit = "BP4", value = BP4_discounted)])
w <- dcast(long, tool + crit + level ~ set, value.var = "value")
w <- merge(w, pub, by = c("tool", "crit", "level"), all.x = TRUE)
w[, level := factor(level, levels = c("Supporting", "Moderate", "Strong", "VeryStrong"))]
setcolorder(w, c("tool", "crit", "level", "published", "r2019_all", "r2019_noHGMD", "r2026_all", "r2026_noHGMD"))
cat("\n== Thresholds with and without HGMD Public positions (published values: REVEL/BayesDel Pejaver 2022; AlphaMissense Bergquist 2025) ==\n")
print(w[order(tool, crit, level)])
fwrite(w, file.path(out, "table_thresholds_hgmd.tsv"), sep = "\t")
