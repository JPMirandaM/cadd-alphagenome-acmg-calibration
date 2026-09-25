# Verifies that pejaver_port_fast.R is numerically identical to the validated pejaver_port.R,
# on the CADD v1.6 data of Pejaver et al. (2022):
#   1. point-estimate posteriors: local_posteriors_counts() vs local_posteriors()
#   2. 25 bootstrap replicates with the same seeds: bootstrap_once_fast() vs bootstrap_once()
#   3. full calibration, B = 10,000, seed 20260925, 10 workers: bootstrap thresholds and final
#      thresholds vs the stored results of the original run (step 5)
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/pejaver_port_fast.R")
d <- "C:/cadd/pj/data/data_combined_tools_split"
D <- fread(file.path(d, "CADDv1.6_PHRED_PLP_BLB_predictions.txt"), header = FALSE); x <- D$V1; y <- D$V2
U <- fread(file.path(d, "CADDv1.6_PHRED_PLP_U_predictions.txt"), header = FALSE); g <- U[V2 == 0, V1]
w <- negative_weight(y, 0.0441); thrs <- score_grid(x, g); pc <- prepare_codes(x, g); k <- length(pc$v)

a <- local_posteriors(x, y, g, thrs, w)
b <- local_posteriors_counts(tabulate(pc$xc[y == 1], k), tabulate(pc$xc[y == 0], k), tabulate(pc$gc, k),
                             pc$v, thrs, w, lg = length(g))
cat("1. point estimate identical:", identical(a$post, b$post), "| windows identical:", identical(a$finalwindow, b$finalwindow), "\n")

same <- vapply(1:25, function(s) {
  set.seed(s); o <- bootstrap_once(x, y, g, thrs, w, 100, 0.03, 0.1)
  set.seed(s); f <- bootstrap_once_fast(pc$xc, y, pc$gc, pc$v, thrs, w, 100, 0.03, 0.1)
  identical(o, f)
}, logical(1))
cat("2. bootstrap replicates identical:", sum(same), "of", length(same), "\n")

t0 <- Sys.time()
r <- calibrate_fast(x, y, g, B = 10000, workers = 10, seed = 20260925)
cat(sprintf("3. fast calibration B = 10,000: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
old_bt <- fread("C:/cadd/calibration/results/cadd16_reproduction_B10000/bootstrap_thresholds.tsv")
old_th <- fread("C:/cadd/calibration/results/cadd16_reproduction_B10000/thresholds.tsv")
# Element-wise comparison (all-NA columns are read back by fread as logical, so compare as doubles)
om <- as.matrix(old_bt); storage.mode(om) <- "double"; nm <- as.matrix(r$boot_thresholds)
cat("   bootstrap thresholds: same NA pattern:", identical(is.na(om), is.na(nm)),
    "| max abs difference:", max(abs(om - nm), na.rm = TRUE), "\n")
cat("   final thresholds identical to original run:",
    isTRUE(all.equal(as.data.table(r$thresholds)[, .(PP3_point, PP3_discounted, BP4_point, BP4_discounted)],
                     old_th[, .(PP3_point, PP3_discounted, BP4_point, BP4_discounted)], check.attributes = FALSE)), "\n")
print(r$thresholds, digits = 4)
