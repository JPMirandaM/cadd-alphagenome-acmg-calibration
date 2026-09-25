# Monte Carlo (seed) sensitivity of the bootstrap-based thresholds: reruns the CADD v1.6
# reproduction (same data and parameters as 01_reproduce_pejaver_cadd16.R, B = 10,000) with
# additional random seeds and records every threshold per seed. The original MATLAB code used
# rng('shuffle'), so the published thresholds are themselves a single random draw.
# Usage: Rscript 05b_seed_sensitivity_cadd16.R <seed1,seed2,...> [B] [workers]
suppressPackageStartupMessages(library(data.table))
args    <- commandArgs(trailingOnly = TRUE)
seeds   <- as.integer(strsplit(args[1], ",")[[1]])
B       <- if (length(args) >= 2) as.integer(args[2]) else 10000L
workers <- if (length(args) >= 3) as.integer(args[3]) else 10L
source("C:/cadd/calibration/R/pejaver_port.R")

indir  <- "C:/cadd/pj/data/data_combined_tools_split"
outdir <- sprintf("C:/cadd/calibration/results/cadd16_seed_sensitivity_B%d", B)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

D <- fread(file.path(indir, "CADDv1.6_PHRED_PLP_BLB_predictions.txt"), header = FALSE)
U <- fread(file.path(indir, "CADDv1.6_PHRED_PLP_U_predictions.txt"), header = FALSE)
x <- D$V1; y <- D$V2; g <- U[V2 == 0, V1]

for (s in seeds) {
  f <- file.path(outdir, sprintf("thresholds_seed%d.tsv", s))
  if (file.exists(f)) { cat("seed", s, "already done\n"); next }
  t0 <- Sys.time()
  res <- calibrate(x, y, g, B = B, workers = workers, increment = 0.1, seed = s)
  th <- data.table(seed = s, res$thresholds)
  # NA fraction per level, needed to interpret the one-sided bound
  na <- colMeans(is.na(res$boot_thresholds))
  th[, PP3_boot_NA_frac := na[paste0("PP3_", level)]][, BP4_boot_NA_frac := na[paste0("BP4_", level)]]
  fwrite(th, f, sep = "\t")
  fwrite(res$boot_thresholds, file.path(outdir, sprintf("bootstrap_thresholds_seed%d.tsv", s)), sep = "\t")
  cat(sprintf("seed %d done in %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(th, digits = 4)
}

all <- rbindlist(lapply(list.files(outdir, "^thresholds_seed\\d+\\.tsv$", full.names = TRUE), fread))
fwrite(all, file.path(outdir, "thresholds_all_seeds.tsv"), sep = "\t")
cat("\nDiscounted thresholds across seeds:\n")
print(dcast(melt(all, id.vars = c("seed", "level"), measure.vars = c("PP3_discounted", "BP4_discounted")),
            variable + level ~ seed, value.var = "value"))
