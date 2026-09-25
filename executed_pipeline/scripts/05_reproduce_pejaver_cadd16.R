# Validation of the R port: reproduces the CADD v1.6 PHRED calibration of Pejaver et al. (2022)
# using their original data (Zenodo 10.5281/zenodo.8347415, data.zip, md5 a25359e5bbb3de9ea47dd21f03950ffb).
# Usage: Rscript 01_reproduce_pejaver_cadd16.R [B] [workers]
suppressPackageStartupMessages(library(data.table))
args    <- commandArgs(trailingOnly = TRUE)
B       <- if (length(args) >= 1) as.integer(args[1]) else 10000L
workers <- if (length(args) >= 2) as.integer(args[2]) else 10L
source("C:/cadd/calibration/R/pejaver_port.R")

indir  <- "C:/cadd/pj/data/data_combined_tools_split"
outdir <- sprintf("C:/cadd/calibration/results/cadd16_reproduction_B%d", B)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# main.m, lines 134-144: x, y from PLP_BLB; g = label 0 of the PLP_U file
D <- fread(file.path(indir, "CADDv1.6_PHRED_PLP_BLB_predictions.txt"), header = FALSE)
x <- D$V1; y <- D$V2
U <- fread(file.path(indir, "CADDv1.6_PHRED_PLP_U_predictions.txt"), header = FALSE)
g <- U[V2 == 0, V1]
cat(sprintf("P/LP = %d, B/LB = %d, gnomAD = %d, B = %d\n", sum(y == 1), sum(y == 0), length(g), B))

t0 <- Sys.time()
res <- calibrate(x, y, g, B = B, workers = workers, increment = 0.1)   # CADDphred: increment 0.1
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(res$levels, digits = 4)
print(res$thresholds, digits = 4)

fwrite(res$thresholds, file.path(outdir, "thresholds.tsv"), sep = "\t")
fwrite(res$levels, file.path(outdir, "posterior_levels.tsv"), sep = "\t")
fwrite(res$curves, file.path(outdir, "posterior_curves.tsv"), sep = "\t")
if (!is.null(res$boot_thresholds)) fwrite(res$boot_thresholds, file.path(outdir, "bootstrap_thresholds.tsv"), sep = "\t")
writeLines(capture.output(str(res$params)), file.path(outdir, "parameters.txt"))
saveRDS(res, file.path(outdir, "full_result.rds"))
cat("Results in", outdir, "\n")
