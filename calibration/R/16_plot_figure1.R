# Figure 1: local posterior probability curves of CADD v1.7 PHRED and AVI Phred (main analysis,
# seed 20260925), pathogenic and benign panels, in the style of Pejaver et al. (2022) Figure 3.
# Built only from figures/data/Fig1_cadd17_avi_calibration/<tool>/*.tsv (see plot_posteriors.R).
# Usage: Rscript 16_plot_figure1.R [output stem]
# The output stem is optional and defaults to the path used for the original submission, so a call
# without arguments behaves exactly as before (2026-11-05).
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("C:/cadd/calibration/R/plot_posteriors.R")
a <- commandArgs(trailingOnly = TRUE)
out <- if (length(a) >= 1L && nzchar(a[1])) a[1] else "C:/cadd/figures/Fig1_cadd17_avi_calibration"
d <- "C:/cadd/figures/data/Fig1_cadd17_avi_calibration"
p1 <- plot_posteriors(file.path(d, "CADD"), "CADD v1.7 PHRED")
p2 <- plot_posteriors(file.path(d, "AVI"), "AlphaGenome AVI Phred")
fig <- (p1 / p2) + plot_annotation(tag_levels = "A")   # notes moved to the figure legend (2026-09-27)
ggsave(paste0(out, ".png"), fig, width = 11, height = 11.5, dpi = 300, device = ragg::agg_png)
ggsave(paste0(out, ".pdf"), fig, width = 11, height = 11.5, device = cairo_pdf)
cat("Figure 1 written\n")
