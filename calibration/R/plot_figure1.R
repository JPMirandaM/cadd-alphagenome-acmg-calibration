# Figure 1: local posterior probability curves of CADD v1.7 PHRED and AVI Phred (main analysis,
# seed 20260925), pathogenic and benign panels, in the style of Pejaver et al. (2022) Figure 3.
# Built only from figures/data/Fig1_cadd17_avi_calibration/<tool>/*.tsv (see plot_posteriors.R).
# Usage: Rscript plot_figure1.R
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("C:/cadd/calibration/R/plot_posteriors.R")
d <- "C:/cadd/figures/data/Fig1_cadd17_avi_calibration"
p1 <- plot_posteriors(file.path(d, "CADD"), "CADD v1.7 PHRED")
p2 <- plot_posteriors(file.path(d, "AVI"), "AlphaGenome AVI Phred")
fig <- (p1 / p2) + plot_annotation(tag_levels = "A",
  caption = "ClinVar (2026-09-23) rare missense variants, MANE Select; 10,000 bootstrap replicates; seed 20260925. Dashed vertical lines: final thresholds (one-sided 95% bound).")
ggsave("C:/cadd/figures/Fig1_cadd17_avi_calibration.png", fig, width = 11, height = 11.5, dpi = 300, device = ragg::agg_png)
ggsave("C:/cadd/figures/Fig1_cadd17_avi_calibration.pdf", fig, width = 11, height = 11.5, device = cairo_pdf)
cat("Figure 1 written\n")
