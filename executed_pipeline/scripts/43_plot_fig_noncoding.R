# Figure (noncoding): how the missense-calibrated thresholds of CADD v1.7 and AVI behave outside the
# MANE missense context (Aim 9, amendment A10; results R13).
# (A) Likelihood ratio, P(interval | P/LP) / P(interval | B/LB), of the missense BP4 intervals in each
#     ClinVar region class; LR > 1 means the "benign" interval is enriched for pathogenic variants.
# (B) AUROC (P/LP vs B/LB) of each predictor in each class.
# Built only from evaluation_missense_thresholds_by_class.tsv (script 30). Grayscale (journal colour charge).
# Usage: Rscript 43_plot_fig_noncoding.R <data_dir> <output_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(trailingOnly = TRUE)
d <- fread(file.path(a[1], "evaluation_missense_thresholds_by_class.tsv"))
lab <- c(canonical_splice = "Canonical splice site", splice_region_intronic = "Intronic splice region (3-8 nt)",
         near_intronic = "Near intronic (9-50 nt)", deep_intronic = "Deep intronic (> 50 nt)",
         `5UTR` = "5' UTR", `3UTR` = "3' UTR", upstream_5kb = "Upstream (< 5 kb)", downstream_5kb = "Downstream (< 5 kb)",
         intergenic = "Intergenic", synonymous_exon_edge = "Synonymous, exon edge", synonymous_other = "Synonymous, other")
d[, cls := factor(lab[class], levels = rev(lab))]
d[, tool := factor(tool, levels = c("CADD", "AVI"), labels = c("CADD v1.7", "AVI"))]
pal <- c("CADD v1.7" = "#404040", "AVI" = "#A0A0A0")
b <- d[grepl("^BP4", interval) & n_PLP + n_BLB >= 10 & is.finite(LR) & LR > 0]
b[, interval := factor(sub("BP4_", "BP4 ", interval), levels = c("BP4 Supporting", "BP4 Moderate", "BP4 Strong"))]
b[, `:=`(lo = pmax(lo, 1e-3), hi = pmin(hi, 1e3))]
pa <- ggplot(b, aes(LR, cls, colour = tool, shape = interval)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 1 / 2.406, linetype = "dotted", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.3, position = position_dodge(width = 0.7)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.7)) +
  scale_x_log10(breaks = c(0.01, 0.1, 0.415, 1, 10, 100), labels = c("0.01", "0.1", "0.42", "1", "10", "100")) +
  scale_colour_manual(values = pal, name = NULL) + scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(x = "Likelihood ratio of the missense BP4 interval (log scale)", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank())
u <- unique(d[, .(cls, tool, AUROC)])
pb <- ggplot(u, aes(AUROC, cls, fill = tool)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = pal, name = NULL) + coord_cartesian(xlim = c(0.5, 1)) +
  labs(x = "AUROC (P/LP vs B/LB)", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.y = element_blank())
p <- (pa | pb) + plot_layout(widths = c(1.4, 1)) + plot_annotation(tag_levels = "A")
ggsave(paste0(a[2], ".png"), p, width = 10, height = 4.8, dpi = 300, device = ragg::agg_png)
ggsave(paste0(a[2], ".pdf"), p, width = 10, height = 4.8, device = cairo_pdf)
cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
