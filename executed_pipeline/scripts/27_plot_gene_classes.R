# Figure 2: calibrated thresholds by gene class (constraint, inheritance, mechanism) for the five
# predictors, with 2.5th-97.5th percentile intervals of the 10,000 bootstrap threshold values.
# Built only from thresholds_by_class.tsv (26_gene_class_stratification.R).
# Usage: Rscript plot_gene_classes.R <data_dir> <output_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
d <- fread(file.path(a[1], "thresholds_by_class.tsv"))[level %in% c("Supporting", "Moderate") & !is.na(threshold)]
d[, class := factor(class, levels = c("LOEUF_0-2", "LOEUF_3-6", "LOEUF_7-9", "AD_only", "AR_only", "LoF", "GoF_DN"),
                    labels = c("LOEUF 0-2\n(constrained)", "LOEUF 3-6", "LOEUF 7-9\n(tolerant)", "AD only", "AR only", "LoF", "GoF / DN"))]
d[, axis := factor(axis, levels = c("constraint", "inheritance", "mechanism"), labels = c("Constraint", "Inheritance", "Mechanism"))]
d[, tool := factor(tool, levels = c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD"), labels = c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD v1.7"))]
d[, criterion := factor(paste(crit, level), levels = c("PP3 Moderate", "PP3 Supporting", "BP4 Supporting", "BP4 Moderate"))]
p <- ggplot(d, aes(class, threshold, colour = criterion)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.55), size = 0.25) +
  facet_grid(tool ~ axis, scales = "free", space = "free_x") +
  scale_colour_manual(values = c("PP3 Moderate" = "#B2182B", "PP3 Supporting" = "#EF8A62",
                                 "BP4 Supporting" = "#67A9CF", "BP4 Moderate" = "#2166AC"), name = NULL) +
  labs(x = NULL, y = "Calibrated threshold (score)") +   # title and notes moved to the figure legend (2026-09-27)
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                  strip.text.y = element_text(angle = 0), axis.text.x = element_text(size = 7))
ggsave(paste0(a[2], ".png"), p, width = 10, height = 11, dpi = 300, device = ragg::agg_png)
ggsave(paste0(a[2], ".pdf"), p, width = 10, height = 11, device = cairo_pdf)
cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
