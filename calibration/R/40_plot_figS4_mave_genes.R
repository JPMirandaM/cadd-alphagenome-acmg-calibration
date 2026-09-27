# Figure S4: AUROC of five predictors for functionally abnormal vs normal MANE missense SNVs per gene
# (MAVE validation, amendment A13; genes with >= 30 abnormal and >= 30 normal SNVs).
# Built only from auroc_per_gene.tsv (script 33c).
# Usage: Rscript 40_plot_figS4_mave_genes.R <data_dir> <output_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
g <- fread(file.path(a[1], "auroc_per_gene.tsv"))
l <- melt(g, id.vars = c("gene", "n", "n_abnormal"), measure.vars = c("AUROC_CADD", "AUROC_AVI", "AUROC_REVEL", "AUROC_BayesDel", "AUROC_AM"),
          variable.name = "tool", value.name = "AUROC")
l[, tool := factor(tool, levels = c("AUROC_REVEL", "AUROC_BayesDel", "AUROC_AM", "AUROC_AVI", "AUROC_CADD"),
                   labels = c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD v1.7"))]
ordg <- g[order(AUROC_AVI), gene]; l[, gene := factor(gene, levels = ordg)]
lab <- g[, .(gene = factor(gene, levels = ordg), txt = sprintf("%s (%d/%d)", gene, n_abnormal, n - n_abnormal))]
pal <- c("REVEL" = "#E69F00", "BayesDel" = "#56B4E9", "AlphaMissense" = "#CC79A7", "AVI" = "#1B9E77", "CADD v1.7" = "#4D4D4D")
p <- ggplot(l, aes(AUROC, gene, colour = tool, shape = tool)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey60") + geom_point(size = 1.8, alpha = 0.9) +
  scale_colour_manual(values = pal, name = NULL) + scale_shape_manual(values = c(16, 17, 15, 18, 4), name = NULL) +
  scale_y_discrete(labels = setNames(lab$txt, as.character(lab$gene))) +
  labs(x = "AUROC (functionally abnormal vs normal)", y = "Gene (abnormal/normal SNVs)",
       title = "Discrimination of MAVE functional classes by gene",
       caption = "MaveDB score sets with normal/abnormal functional classes; genes with >= 30 abnormal and >= 30 normal MANE missense SNVs; genes ordered by AVI AUROC.") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.caption = element_text(size = 7))
ggsave(paste0(a[2], ".png"), p, width = 7.5, height = 8.5, dpi = 300, device = ragg::agg_png)
ggsave(paste0(a[2], ".pdf"), p, width = 7.5, height = 8.5, device = cairo_pdf)
cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
