# Figure 4: population frequency of MANE missense SNVs with PP3-level scores (Aim 8).
# (A) gnomAD v4.1 group with the highest AF among PP3-level variants common (> 1%) in >= 1 group;
# (B) Mexico City Prospective Study: PP3-level SNVs with AF > 1% overall and in the Indigenous
#     Mexican ancestry component, and how many of them are <= 1% in every gnomAD group.
# Built only from table_grpmax.tsv (script 22) and pp3_levels_mcps_summary.tsv (script 32).
# Revision 2026-09-27 (author's comments): notes moved to the figure legend; long category labels
# wrapped on two lines; wider spacing between panels.
# Usage: Rscript 37_plot_fig4_population.R <data_dir> <output_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(trailingOnly = TRUE)
pal <- c("CADD v1.7" = "#404040", "AVI" = "#C8C8C8")   # grayscale (author, 2026-09-27: colour figures incur a journal charge)
g <- fread(file.path(a[1], "table_grpmax.tsv"))[PP3 == "Supporting+"]
g[, tool := factor(tool, levels = c("CADD", "AVI"), labels = c("CADD v1.7", "AVI"))]
ord <- g[, .(m = mean(pct)), by = grpmax][order(m), grpmax]; g[, grpmax := factor(grpmax, levels = ord)]
pa <- ggplot(g, aes(pct, grpmax, fill = tool)) + geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = "PP3-level variants common in >= 1 group (%)", y = "gnomAD group with highest AF") + theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
m <- fread(file.path(a[1], "pp3_levels_mcps_summary.tsv"))
m <- m[, .(v1 = sum(n_mcps_af_gt_0.01), v2 = sum(n_mcps_gt_0.01_gnomad_grpmax_le_0.01),
           v3 = sum(n_imx_af_gt_0.01), v4 = sum(n_imx_gt_0.01_gnomad_grpmax_le_0.01)), by = predictor]
m <- melt(m, id.vars = "predictor", variable.name = "measure", value.name = "n")
m[, measure := factor(measure, levels = rev(c("v1", "v2", "v3", "v4")),
                      labels = rev(c("MCPS AF > 1%", "MCPS AF > 1%,\ngnomAD group max <= 1%",
                                     "Indigenous Mexican AF > 1%", "Indigenous Mexican AF > 1%,\ngnomAD group max <= 1%")))]
m[, predictor := factor(predictor, levels = c("CADD", "AVI"), labels = c("CADD v1.7", "AVI"))]
pb <- ggplot(m, aes(n, measure, fill = predictor)) + geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), hjust = -0.15, size = 2.5) +
  scale_fill_manual(values = pal, name = NULL) + scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "PP3-level MANE missense SNVs", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
p <- (pa | plot_spacer() | pb) + plot_layout(widths = c(1, 0.08, 1.1)) + plot_annotation(tag_levels = "A")
ggsave(paste0(a[2], ".png"), p, width = 11.5, height = 4.5, dpi = 300, device = ragg::agg_png)
ggsave(paste0(a[2], ".pdf"), p, width = 11.5, height = 4.5, device = cairo_pdf)
cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
