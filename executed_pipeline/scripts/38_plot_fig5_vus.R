# Figure 5: computational evidence for ClinVar missense VUS (Aim 7): distribution of PP3/BP4 levels
# for each predictor with the intersection-calibrated thresholds (all VUS, no frequency filter).
# Built only from table_vus_evidence.tsv (script 23).
# Usage: Rscript 38_plot_fig5_vus.R <data_dir> <output_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
v <- fread(file.path(a[1], "table_vus_evidence.tsv"))[thresholds == "intersection" & subset == "all" & !grepl("VeryStrong", evidence)]
tools <- c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD"); labs_t <- c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD v1.7")
v <- v[tool %in% tools][, tool := factor(tool, levels = rev(tools), labels = rev(labs_t))]
lv <- c("PP3_Strong", "PP3_Moderate", "PP3_Supporting", "Indeterminate", "BP4_Supporting", "BP4_Moderate", "BP4_Strong")
v[, evidence := factor(evidence, levels = rev(lv), labels = rev(sub("_", " ", lv)))]
cols <- c("PP3 Strong" = "#67001F", "PP3 Moderate" = "#B2182B", "PP3 Supporting" = "#EF8A62", "Indeterminate" = "#D9D9D9",
          "BP4 Supporting" = "#67A9CF", "BP4 Moderate" = "#2166AC", "BP4 Strong" = "#053061")
n_tool <- v[, .(n = sum(N)), by = tool]
p <- ggplot(v, aes(pct, tool, fill = evidence)) + geom_col(width = 0.75) +
  geom_text(data = n_tool, aes(x = 101, y = tool, label = paste0("n = ", format(n, big.mark = ","))), inherit.aes = FALSE, hjust = 0, size = 2.5) +
  scale_fill_manual(values = cols, name = NULL, breaks = sub("_", " ", lv)) +
  scale_x_continuous(limits = c(0, 118), breaks = seq(0, 100, 25)) +
  labs(x = "ClinVar missense VUS (%)", y = NULL) +   # title and notes moved to the figure legend (2026-09-27)
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.caption = element_text(size = 7)) +
  guides(fill = guide_legend(nrow = 1))
ggsave(paste0(a[2], ".png"), p, width = 9.5, height = 3.6, dpi = 300, device = ragg::agg_png)
ggsave(paste0(a[2], ".pdf"), p, width = 9.5, height = 3.6, device = cairo_pdf)
cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
