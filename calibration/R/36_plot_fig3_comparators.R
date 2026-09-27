# Figure 3: head-to-head comparison on the 2026 intersection reference set (Aim 6).
# (A) AUROC with 95% bootstrap intervals; (B) evidence levels assigned to P/LP and B/LB reference
# variants with the calibrated thresholds of each predictor (median of five seeds; in-sample).
# Built only from table_discrimination.tsv and table_evidence_levels.tsv (script 21 outputs).
# Revision 2026-09-27 (author's comments): notes moved to the figure legend; percentages printed in
# segments of >= 6% so that intermediate categories in panel B can be read.
# Usage: Rscript 36_plot_fig3_comparators.R <data_dir> <output_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(trailingOnly = TRUE)
tools <- c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD"); labs_t <- c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD v1.7")
disc <- fread(file.path(a[1], "table_discrimination.tsv"))[, tool := factor(tool, levels = rev(tools), labels = rev(labs_t))]
ev <- fread(file.path(a[1], "table_evidence_levels.tsv"))
ev <- melt(ev, id.vars = c("tool", "evidence"), variable.name = "set", value.name = "pct")
lv <- c("PP3_Strong", "PP3_Moderate", "PP3_Supporting", "Indeterminate", "BP4_Supporting", "BP4_Moderate", "BP4_Strong")
ev[, evidence := factor(evidence, levels = rev(lv), labels = rev(sub("_", " ", lv)))]
ev[, tool := factor(tool, levels = rev(tools), labels = rev(labs_t))]
ev[, set := factor(set, levels = c("P/LP", "B/LB"), labels = c("ClinVar P/LP", "ClinVar B/LB"))]
cols <- c("PP3 Strong" = "#67001F", "PP3 Moderate" = "#B2182B", "PP3 Supporting" = "#EF8A62", "Indeterminate" = "#D9D9D9",
          "BP4 Supporting" = "#67A9CF", "BP4 Moderate" = "#2166AC", "BP4 Strong" = "#053061")
dark <- c("PP3 Strong", "PP3 Moderate", "BP4 Moderate", "BP4 Strong")
ev[, txt_col := fifelse(evidence %in% dark, "white", "black")]
pa <- ggplot(disc, aes(AUROC, tool)) + geom_pointrange(aes(xmin = AUROC_lo, xmax = AUROC_hi), size = 0.3) +
  geom_text(aes(label = sprintf("%.3f", AUROC)), nudge_y = 0.3, size = 2.6) +
  labs(x = "AUROC (95% bootstrap interval)", y = NULL) + theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())
pb <- ggplot(ev, aes(pct, tool, fill = evidence)) + geom_col(width = 0.75) +
  # labels drawn for every segment (empty below 6%) so that label stacking matches the bars
  geom_text(aes(label = ifelse(pct >= 6, round(pct), ""), colour = txt_col, group = evidence), position = position_stack(vjust = 0.5), size = 2.3, show.legend = FALSE) +
  scale_colour_identity() + facet_wrap(~set) +
  scale_fill_manual(values = cols, name = NULL, breaks = sub("_", " ", lv)) +
  labs(x = "Reference variants (%)", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(nrow = 2))
p <- pa + pb + plot_layout(widths = c(1, 2.2)) + plot_annotation(tag_levels = "A")
ggsave(paste0(a[2], ".png"), p, width = 11, height = 4.3, dpi = 300, device = ragg::agg_png)
ggsave(paste0(a[2], ".pdf"), p, width = 11, height = 4.3, device = cairo_pdf)
cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
