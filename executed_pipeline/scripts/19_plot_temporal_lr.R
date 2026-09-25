# Figure S3: interval-based likelihood ratios on the temporal test set (variants first classified
# after ClinVar 2024-01-27), for thresholds calibrated on the 2024 release. Points: LR with 95%
# bootstrap interval; horizontal marks: LR required for the evidence strength (PP3: LR must be at
# least the mark; BP4: LR must be at most the mark). Built only from interval_LR_test_set.tsv.
# Usage: Rscript plot_temporal_lr.R <data_dir> <output_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
d <- fread(file.path(a[1], "interval_LR_test_set.tsv"))[interval != "Indeterminate"]
ord <- c("BP4_Strong", "BP4_Moderate", "BP4_Supporting", "PP3_Supporting", "PP3_Moderate")
d[, interval := factor(interval, levels = ord)]
d[, side := ifelse(grepl("^PP3", interval), "PP3 (pathogenic)", "BP4 (benign)")]
d[, tool := factor(tool, levels = c("CADD", "AVI"), labels = c("CADD v1.7", "AVI"))]
p <- ggplot(d, aes(interval, LR, colour = tool)) +
  geom_errorbar(aes(ymin = LR_required, ymax = LR_required), colour = "grey30", width = 0.6, linewidth = 0.8,
                position = position_identity()) +
  geom_pointrange(aes(ymin = lower95, ymax = upper95), position = position_dodge(width = 0.5), size = 0.35) +
  scale_y_log10() + scale_colour_manual(values = c("CADD v1.7" = "#1F77B4", "AVI" = "#D62728"), name = NULL) +
  facet_wrap(~side, scales = "free") +
  labs(x = NULL, y = "Positive likelihood ratio in interval (log scale)",
       title = "Temporal validation: thresholds from ClinVar 2024-01-27 evaluated on later-classified variants",
       caption = paste0("Test set: 22,626 P/LP and 36,294 B/LB rare missense variants absent from the 2024 set.\n",
                        "Grey marks: required LR (PP3: LR >= mark; BP4: LR <= mark). Bars: 95% bootstrap intervals.")) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                   axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(paste0(a[2], ".png"), p, width = 9, height = 5.5, dpi = 300, device = ragg::agg_png)
ggsave(paste0(a[2], ".pdf"), p, width = 9, height = 5.5, device = cairo_pdf)
cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
