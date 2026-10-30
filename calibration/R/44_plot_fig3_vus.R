# Figure 3 (manuscript v5; panels swapped in v9, 2026-10-29, so that they follow the order of citation):
# calibrated computational evidence for missense VUS.
# (A) VUS of the January 2024 release resolved by September 2026 (amendment A15): distribution of the
#     evidence that thresholds calibrated on the January 2024 release would have given, separately for
#     variants resolved to P/LP and to B/LB (gnomAD AF < 0.01).
# (B) VUS of the September 2026 ClinVar release: distribution of PP3/BP4 levels per predictor
#     (intersection thresholds; former Figure 5, script 38).
# Built only from table_vus_evidence.tsv (script 23) and table_interval_LR.tsv (script 42).
# Usage: Rscript 44_plot_fig3_vus.R <data_dir> <output_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(trailingOnly = TRUE)
tools <- c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD"); labs_t <- c("REVEL", "BayesDel", "AlphaMissense", "AVI", "CADD v1.7")
lv <- c("PP3_Strong", "PP3_Moderate", "PP3_Supporting", "Indeterminate", "BP4_Supporting", "BP4_Moderate", "BP4_Strong")
cols <- c("PP3 Strong" = "#67001F", "PP3 Moderate" = "#B2182B", "PP3 Supporting" = "#EF8A62", "Indeterminate" = "#D9D9D9",
          "BP4 Supporting" = "#67A9CF", "BP4 Moderate" = "#2166AC", "BP4 Strong" = "#053061")
dark <- c("PP3 Strong", "PP3 Moderate", "BP4 Moderate", "BP4 Strong")
prep <- function(v) {
  v[, evidence := factor(evidence, levels = rev(lv), labels = rev(sub("_", " ", lv)))]
  v[, txt_col := fifelse(evidence %in% dark, "white", "black")]; v }
bars <- function(v, xlab) {
  ggplot(v, aes(pct, tool, fill = evidence)) + geom_col(width = 0.75) +
    geom_text(aes(label = ifelse(pct >= 3, round(pct), ""), colour = txt_col, group = evidence), position = position_stack(vjust = 0.5),
              size = 2.3, show.legend = FALSE) + scale_colour_identity() +
    geom_text(data = unique(v[, intersect(c("tool", "n", "outcome"), names(v)), with = FALSE]), aes(x = 101, y = tool, label = paste0("n = ", format(n, big.mark = ","))), inherit.aes = FALSE, hjust = 0, size = 2.3) +
    scale_fill_manual(values = cols, name = NULL, breaks = sub("_", " ", lv), drop = FALSE) +
    scale_x_continuous(limits = c(0, 122), breaks = seq(0, 100, 25)) +
    labs(x = xlab, y = NULL) + theme_bw(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")) +
    guides(fill = guide_legend(nrow = 1)) }
# all 2026 VUS (panel B)
v <- fread(file.path(a[1], "table_vus_evidence.tsv"))[thresholds == "intersection" & subset == "all" & !grepl("VeryStrong", evidence) & tool %in% tools]
v[, n := sum(N), by = tool][, tool := factor(tool, levels = rev(tools), labels = rev(labs_t))]
pa <- bars(prep(v), "ClinVar 2026 missense VUS (%)")
# 2024 VUS resolved by 2026 (panel A)
r <- fread(file.path(a[1], "table_interval_LR.tsv"))[thresholds == "r2024" & subset == "AF<0.01" & !grepl("VeryStrong", evidence) & tool %in% tools]
r <- melt(r[, .(tool, evidence, `Resolved to P/LP` = n_PLP, `Resolved to B/LB` = n_BLB)], id.vars = c("tool", "evidence"), variable.name = "outcome", value.name = "N")
r[, n := sum(N), by = .(tool, outcome)][, pct := 100 * N / n][, tool := factor(tool, levels = rev(tools), labels = rev(labs_t))]
pb <- bars(prep(r), "ClinVar 2024 missense VUS resolved by 2026 (%)") + facet_wrap(~outcome, nrow = 1)
p <- (pb / pa) + plot_layout(heights = c(1, 1), guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
ggsave(paste0(a[2], ".png"), p, width = 9.5, height = 6.4, dpi = 300, device = ragg::agg_png)
ggsave(paste0(a[2], ".pdf"), p, width = 9.5, height = 6.4, device = cairo_pdf)
cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
