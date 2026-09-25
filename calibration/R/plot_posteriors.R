# Local posterior probability figure in the style of Pejaver et al. (2022), Figure 3
# (reproduces plotting/plot_both_posteriors_pub.m).
# Built ONLY from the .tsv files, so it can be modified without recomputing:
#   posterior_curves.tsv   score, post_p, post_p_p5, post_b, post_b_p5
#   posterior_levels.tsv   level, LR, Post_p, Post_b
#   thresholds.tsv         level, PP3_point, PP3_discounted, BP4_point, BP4_discounted
# Usage: Rscript plot_posteriors.R <data_dir> <output_path_without_extension> "<title>"
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

plot_posteriors <- function(datadir, title, show_thresholds = TRUE) {
  cur <- fread(file.path(datadir, "posterior_curves.tsv"))
  lev <- fread(file.path(datadir, "posterior_levels.tsv"))
  thr <- fread(file.path(datadir, "thresholds.tsv"))
  lev[, level := factor(level, levels = c("VeryStrong", "Strong", "Moderate", "Supporting"))]
  lty <- c(VeryStrong = "solid", Strong = "dotdash", Moderate = "dashed", Supporting = "dotted")
  base <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.box = "vertical",
          legend.text = element_text(size = 8), legend.title = element_text(size = 9),
          legend.key.width = unit(1.6, "lines"), plot.title = element_text(size = 11))

  curve_df <- function(p, p5) rbind(
    data.table(score = cur$score, value = cur[[p]], series = "Point estimate"),
    data.table(score = cur$score, value = cur[[p5]], series = "Bootstrap 5th percentile"))
  cols_curve <- c("Point estimate" = "black", "Bootstrap 5th percentile" = "grey50")

  vlines <- function(col) if (show_thresholds) {
    v <- thr[!is.na(get(col))]
    if (nrow(v)) geom_vline(data = v, aes(xintercept = .data[[col]]), colour = "grey35",
                            linewidth = 0.3, linetype = "longdash")
  }

  pp <- ggplot() +
    vlines("PP3_discounted") +
    geom_hline(data = lev, aes(yintercept = Post_p, linetype = level), colour = "#D62728", linewidth = 0.8) +
    geom_line(data = curve_df("post_p", "post_p_p5")[!is.na(value)], aes(score, value, colour = series), linewidth = 0.5) +
    scale_colour_manual(values = cols_curve, name = NULL) +
    scale_linetype_manual(values = lty, name = "PP3", labels = paste0("PP3_", levels(lev$level))) +
    coord_cartesian(ylim = c(0, 1), xlim = range(cur$score), expand = FALSE) +
    labs(x = "Score", y = "Posterior probability of pathogenicity", title = paste(title, "- pathogenic")) +
    base + guides(colour = guide_legend(nrow = 1, order = 1), linetype = guide_legend(nrow = 2, order = 2))

  pb <- ggplot() +
    vlines("BP4_discounted") +
    geom_hline(data = lev, aes(yintercept = Post_b, linetype = level), colour = "#1F77B4", linewidth = 0.8) +
    geom_line(data = curve_df("post_b", "post_b_p5")[!is.na(value)], aes(score, value, colour = series), linewidth = 0.5) +
    scale_colour_manual(values = cols_curve, name = NULL) +
    scale_linetype_manual(values = lty, name = "BP4", labels = paste0("BP4_", levels(lev$level))) +
    coord_cartesian(ylim = c(0.975, 1), xlim = range(cur$score), expand = FALSE) +
    labs(x = "Score", y = "Posterior probability of benignity", title = paste(title, "- benign")) +
    base + guides(colour = guide_legend(nrow = 1, order = 1), linetype = guide_legend(nrow = 2, order = 2))

  pp + pb + plot_annotation(caption = "Dashed vertical lines: final thresholds (one-sided 95% bootstrap confidence bound).")
}

save_figure <- function(p, out, width = 11, height = 6) {
  ggsave(paste0(out, ".png"), p, width = width, height = height, dpi = 300, device = ragg::agg_png)
  ggsave(paste0(out, ".pdf"), p, width = width, height = height, device = cairo_pdf)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 3) {
  a <- commandArgs(trailingOnly = TRUE)
  save_figure(plot_posteriors(a[1], a[3]), a[2])
  cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
}
