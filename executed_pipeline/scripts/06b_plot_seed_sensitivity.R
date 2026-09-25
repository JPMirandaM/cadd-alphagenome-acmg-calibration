# Figure S2: seed (Monte Carlo) sensitivity of the bootstrap-based BP4 Strong threshold.
# Built ONLY from the files in <data_dir>:
#   bootstrap_thresholds_seed<seed>.tsv  one row per bootstrap replicate (column BP4_Strong)
#   thresholds_all_seeds.tsv             final (discounted) thresholds per seed
# Usage: Rscript plot_seed_sensitivity.R <data_dir> <output_path_without_extension>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

plot_seed_sensitivity <- function(datadir, xmax = 1, main_seed = "20260925") {
  files <- list.files(datadir, "^bootstrap_thresholds_seed\\d+\\.tsv$", full.names = TRUE)
  boot <- rbindlist(lapply(files, function(f)
    data.table(seed = sub("^bootstrap_thresholds_seed(\\d+)\\.tsv$", "\\1", basename(f)),
               value = fread(f)$BP4_Strong)))
  fin <- fread(file.path(datadir, "thresholds_all_seeds.tsv"))[level == "Strong",
                 .(seed = as.character(seed), BP4_discounted, BP4_boot_NA_frac)]
  boot <- merge(boot, fin, by = "seed")
  boot[, seed_lab := sprintf("seed %s%s (%.1f%% of replicates fail; bound = %.3f)",
                             seed, ifelse(seed == main_seed, " (main run, Fig. S1)", ""),
                             100 * BP4_boot_NA_frac, BP4_discounted)]
  # Legend order: main run first, then the additional seeds in numeric order
  ord <- unique(boot[order(seed != main_seed, as.numeric(seed)), seed_lab])
  boot[, seed_lab := factor(seed_lab, levels = ord)]
  # Cumulative fraction over ALL replicates, with failed replicates (level not reached) placed
  # first, i.e. treated as the most conservative outcome, exactly as the one-sided rule of
  # get_discounted_thresholds.m does. The curve therefore starts at the failure fraction and
  # crosses 5% at the reported bound.
  ecdf_dt <- boot[, {
    v <- sort(value[!is.na(value)]); n <- .N; n_na <- sum(is.na(value))
    list(x = c(0, v), F = c(n_na, n_na + seq_along(v)) / n)
  }, by = seed_lab]
  ggplot(ecdf_dt[x <= xmax], aes(x, F, colour = seed_lab)) +
    geom_step(linewidth = 0.5) +
    geom_hline(aes(yintercept = 0.05), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = 0.15, linetype = "dotted", colour = "black") +
    annotate("text", x = 0.15, y = 0.015, label = "published 0.15", angle = 90, vjust = -0.4, size = 3) +
    scale_colour_brewer(palette = "Dark2", name = NULL) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 0.10)) +
    labs(x = "BP4 Strong threshold in bootstrap replicate (CADD v1.6 PHRED)",
         y = "Cumulative fraction of all replicates",
         title = "Seed sensitivity of the one-sided BP4 Strong bound (B = 10,000 per seed)",
         caption = paste0("Each curve starts at the fraction of replicates that do not reach the level; ",
                          "the bound is where it crosses 5% (dashed line).\n",
                          "Seeds are arbitrary integers fixed in advance: ", main_seed,
                          " (the run date) for the main run, and 101, 202, 303, 404 for the repetitions.")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", legend.direction = "vertical", panel.grid.minor = element_blank())
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 2) {
  a <- commandArgs(trailingOnly = TRUE)
  p <- plot_seed_sensitivity(a[1])
  ggsave(paste0(a[2], ".png"), p, width = 8, height = 6, dpi = 300, device = ragg::agg_png)
  ggsave(paste0(a[2], ".pdf"), p, width = 8, height = 6, device = cairo_pdf)
  cat("Figure:", paste0(a[2], c(".png", ".pdf")), sep = "\n")
}
