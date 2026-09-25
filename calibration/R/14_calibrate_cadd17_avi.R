# Calibration of CADD v1.7 PHRED and AVI Phred as pre-specified in analysis plan v1
# (STUDY_PLAN.md, Git tag analysis-plan-v1). Uses the validated R port (pejaver_port.R) unchanged;
# adds (a) the gene-level cluster bootstrap for sensitivity analysis S3 and (b) AUROC/AUPRC.
# Runs (resumable; a run is skipped if its thresholds.tsv exists):
#   main : 5 seeds x {CADD, AVI}
#   S1 (>= 2 stars), S2 (AF < 0.001), S3 (gene bootstrap), S4 (no exon-edge): seed 20260925 x {CADD, AVI}
# Output: C:/cadd/calibration/results/v17_avi/<analysis>/<tool>/seed<seed>/
# Usage: Rscript 14_calibrate_cadd17_avi.R [workers] [analyses, comma-separated; default all] [B]
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R")
source("C:/cadd/calibration/R/pejaver_port_fast.R")   # bit-identical to pejaver_port.R (auxiliary A14)
args <- commandArgs(trailingOnly = TRUE)
workers  <- if (length(args) >= 1) as.integer(args[1]) else 10L
analyses <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("main", "S1", "S2", "S3", "S4")
B <- if (length(args) >= 3) as.integer(args[3]) else 10000L      # B < 10000 only for smoke tests
main_seeds <- c(20260925L, 101L, 202L, 303L, 404L)
root <- if (B == 10000L) "C:/cadd/calibration/results/v17_avi" else sprintf("C:/cadd/calibration/results/v17_avi_test_B%d", B)
clin_all <- readRDS("C:/cadd/data/refsets/clinvar_ref.rds")
gno_all  <- readRDS("C:/cadd/data/refsets/gnomad_ref.rds")
tools <- c(CADD = "CADD_PHRED", AVI = "AVI_PHRED")

subset_for <- function(a) switch(a,
  main = list(c = clin_all, g = gno_all),
  S1   = list(c = clin_all[stars >= 2], g = gno_all),
  S2   = list(c = clin_all[AF_max < 0.001], g = gno_all),
  S3   = list(c = clin_all, g = gno_all),
  S4   = list(c = clin_all[exon_edge == FALSE], g = gno_all[exon_edge == FALSE]))

# Gene-level (cluster) bootstrap replicate: resample genes with replacement, keep all their variants
bootstrap_once_cluster <- function(x, y, g, thrs, w, minpoints, gft, increment, gx, gg) {
  ux <- unique(gx); ug <- unique(gg)
  sx <- sample(ux, length(ux), replace = TRUE); sg <- sample(ug, length(ug), replace = TRUE)
  ix <- unlist(split(seq_along(gx), gx)[sx], use.names = FALSE)
  ig <- unlist(split(seq_along(gg), gg)[sg], use.names = FALSE)
  local_posteriors(x[ix], y[ix], g[ig], thrs, w, minpoints, gft, increment)$post
}

# calibrate() from pejaver_port.R with a pluggable bootstrap (identical otherwise)
calibrate_ext <- function(x, y, g, seed, cluster = FALSE, gx = NULL, gg = NULL,
                          alpha = 0.0441, c = 1124, discountonesided = 0.05, minpoints = 100,
                          gft = 0.03, increment = 0.1, pct_band = 5) {
  pt <- posterior_thresholds(alpha, c); w <- negative_weight(y, alpha); thrs <- score_grid(x, g)
  lp <- local_posteriors(x, y, g, thrs, w, minpoints, gft, increment)
  post_p <- lp$post; post_b <- 1 - rev(post_p)
  cl <- parallel::makeCluster(workers); on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, seed)
  parallel::clusterExport(cl, c("local_posteriors", "bootstrap_once", "bootstrap_once_cluster"), envir = globalenv())
  boot_p <- do.call(rbind, parallel::parLapply(cl, seq_len(B), function(i, x, y, g, thrs, w, mp, gf, inc, cluster, gx, gg)
    if (cluster) bootstrap_once_cluster(x, y, g, thrs, w, mp, gf, inc, gx, gg) else bootstrap_once(x, y, g, thrs, w, mp, gf, inc),
    x, y, g, thrs, w, minpoints, gft, increment, cluster, gx, gg))
  boot_b <- 1 - boot_p[, ncol(boot_p):1, drop = FALSE]
  bthr_p <- all_thresholds(boot_p, thrs, pt$Post_p); bthr_b <- all_thresholds(boot_b, rev(thrs), pt$Post_b)
  band <- function(m) apply(m, 2, quantile, probs = pct_band / 100, type = 5, names = FALSE)
  lev <- c("VeryStrong", "Strong", "Moderate", "Supporting")
  list(
    levels = data.frame(level = lev, LR = pt$LR, Post_p = pt$Post_p, Post_b = pt$Post_b),
    thresholds = data.frame(level = lev,
      PP3_point = all_thresholds(post_p, thrs, pt$Post_p)[1, ], PP3_discounted = discounted_thresholds(bthr_p, B, discountonesided, "pathogenic"),
      BP4_point = all_thresholds(post_b, rev(thrs), pt$Post_b)[1, ], BP4_discounted = discounted_thresholds(bthr_b, B, discountonesided, "benign"),
      PP3_boot_NA_frac = colMeans(is.na(bthr_p)), BP4_boot_NA_frac = colMeans(is.na(bthr_b))),
    curves = data.frame(score = thrs, post_p = post_p, post_p_p5 = band(boot_p), post_b = rev(post_b),
                        post_b_p5 = rev(band(boot_b)), halfwindow = lp$finalwindow),
    boot_thresholds = setNames(as.data.frame(cbind(bthr_p, bthr_b)), c(paste0("PP3_", lev), paste0("BP4_", lev))),
    params = list(alpha = alpha, c = c, B = B, seed = seed, cluster_bootstrap = cluster, minpoints = minpoints,
                  gft = gft, increment = increment, n_plp = sum(y == 1), n_blb = sum(y == 0), n_gnomad = length(g), w = w))
}

# AUROC (Mann-Whitney, ties = 0.5) and AUPRC (average precision over distinct thresholds)
auroc <- function(s, y) {   # counts as double: n1 * n0 overflows R integers with ~1e5 variants per class
  r <- rank(s); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0))
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auprc <- function(s, y) {
  d <- data.table(s, y)[, .(tp = sum(y == 1), fp = sum(y == 0)), by = s][order(-s)]
  d[, `:=`(TP = cumsum(tp), FP = cumsum(fp))]; d[, `:=`(prec = TP / (TP + FP), rec = TP / sum(tp))]
  sum(diff(c(0, d$rec)) * d$prec)
}
discrimination <- function(s, y, seed, R = 1000L) {
  set.seed(seed)
  bs <- replicate(R, { i <- sample.int(length(y), replace = TRUE); c(auroc(s[i], y[i]), auprc(s[i], y[i])) })
  data.table(metric = c("AUROC", "AUPRC"), estimate = c(auroc(s, y), auprc(s, y)),
             lower95 = apply(bs, 1, quantile, 0.025), upper95 = apply(bs, 1, quantile, 0.975),
             prevalence = mean(y == 1), n = length(y))
}

for (a in analyses) {
  d <- subset_for(a)
  seeds <- if (a == "main") main_seeds else main_seeds[1]
  for (tn in names(tools)) for (s in seeds) {
    out <- file.path(root, a, tn, paste0("seed", s))
    if (file.exists(file.path(out, "thresholds.tsv"))) { cat("skip", a, tn, s, "\n"); next }
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    x <- d$c[[tools[[tn]]]]; y <- d$c$y; g <- d$g[[tools[[tn]]]]
    t0 <- Sys.time()
    res <- calibrate_fast(x, y, g, B = B, workers = workers, seed = s,
                          cluster = (a == "S3"), gx = d$c$gene, gg = d$g$gene)
    fwrite(res$thresholds, file.path(out, "thresholds.tsv"), sep = "\t")
    fwrite(res$levels, file.path(out, "posterior_levels.tsv"), sep = "\t")
    fwrite(res$curves, file.path(out, "posterior_curves.tsv"), sep = "\t")
    fwrite(res$boot_thresholds, file.path(out, "bootstrap_thresholds.tsv"), sep = "\t")
    writeLines(capture.output(str(res$params)), file.path(out, "parameters.txt"))
    if (a == "main" && s == main_seeds[1]) fwrite(discrimination(x, y, s), file.path(out, "discrimination.tsv"), sep = "\t")
    cat(sprintf("[%s] %s %s seed %d: P/LP %d, B/LB %d, gnomAD %d, %.1f min\n", format(Sys.time(), "%H:%M"), a, tn, s,
                sum(y == 1), sum(y == 0), length(g), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    print(res$thresholds, digits = 4)
  }
}

# Summary: main analysis median and range over seeds; sensitivity analyses
all <- rbindlist(lapply(list.files(root, "thresholds.tsv", recursive = TRUE, full.names = TRUE), function(f) {
  p <- strsplit(sub(paste0(root, "/"), "", f), "/")[[1]]
  num <- c("PP3_point", "PP3_discounted", "BP4_point", "BP4_discounted", "PP3_boot_NA_frac", "BP4_boot_NA_frac")
  data.table(analysis = p[1], tool = p[2], seed = sub("seed", "", p[3]),
             fread(f, colClasses = c(level = "character", setNames(rep("numeric", length(num)), num))))
}), fill = TRUE)
fwrite(all, file.path(root, "thresholds_all_runs.tsv"), sep = "\t")
cat("\nMain analysis, discounted thresholds (median [min-max] over seeds):\n")
long <- melt(all[analysis == "main"], id.vars = c("tool", "seed", "level"), measure.vars = c("PP3_discounted", "BP4_discounted"))
print(long[, .(median = median(value), min = min(value), max = max(value), n_NA = sum(is.na(value))), by = .(tool, variable, level)])
