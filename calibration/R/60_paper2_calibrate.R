# Paper 2 (plan v1 "Calibration"; amendment A-P2-1): ClinGen local posterior calibration of AVI Phred per
# mechanism stratum, per region class and for all classes pooled (the single genome-wide AVI set of H3/H-P).
# Same R port and parameters as the main article (pejaver_port_fast.R; B = 10,000; minpoints 100; gnomAD
# fraction 0.03; increment 0.1; one-sided 95% discount). A set is calibrated when its labelled data hold
# >= 1,000 P/LP and >= 1,000 benign variants (criterion of A10).
# AVI Phred = -10 log10(1 - AVI CDF quantile from the Atlas API), rounded to 2 decimals as in script 30; it agrees
# with the local Atlas file (max |difference| 0.09, 99th percentile 0.01, 137,102 SNVs of the 2026 set).
# Unlabelled set for window sizing (gft): the rare gnomAD sample (script 54) in the same stratum or class.
# Analyses (resumable: a run is skipped when its thresholds.tsv exists), in this order:
#   main     : 2026-09-23 set, seed 20260925, then seeds 101, 202, 303, 404 (plan: five seeds)
#   holdout  : 2024-01-27 set, seed 20260925 (thresholds for the temporal hold-out, as in A2)
#   popctrl  : 2026 P/LP against the gnomAD sample as the benign class (A-P2-1 secondary question 3)
#   prior001, prior010 : 2026 set at priors 0.01 and 0.10 (plan v1 sensitivity analyses)
# Output: C:/cadd/calibration/results/paper2/calibration/<analysis>/<set>/seed<seed>/
# Usage: Rscript 60_paper2_calibrate.R [workers] [B]
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/pejaver_port_fast.R")
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 20L
B <- if (length(args) >= 2) as.integer(args[2]) else 10000L
root <- if (B == 10000L) "C:/cadd/calibration/results/paper2/calibration" else sprintf("C:/cadd/calibration/results/paper2/calibration_test_B%d", B)
k <- c("chrom", "pos", "ref", "alt")
phred <- function(cdf) round(-10 * log10(1 - cdf), 2)

## Labelled sets (script 59) and the gnomAD sample with strata (same rule as script 59) ------------------
rs <- readRDS("C:/cadd/calibration/results/paper2/refsets/refsets.rds")[, s := phred(avi_cdf)]
gs <- readRDS("C:/cadd/calibration/results/paper2/universe/gnomad_sample.rds")
fi <- rbindlist(lapply(list.files("C:/cadd/calibration/results/paper2/atlas_gnomad", "^chunk_.*\\.tsv\\.gz$", full.names = TRUE),
                       fread, sep = "\t", colClasses = list(character = "error")))
fi <- fi[error == "" | is.na(error)]
pos_sum <- function(cols) Reduce(`+`, lapply(cols, function(cn) pmax(fi[[cn]], 0)))
fi[, `:=`(g_spl = pos_sum("fi_MERGED_SPLICING"), g_pro = pos_sum(c("fi_ALPHAMISSENSE", "fi_PROTEIN_TERMINATION", "fi_START_LOST", "fi_STOP_LOST")),
          g_reg = pos_sum(grep("^fi_MAX_ABS_", names(fi), value = TRUE)), g_con = pos_sum(c("fi_CACTUS_241_WAY", "fi_PHASTCONS_470_WAY")))]
fi[, mech := g_spl + g_pro + g_reg]
fi[, stratum := fifelse(mech >= 0.05,
                        fifelse(g_spl >= 0.5 * mech & g_spl >= g_pro & g_spl >= g_reg, "splicing",
                        fifelse(g_pro >= 0.5 * mech & g_pro >= g_reg, "protein",
                        fifelse(g_reg >= 0.5 * mech, "regulatory", "mixed"))),
                        fifelse(g_con >= 0.05, "conservation_only", "no_predicted_effect"))]
v <- tstrsplit(sub("^chr", "", fi$variant), "[:>]")
fi[, `:=`(chrom = fifelse(v[[1]] == "M", "MT", v[[1]]), pos = as.integer(v[[2]]), ref = v[[3]], alt = v[[4]])]
gs <- merge(gs, fi[, .(chrom, pos, ref, alt, avi_cdf, stratum)], by = k)[, s := phred(avi_cdf)]
cat(sprintf("gnomAD sample with attributions: %d of %d\n", nrow(gs), nrow(readRDS("C:/cadd/calibration/results/paper2/universe/gnomad_sample.rds"))))
print(dcast(gs[, .N, by = .(stratum)], . ~ stratum, value.var = "N"))

## Sets -------------------------------------------------------------------------------------------------
strata <- c("splicing", "protein", "regulatory", "mixed", "conservation_only", "no_predicted_effect")
classes <- c("canonical_splice", "splice_region_intronic", "5UTR", "3UTR", "near_intronic", "deep_intronic", "upstream_5kb",
             "downstream_5kb", "missense", "synonymous_exon_edge", "synonymous_other", "stop_gained", "start_lost", "stop_lost")
pick <- function(d, set_name) if (set_name == "all") d else if (set_name %in% strata) d[stratum == set_name] else d[class == set_name]
build <- function(analysis, set_name) {
  rel <- if (analysis == "holdout") "cal_2024" else "cal_2026"
  lab <- pick(rs[set == rel], set_name); g <- pick(gs, set_name)$s
  if (analysis == "popctrl") { x <- c(lab[y == 1, s], g); y <- c(rep(1L, lab[y == 1, .N]), rep(0L, length(g))) }
  else { x <- lab$s; y <- lab$y }
  list(x = x, y = y, g = g, gene = if (analysis == "popctrl") NULL else lab$gene)
}
prior_of <- function(a) switch(a, prior001 = 0.01, prior010 = 0.10, 0.0441)

## Runs -------------------------------------------------------------------------------------------------
sets <- c(strata, classes, "all")
runs <- rbind(data.table(analysis = "main", seed = 20260925L), data.table(analysis = "holdout", seed = 20260925L),
              data.table(analysis = "main", seed = c(101L, 202L, 303L, 404L)), data.table(analysis = "popctrl", seed = 20260925L),
              data.table(analysis = c("prior001", "prior010"), seed = 20260925L))
sizes <- rbindlist(lapply(c("main", "holdout", "popctrl"), function(a) rbindlist(lapply(sets, function(st) {
  d <- build(a, st); data.table(analysis = a, set = st, n_plp = sum(d$y == 1), n_benign = sum(d$y == 0), n_gnomad = length(d$g)) }))))
sizes[, calibrate := n_plp >= 1000 & n_benign >= 1000 & n_gnomad > 0]
dir.create(root, recursive = TRUE, showWarnings = FALSE); fwrite(sizes, file.path(root, "set_sizes.tsv"), sep = "\t")
cat("\n== Sets calibrated ==\n"); print(sizes[calibrate == TRUE])
for (i in seq_len(nrow(runs))) {
  a <- runs$analysis[i]; sd <- runs$seed[i]
  todo <- sizes[analysis == fifelse(a %in% c("prior001", "prior010"), "main", a) & calibrate == TRUE, set]
  for (st in todo) {
    od <- file.path(root, a, st, paste0("seed", sd))
    if (file.exists(file.path(od, "thresholds.tsv"))) next
    dir.create(od, recursive = TRUE, showWarnings = FALSE); d <- build(a, st); t0 <- Sys.time()
    r <- calibrate_fast(d$x, d$y, d$g, alpha = prior_of(a), B = B, workers = workers, seed = sd)
    fwrite(r$thresholds, file.path(od, "thresholds.tsv"), sep = "\t"); fwrite(r$curves, file.path(od, "curves.tsv"), sep = "\t")
    fwrite(r$boot_thresholds, file.path(od, "boot_thresholds.tsv.gz"), sep = "\t"); saveRDS(r$params, file.path(od, "params.rds"))
    cat(sprintf("[%s] %s / %s / seed %d: %d P/LP, %d benign, %d gnomAD; %.1f min\n", format(Sys.time(), "%m-%d %H:%M"), a, st, sd,
                sum(d$y == 1), sum(d$y == 0), length(d$g), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
cat("All runs finished\n")
