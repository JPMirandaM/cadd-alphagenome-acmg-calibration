# Paper 2, amendment A-P2-3: evaluation of the primary question H-P and the prespecified secondary outputs.
# Strategies (2024 hold-out calibration, seed 20260925, discounted thresholds): M mechanism, R region (both with the
# pooled genome-wide thresholds as fallback, and secondarily without fallback), G genome-wide, X missense everywhere.
# (a) Temporal hold-out: per threshold group x interval LRs in the test set (2026 variants absent from 2024); evaluable
#     if the group has >= 100 test P/LP and >= 100 test B/LB and the interval holds >= 10 test variants; proportion of
#     evaluable intervals meeting the required LR; bootstrap (1,000, seed 20261101) of differences between strategies.
# (b) Resolved VUS of 2024-01-27 (all classes): (final direction - opposite direction) / resolved, with bootstrap.
# Sensitivity (A-P2-2a): both criteria without the AVI training chromosomes (1, 4, 7, 8, 10, 13, 15).
# Other outputs: H1/H2 replication, conservation-only check, descriptive bins (regulatory stratum and UTR/flank
# classes), seed and prior agreement of the 2026 thresholds, population-control thresholds.
# Output: C:/cadd/calibration/results/paper2/evaluation/
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/paper2_regions.R")
cal <- "C:/cadd/calibration/results/paper2/calibration"; out <- "C:/cadd/calibration/results/paper2/evaluation"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
RB <- if (length(args) >= 1) as.integer(args[1]) else 1000L      # RB < 1000 only for smoke tests
if (RB != 1000L) out <- paste0(out, "_test_RB", RB)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt"); SEED <- 20261101L
phred <- function(cdf) round(-10 * log10(1 - cdf), 2)
train_chr <- c("1", "4", "7", "8", "10", "13", "15")
lev <- c("Supporting", "Moderate", "Strong", "VeryStrong")
req <- setNames(posterior_thresholds()$LR, c("VeryStrong", "Strong", "Moderate", "Supporting"))
log <- function(...) cat(sprintf(...), "\n")
strata <- c("splicing", "protein", "regulatory", "mixed", "conservation_only", "no_predicted_effect")

## Thresholds ------------------------------------------------------------------------------------------
rd_th <- function(a, st, seed = 20260925L) {
  f <- file.path(cal, a, st, paste0("seed", seed), "thresholds.tsv"); if (!file.exists(f)) return(NULL)
  th <- fread(f); list(P = setNames(th$PP3_discounted, th$level), B = setNames(th$BP4_discounted, th$level))
}
ho_sets <- basename(list.dirs(file.path(cal, "holdout"), recursive = FALSE))
TH <- setNames(lapply(ho_sets, function(s) rd_th("holdout", s)), ho_sets)
log("2024 hold-out thresholds available for: %s", paste(ho_sets, collapse = ", "))
interval <- function(s, th) {
  lv <- rep("Indeterminate", length(s))
  for (L in lev) if (!is.na(th$P[L])) lv[s >= th$P[L]] <- paste0("PP3_", L)
  for (L in lev) if (!is.na(th$B[L])) lv[s <= th$B[L]] <- paste0("BP4_", L)
  lv
}
# Threshold group of each variant under a strategy; NA = no evidence (no-fallback variants)
assign_group <- function(d, strategy) switch(strategy,
  M = fifelse(d$stratum %chin% ho_sets, d$stratum, "all"), R = fifelse(d$class %chin% ho_sets, d$class, "all"),
  M_nofb = fifelse(d$stratum %chin% ho_sets, d$stratum, NA_character_), R_nofb = fifelse(d$class %chin% ho_sets, d$class, NA_character_),
  G = rep("all", nrow(d)), X = rep("missense", nrow(d)))
apply_strategy <- function(d, strategy) {
  g <- assign_group(d, strategy); iv <- rep("None", nrow(d))
  for (gr in unique(na.omit(g))) { i <- which(g == gr); iv[i] <- interval(d$s[i], TH[[gr]]) }
  data.table(group = fifelse(is.na(g), "none", fifelse(g == "all" & strategy %chin% c("M", "R"), "fallback_all", g)), iv = iv)
}
strategies <- c("M", "R", "G", "X", "M_nofb", "R_nofb")

## Data: test set, 2026 calibration set ---------------------------------------------------------------
rs <- readRDS("C:/cadd/calibration/results/paper2/refsets/refsets.rds")[, s := phred(avi_cdf)]
test <- rs[set == "test_2026_new"]; c26 <- rs[set == "cal_2026"]
log("Test set: %d P/LP, %d B/LB", test[y == 1, .N], test[y == 0, .N])

## (a) Interval LRs per strategy ------------------------------------------------------------------------
lr_tab <- function(d, a) {   # d: y; a: group, iv
  x <- data.table(y = d$y, group = a$group, iv = a$iv)[!iv %chin% c("Indeterminate", "None")]
  tot <- data.table(y = d$y, group = a$group)[, .(nP = sum(y == 1), nB = sum(y == 0)), by = group]
  z <- x[, .(p = sum(y == 1), b = sum(y == 0)), by = .(group, iv)][tot, on = "group", nomatch = NULL]
  z[, LR := (pmax(p, 0.5) / nP) / (pmax(b, 0.5) / nB)]
  z[, level := sub("^(PP3|BP4)_", "", iv)][, dir := sub("_.*", "", iv)]
  z[, LR_req := fifelse(dir == "PP3", req[level], 1 / req[level])]
  z[, evaluable := nP >= 100 & nB >= 100 & (p + b) >= 10]
  z[, meets := fifelse(dir == "PP3", LR >= LR_req, LR <= LR_req)]
  z
}
summ_a <- function(z) z[evaluable == TRUE, .(n_eval = .N, n_meet = sum(meets), prop = mean(meets),
                                             n_var_evid = sum(p + b), n_var_meet = sum((p + b) * meets))]
eval_a <- function(d, tag) {
  A <- lapply(strategies, function(st) apply_strategy(d, st)); names(A) <- strategies
  Z <- rbindlist(lapply(strategies, function(st) lr_tab(d, A[[st]])[, strategy := st]))
  # Bootstrap: resample variants, recompute LRs, CI per interval and for strategy differences
  set.seed(SEED)
  bs <- rbindlist(lapply(seq_len(RB), function(b) {
    i <- sample.int(nrow(d), replace = TRUE); db <- d[i]
    rbindlist(lapply(strategies, function(st) { ab <- A[[st]][i]; zb <- lr_tab(db, ab)
      cbind(data.table(rep = b, strategy = st), summ_a(zb), data.table(ivs = list(setNames(zb$LR, paste(zb$group, zb$iv))))) }))
  }))
  ci <- rbindlist(lapply(strategies, function(st) {
    m <- rbindlist(lapply(bs[strategy == st, ivs], function(v) as.data.table(as.list(v))), fill = TRUE)
    data.table(strategy = st, ivkey = names(m), lo = sapply(m, quantile, 0.025, na.rm = TRUE), hi = sapply(m, quantile, 0.975, na.rm = TRUE))
  }))
  Z[, ivkey := paste(group, iv)]; Z <- merge(Z, ci, by = c("strategy", "ivkey"), all.x = TRUE)[, ivkey := NULL]
  Z[, meets_bound := fifelse(dir == "PP3", lo >= LR_req, hi <= LR_req)]
  S <- Z[evaluable == TRUE, .(n_eval = .N, n_meet = sum(meets), prop = mean(meets), n_meet_bound = sum(meets_bound, na.rm = TRUE),
                              n_var_evid = sum(p + b), n_var_meet = sum((p + b) * meets)), keyby = strategy]
  S[, prop_bound := n_meet_bound / n_eval][, share_var_meet := n_var_meet / n_var_evid]
  w <- dcast(bs[, .(rep, strategy, prop)], rep ~ strategy, value.var = "prop")
  D <- rbindlist(lapply(c("R", "G", "X"), function(o) data.table(contrast = paste0("M - ", o),
    diff = S[strategy == "M", prop] - S[strategy == o, prop],
    lo = quantile(w$M - w[[o]], 0.025, na.rm = TRUE), hi = quantile(w$M - w[[o]], 0.975, na.rm = TRUE))))
  fwrite(Z[order(strategy, group, dir, match(level, lev))], file.path(out, sprintf("a_intervals_%s.tsv", tag)), sep = "\t")
  fwrite(S, file.path(out, sprintf("a_summary_%s.tsv", tag)), sep = "\t"); fwrite(D, file.path(out, sprintf("a_contrasts_%s.tsv", tag)), sep = "\t")
  cat(sprintf("\n== Criterion (a), %s ==\n", tag)); print(S); print(D)
  list(Z = Z, A = A)
}
ra <- eval_a(test, "all_chromosomes")
invisible(eval_a(test[!chrom %chin% train_chr], "no_training_chromosomes"))

## (b) Resolved VUS --------------------------------------------------------------------------------------
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"); blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")
rd <- function(f) { cv <- fread(f, sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5, 8), col.names = c("chrom", "pos", "ref", "alt", "info"), colClasses = list(character = 1))
  cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
  cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                      grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                      grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
  cv }
cv24 <- rd("C:/cadd/annot/clinvar_archive/clinvar_20240127.vcf.gz")
dg24 <- cv24[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo), unique(sub(":.*", "", unlist(strsplit(geneinfo, "|", fixed = TRUE))))]
vus <- unique(cv24[grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt) & stars >= 1 & clnsig == "Uncertain_significance", ..k]); rm(cv24); invisible(gc())
cv26 <- unique(rd("C:/cadd/annot/clinvar.vcf.gz")[, .(chrom, pos, ref, alt, clnsig_2026 = clnsig)], by = k)
vus <- merge(vus, cv26, by = k, all.x = TRUE); rm(cv26); invisible(gc())
vus[, outcome := fcase(clnsig_2026 %chin% plp, 1L, clnsig_2026 %chin% blb, 0L, default = NA_integer_)]
log("2024 SNV VUS (>= 1 star): %d; resolved P/LP %d, B/LB %d", nrow(vus), vus[outcome == 1, .N], vus[outcome == 0, .N])
vus <- vus[!is.na(outcome)]
# Region class, AF and attributions (rules of script 59)
coding_pos <- unique(readRDS(file.path(CDIR, "mane_cds_positions.rds"))[, .(chrom, pos)])
uc <- vus[coding_pos, on = .(chrom, pos), nomatch = NULL][, ..k]; un <- vus[!coding_pos, on = .(chrom, pos)][, ..k]
edge <- exon_edge_positions()
cc <- rbindlist(lapply(split(uc, uc$chrom), function(x) classify_coding_chr(x, edge)))
nc <- merge(un, classify_noncoding(un, mane_structure()), by = c("chrom", "pos"))[!class %chin% c("exonic_other", "intergenic")]
vus <- merge(vus, rbind(cc, nc[, .(chrom, pos, ref, alt, gene, class)]), by = k)[gene %chin% dg24]
af <- unique(rbind(readRDS(file.path(CDIR, "gnomad_af_mane_coding.rds"))[, .(chrom, pos, ref, alt, AF_max)],
                   rbindlist(lapply(list.files("C:/cadd/calibration/results/paper2/refsets/tmp", "^af_.*\\.rds$", full.names = TRUE), readRDS))[, .(chrom, pos, ref, alt, AF_max)]))
af <- af[, .(AF_max = max(AF_max)), by = k]
vus <- merge(vus, af, by = k, all.x = TRUE)
# Noncoding VUS not in the caches: gnomAD AF from the VEP 116 cache, as in scripts 54 and 59 (absent = 0)
vep <- "C:/cadd/vep/homo_sapiens/116_GRCh38"
vcols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", readLines(file.path(vep, "info.txt")), value = TRUE)), ",")[[1]]
af_chr <- function(x) {
  f <- file.path(vep, x$chrom[1], "all_vars.gz"); if (!file.exists(f)) return(NULL)
  p <- sort(unique(x$pos)); grp <- cumsum(c(1L, diff(p) > 1L))
  rg <- data.table(s = p, g = grp)[, .(start = min(s), end = max(s)), by = g]; tf <- open(TabixFile(f, index = paste0(f, ".csi"))); on.exit(close(tf))
  r <- rbindlist(lapply(split(rg, ceiling(seq_len(nrow(rg)) / 5000)), function(b) {
    l <- unlist(scanTabix(tf, param = GenomicRanges::GRanges(x$chrom[1], IRanges::IRanges(b$start, b$end))), use.names = FALSE)
    if (!length(l)) return(NULL); sel <- which(vcols %in% c("start", "allele_string", "gnomADe", "gnomADg"))
    d <- fread(text = l, sep = "\t", header = FALSE, colClasses = "character", select = sel, col.names = vcols[sel], nThread = 1)
    d <- d[grepl("^[ACGT](/[ACGT])+$", allele_string) & (gnomADe != "." | gnomADg != ".")]; if (!nrow(d)) return(NULL)
    lg <- rbindlist(lapply(c("gnomADe", "gnomADg"), function(kk) { y <- d[get(kk) != ".", .(pos = as.integer(start), allele_string, v = get(kk))]
      if (!nrow(y)) return(NULL); y <- y[, .(e = unlist(strsplit(v, ",", fixed = TRUE))), by = .(pos, allele_string)]
      y[, .(pos, ref = sub("/.*", "", allele_string), alt = sub(":.*", "", e), af = as.numeric(sub(".*:", "", e)))] }))
    lg[, .(AF_max = max(af)), by = .(pos, ref, alt)] }))
  if (is.null(r) || !nrow(r)) return(NULL)
  r[, chrom := x$chrom[1]]
}
miss <- vus[is.na(AF_max), ..k]
if (nrow(miss)) {
  cl <- parallel::makeCluster(6L)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) }))
  parallel::clusterExport(cl, c("vep", "vcols", "af_chr"))
  afm <- rbindlist(parallel::parLapplyLB(cl, split(miss, miss$chrom), af_chr)); parallel::stopCluster(cl)
  if (nrow(afm)) vus[afm, on = k, AF_max := i.AF_max]
}
log("Resolved VUS with a region class in a 2024 disease gene: %d; AF looked up in the VEP cache: %d; absent from gnomAD (AF = 0): %d",
    nrow(vus), nrow(miss), vus[is.na(AF_max), .N])
vus[is.na(AF_max), AF_max := 0]; vus <- vus[AF_max < 0.01]
fi <- rbindlist(lapply(list.files(c("C:/cadd/calibration/results/paper2/atlas_clinvar", "C:/cadd/calibration/results/paper2/atlas_clinvar_batch2"),
                                  "^chunk_.*\\.tsv\\.gz$", full.names = TRUE), fread, sep = "\t", colClasses = list(character = "error")))
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
v <- tstrsplit(sub("^chr", "", fi$variant), "[:>]"); fi[, `:=`(chrom = fifelse(v[[1]] == "M", "MT", v[[1]]), pos = as.integer(v[[2]]), ref = v[[3]], alt = v[[4]])]
n0 <- nrow(vus); vus <- merge(vus, fi[, .(chrom, pos, ref, alt, avi_cdf, stratum)], by = k)[, s := phred(avi_cdf)]
rm(fi); invisible(gc())
log("Resolved rare VUS: %d (P/LP %d, B/LB %d); without Atlas attributions (excluded): %d", nrow(vus), vus[outcome == 1, .N], vus[outcome == 0, .N], n0 - nrow(vus))
fwrite(vus[, .N, by = .(stratum, class, outcome)][order(stratum, class, outcome)], file.path(out, "b_resolved_vus_counts.tsv"), sep = "\t")
eval_b <- function(d, tag) {
  A <- lapply(strategies, function(st) apply_strategy(d, st)); names(A) <- strategies
  sc <- function(iv, y) { pp <- grepl("^PP3", iv); bb <- grepl("^BP4", iv)
    c(final = sum(pp & y == 1) + sum(bb & y == 0), opposite = sum(pp & y == 0) + sum(bb & y == 1)) }
  S <- rbindlist(lapply(strategies, function(st) { x <- sc(A[[st]]$iv, d$y)
    data.table(strategy = st, n = nrow(d), final = x[["final"]], opposite = x[["opposite"]], net = (x[["final"]] - x[["opposite"]]) / nrow(d)) }))
  set.seed(SEED)
  bs <- t(replicate(RB, { i <- sample.int(nrow(d), replace = TRUE)
    sapply(strategies, function(st) { x <- sc(A[[st]]$iv[i], d$y[i]); (x[["final"]] - x[["opposite"]]) / length(i) }) }))
  S[, `:=`(lo = apply(bs, 2, quantile, 0.025), hi = apply(bs, 2, quantile, 0.975))]
  D <- rbindlist(lapply(c("R", "G", "X"), function(o) data.table(contrast = paste0("M - ", o), diff = S[strategy == "M", net] - S[strategy == o, net],
    lo = quantile(bs[, "M"] - bs[, o], 0.025), hi = quantile(bs[, "M"] - bs[, o], 0.975))))
  L <- rbindlist(lapply(strategies, function(st) data.table(strategy = st, iv = A[[st]]$iv, y = d$y)[, .N, by = .(strategy, iv, outcome = fifelse(y == 1, "P/LP", "B/LB"))]))
  fwrite(S, file.path(out, sprintf("b_summary_%s.tsv", tag)), sep = "\t"); fwrite(D, file.path(out, sprintf("b_contrasts_%s.tsv", tag)), sep = "\t")
  fwrite(dcast(L, strategy + iv ~ outcome, value.var = "N", fill = 0), file.path(out, sprintf("b_levels_%s.tsv", tag)), sep = "\t")
  cat(sprintf("\n== Criterion (b), %s ==\n", tag)); print(S); print(D)
}
vus[, y := outcome]
eval_b(vus, "all_chromosomes"); eval_b(vus[!chrom %chin% train_chr], "no_training_chromosomes")

## H1/H2 replication and conservation-only check ----------------------------------------------------------
th26 <- rd_th("main", "splicing")
h1 <- data.table(level = lev, PP3_splicing_2026 = th26$P[lev], BP4_splicing_2026 = th26$B[lev])
fwrite(h1, file.path(out, "h1_splicing_thresholds_2026.tsv"), sep = "\t")
aM <- ra$A$M; h2 <- lr_tab(list(y = test$y), data.table(group = test$class, iv = aM$iv))
fwrite(h2, file.path(out, "h2_strategyM_intervals_by_class.tsv"), sep = "\t")
cat("\n== H2: BP4 intervals of strategy M by region class (evaluable) ==\n"); print(h2[dir == "BP4" & evaluable == TRUE, .(group, iv, b, p, LR, LR_req, meets)])
co <- test[stratum == "conservation_only"]
cons <- lr_tab(co, data.table(group = rep("conservation_only", nrow(co)), iv = interval(co$s, TH[["all"]])))
fwrite(cons, file.path(out, "secondary2_conservation_only.tsv"), sep = "\t")
cat("\n== Conservation-only stratum, genome-wide thresholds, hold-out ==\n"); print(cons)

## Descriptive bins (secondary question 1) ------------------------------------------------------------------
br <- c(0, 10, 15, 20, 25, 30, Inf)
bin_lr <- function(d) { d <- copy(d)[, bin := cut(s, br, right = FALSE, include.lowest = TRUE)]
  f <- function(dd) { nP <- sum(dd$y == 1); nB <- sum(dd$y == 0)
    dd[, .(p = sum(y == 1), b = sum(y == 0)), keyby = bin][, LR := (pmax(p, 0.5) / nP) / (pmax(b, 0.5) / nB)] }
  z <- f(d); set.seed(SEED)
  bb <- sapply(seq_len(RB), function(r) { x <- f(d[sample.int(nrow(d), replace = TRUE)]); x$LR[match(z$bin, x$bin)] })
  z[, `:=`(lo = apply(matrix(bb, nrow = nrow(z)), 1, quantile, 0.025, na.rm = TRUE), hi = apply(matrix(bb, nrow = nrow(z)), 1, quantile, 0.975, na.rm = TRUE))]
  z[, `:=`(nP = nrow(d[y == 1]), nB = nrow(d[y == 0]))] }
desc <- rbind(bin_lr(c26[stratum == "regulatory"])[, set := "stratum_regulatory"],
              rbindlist(lapply(c("5UTR", "3UTR", "upstream_5kb", "downstream_5kb"), function(cl) bin_lr(c26[class == cl])[, set := cl])))
fwrite(desc, file.path(out, "secondary1_descriptive_bins.tsv"), sep = "\t")
cat("\n== Descriptive Phred bins, 2026 set ==\n"); print(desc)

## Seeds, priors and population controls of the 2026 thresholds ------------------------------------------------
main_sets <- basename(list.dirs(file.path(cal, "main"), recursive = FALSE))
grab <- function(a, st, sd) { f <- file.path(cal, a, st, paste0("seed", sd), "thresholds.tsv"); if (!file.exists(f)) return(NULL)
  fread(f)[, .(level, PP3_discounted, BP4_discounted)][, `:=`(analysis = a, set = st, seed = sd)] }
allth <- rbind(rbindlist(lapply(main_sets, function(st) rbindlist(lapply(c(20260925L, 101L, 202L, 303L, 404L), function(sd) grab("main", st, sd))))),
               rbindlist(lapply(c("prior001", "prior010", "popctrl", "holdout"), function(a)
                 rbindlist(lapply(basename(list.dirs(file.path(cal, a), recursive = FALSE)), function(st) grab(a, st, 20260925L))))))
fwrite(allth, file.path(out, "thresholds_all_calibrations.tsv"), sep = "\t")
seedr <- allth[analysis == "main", .(PP3_min = min(PP3_discounted), PP3_max = max(PP3_discounted), BP4_min = min(BP4_discounted), BP4_max = max(BP4_discounted), n_seeds = .N), by = .(set, level)]
fwrite(seedr, file.path(out, "seed_range_2026.tsv"), sep = "\t")
cat("\n== Seed range of the 2026 discounted thresholds ==\n"); print(seedr)
cat("Evaluation finished\n")
