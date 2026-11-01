# Paper 2, amendment A-P2-4: functional validation (SGE, MPRA), clinical outcomes on current VUS, and the exploratory
# hybrid strategy. Thresholds: 2026 calibration (seed 20260925, discounted) for functional data and current VUS; the
# 2024 hold-out calibration for the exploratory re-evaluation of criteria (a) and (b) of A-P2-3.
# Strategies: M mechanism, R region, G genome-wide, X missense everywhere, M/R without fallback; exploratory H (stratum
# where calibrated, else class where calibrated, else no evidence) and H_fb (else pooled thresholds).
# Usage: Rscript 63_paper2_functional_evaluate.R [RB]   (RB < 1000 only for smoke tests)
# Output: C:/cadd/calibration/results/paper2/functional/
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/paper2_regions.R")
cal <- "C:/cadd/calibration/results/paper2/calibration"; out <- "C:/cadd/calibration/results/paper2/functional"
args <- commandArgs(trailingOnly = TRUE)
RB <- if (length(args) >= 1) as.integer(args[1]) else 1000L
odir <- if (RB == 1000L) out else paste0(out, "_test_RB", RB); dir.create(odir, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt"); SEED <- 20261101L
phred <- function(cdf) round(-10 * log10(1 - cdf), 2)
train_chr <- c("1", "4", "7", "8", "10", "13", "15"); valid_chr <- c("2", "5", "11", "14", "17", "20", "22", "X")
lev <- c("Supporting", "Moderate", "Strong", "VeryStrong")
req <- setNames(posterior_thresholds()$LR, c("VeryStrong", "Strong", "Moderate", "Supporting"))
log <- function(...) cat(sprintf(...), "\n")
strata <- c("splicing", "protein", "regulatory", "mixed", "conservation_only", "no_predicted_effect")

## Thresholds and strategies --------------------------------------------------------------------------------
rd_th <- function(a, st, seed = 20260925L) {
  f <- file.path(cal, a, st, paste0("seed", seed), "thresholds.tsv"); if (!file.exists(f)) return(NULL)
  th <- fread(f); list(P = setNames(th$PP3_discounted, th$level), B = setNames(th$BP4_discounted, th$level))
}
load_th <- function(a) { s <- basename(list.dirs(file.path(cal, a), recursive = FALSE)); setNames(lapply(s, function(x) rd_th(a, x)), s) }
interval <- function(s, th) {
  lv <- rep("Indeterminate", length(s))
  for (L in lev) if (!is.na(th$P[L])) lv[s >= th$P[L]] <- paste0("PP3_", L)
  for (L in lev) if (!is.na(th$B[L])) lv[s <= th$B[L]] <- paste0("BP4_", L)
  lv
}
assign_group <- function(d, strategy, sets) {
  st_ok <- d$stratum %chin% intersect(sets, strata); cl_ok <- d$class %chin% setdiff(sets, c(strata, "all"))
  switch(strategy,
    M = fifelse(st_ok, d$stratum, "all"), R = fifelse(cl_ok, d$class, "all"),
    M_nofb = fifelse(st_ok, d$stratum, NA_character_), R_nofb = fifelse(cl_ok, d$class, NA_character_),
    G = rep("all", nrow(d)), X = rep("missense", nrow(d)),
    H = fifelse(st_ok, d$stratum, fifelse(cl_ok, d$class, NA_character_)),
    H_fb = fifelse(st_ok, d$stratum, fifelse(cl_ok, d$class, "all")))
}
apply_strategy <- function(d, strategy, TH) {
  g <- assign_group(d, strategy, names(TH)); iv <- rep("None", nrow(d))
  for (gr in unique(na.omit(g))) { i <- which(g == gr); iv[i] <- interval(d$s[i], TH[[gr]]) }
  data.table(group = fifelse(is.na(g), "none", fifelse(g == "all" & strategy %chin% c("M", "R", "H_fb"), "fallback_all", g)), iv = iv)
}
strategies <- c("M", "R", "G", "X", "M_nofb", "R_nofb", "H", "H_fb")
TH26 <- load_th("main"); THho <- load_th("holdout")
log("2026 threshold sets: %s; 2024 hold-out sets: %s", paste(names(TH26), collapse = ", "), paste(names(THho), collapse = ", "))

## Attributions and strata (rule of script 59) ------------------------------------------------------------------
fi <- rbindlist(lapply(list.files(file.path("C:/cadd/calibration/results/paper2", c("atlas_clinvar", "atlas_clinvar_batch2", "atlas_functional")),
                                  "^chunk_.*\\.tsv\\.gz$", full.names = TRUE), fread, sep = "\t", colClasses = list(character = "error")))
fi <- unique(fi[error == "" | is.na(error)], by = "variant")
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
fi <- fi[, .(chrom, pos, ref, alt, avi_cdf, stratum)]; rm(v); invisible(gc())
attach_fi <- function(d, what) { n0 <- nrow(d); d <- merge(d, fi, by = k)[, s := phred(avi_cdf)]
  log("%s: %d of %d SNVs with Atlas attributions", what, nrow(d), n0); d }

## Shared evaluation helpers -----------------------------------------------------------------------------------
net_eval <- function(d, TH, tag, sts = strategies) {   # y = 1 abnormal / functional / P/LP
  A <- setNames(lapply(sts, function(st) apply_strategy(d, st, TH)), sts)
  sc <- function(iv, y) { pp <- grepl("^PP3", iv); bb <- grepl("^BP4", iv)
    c(final = sum(pp & y == 1) + sum(bb & y == 0), opposite = sum(pp & y == 0) + sum(bb & y == 1)) }
  S <- rbindlist(lapply(sts, function(st) { x <- sc(A[[st]]$iv, d$y)
    data.table(strategy = st, n = nrow(d), agree = x[["final"]], opposite = x[["opposite"]], net = (x[["final"]] - x[["opposite"]]) / nrow(d)) }))
  set.seed(SEED)
  bs <- t(replicate(RB, { i <- sample.int(nrow(d), replace = TRUE)
    sapply(sts, function(st) { x <- sc(A[[st]]$iv[i], d$y[i]); (x[["final"]] - x[["opposite"]]) / length(i) }) }))
  S[, `:=`(lo = apply(bs, 2, quantile, 0.025), hi = apply(bs, 2, quantile, 0.975))]
  D <- rbindlist(lapply(c("R", "G", "X"), function(o) data.table(contrast = paste0("M - ", o), diff = S[strategy == "M", net] - S[strategy == o, net],
    lo = quantile(bs[, "M"] - bs[, o], 0.025), hi = quantile(bs[, "M"] - bs[, o], 0.975))))
  fwrite(S, file.path(odir, sprintf("%s_net.tsv", tag)), sep = "\t"); fwrite(D, file.path(odir, sprintf("%s_contrasts.tsv", tag)), sep = "\t")
  cat(sprintf("\n== %s ==\n", tag)); print(S, digits = 3); print(D, digits = 3)
  invisible(A)
}
lr_func <- function(d, A, tag, min_group = 30L) {   # LR per group x interval with direction check
  z <- rbindlist(lapply(names(A), function(st) {
    x <- data.table(y = d$y, group = A[[st]]$group, iv = A[[st]]$iv)
    tot <- x[, .(nA = sum(y == 1), nN = sum(y == 0)), by = group]
    x[!iv %chin% c("Indeterminate", "None"), .(a = sum(y == 1), n = sum(y == 0)), by = .(group, iv)][tot, on = "group", nomatch = NULL][, strategy := st] }))
  z[, LR := (pmax(a, 0.5) / nA) / (pmax(n, 0.5) / nN)][, dir := sub("_.*", "", iv)][, level := sub("^(PP3|BP4)_", "", iv)]
  z[, evaluable := nA >= min_group & nN >= min_group & (a + n) >= 10]
  z[, direction_ok := fifelse(dir == "PP3", LR > 1, LR < 1)][, meets_clingen := fifelse(dir == "PP3", LR >= req[level], LR <= 1 / req[level])]
  setcolorder(z, c("strategy", "group", "iv")); fwrite(z[order(strategy, group, iv)], file.path(odir, sprintf("%s_intervals.tsv", tag)), sep = "\t")
  s <- z[evaluable == TRUE, .(n_eval = .N, n_direction_ok = sum(direction_ok), n_meets_clingen = sum(meets_clingen)), by = strategy]
  fwrite(s, file.path(odir, sprintf("%s_intervals_summary.tsv", tag)), sep = "\t"); cat(sprintf("\n-- %s intervals --\n", tag)); print(s)
}

## 1-2. Functional validation ----------------------------------------------------------------------------------------
c26 <- readRDS("C:/cadd/calibration/results/paper2/refsets/refsets.rds")[set == "cal_2026", ..k]
for (nm in c("sge", "mpra")) {
  d <- attach_fi(readRDS(file.path(out, sprintf("%s_labels.rds", nm))), toupper(nm))
  if (nm == "mpra") d[is.na(class), class := "intergenic"]
  log("%s analysis set: %d (y = 1: %d, y = 0: %d)", toupper(nm), nrow(d), d[y == 1, .N], d[y == 0, .N])
  print(dcast(d[, .N, by = .(stratum, y)], stratum ~ y, value.var = "N", fill = 0))
  A <- net_eval(d, TH26, nm); lr_func(d, A, nm)
  A1 <- net_eval(d[!c26, on = k], TH26, paste0(nm, "_sens_no_clinvar_ref")); lr_func(d[!c26, on = k], A1, paste0(nm, "_sens_no_clinvar_ref"))
  d2 <- d[!chrom %chin% c(train_chr, valid_chr)]
  log("%s without AVI training and validation chromosomes: %d SNVs (y = 1: %d)", toupper(nm), nrow(d2), d2[y == 1, .N])
  if (d2[y == 1, .N] >= 30 && d2[y == 0, .N] >= 30) { A2 <- net_eval(d2, TH26, paste0(nm, "_sens_no_train_valid_chr")); lr_func(d2, A2, paste0(nm, "_sens_no_train_valid_chr")) }
  if (nm == "sge") {   # per stratum and per gene (net agreement, no bootstrap)
    per <- function(by) rbindlist(lapply(strategies, function(st) data.table(d[, c(by, "y"), with = FALSE], iv = A[[st]]$iv, strategy = st)))[
      , .(n = .N, agree = sum((grepl("^PP3", iv) & y == 1) | (grepl("^BP4", iv) & y == 0)),
          opposite = sum((grepl("^PP3", iv) & y == 0) | (grepl("^BP4", iv) & y == 1))), by = c("strategy", by)][, net := (agree - opposite) / n]
    fwrite(per("stratum"), file.path(odir, "sge_by_stratum.tsv"), sep = "\t"); fwrite(per("set_gene"), file.path(odir, "sge_by_gene.tsv"), sep = "\t")
  }
  if (nm == "mpra") {
    per <- rbindlist(lapply(strategies, function(st) data.table(d[, .(element, y)], iv = A[[st]]$iv, strategy = st)))[
      , .(n = .N, agree = sum((grepl("^PP3", iv) & y == 1) | (grepl("^BP4", iv) & y == 0)),
          opposite = sum((grepl("^PP3", iv) & y == 0) | (grepl("^BP4", iv) & y == 1))), by = .(strategy, element)][, net := (agree - opposite) / n]
    fwrite(per, file.path(odir, "mpra_by_element.tsv"), sep = "\t")
  }
}

## 3. Clinical outcomes: current VUS -----------------------------------------------------------------------------------
vus <- attach_fi(readRDS(file.path(out, "vus_current.rds")), "Current rare VUS")
cnt <- rbindlist(lapply(c("M", "M_nofb", "R", "G", "H"), function(st) data.table(vus[, .(stratum, class, ep = stars >= 3)], strategy = st, iv = apply_strategy(vus, st, TH26)$iv)))
fwrite(cnt[, .N, by = .(strategy, iv)][order(strategy, iv)], file.path(odir, "vus_current_levels.tsv"), sep = "\t")
fwrite(cnt[ep == TRUE, .N, by = .(strategy, iv)][order(strategy, iv)], file.path(odir, "vus_current_levels_expert_panel.tsv"), sep = "\t")
fwrite(cnt[strategy %chin% c("M", "M_nofb"), .N, by = .(strategy, stratum, class, iv)][order(strategy, stratum, class, iv)], file.path(odir, "vus_current_levels_by_stratum_class.tsv"), sep = "\t")
cat("\n== Current VUS: evidence per strategy ==\n"); print(dcast(cnt[, .N, by = .(strategy, dir = fifelse(grepl("^PP3", iv), "PP3", fifelse(grepl("^BP4", iv), "BP4", "none")))], strategy ~ dir, value.var = "N", fill = 0))
vus[, ivM := apply_strategy(vus, "M", TH26)$iv]
cand <- vus[stratum == "splicing" & ivM %chin% c("PP3_Moderate", "PP3_Strong", "PP3_VeryStrong")]
log("Splicing-stratum current VUS with >= PP3 Moderate (RNA-study candidates): %d in %d genes; expert panel: %d", nrow(cand), uniqueN(cand$gene), cand[stars >= 3, .N])
fwrite(cand[, .(n = .N, n_expert_panel = sum(stars >= 3), n_strong = sum(ivM != "PP3_Moderate")), by = gene][order(-n)], file.path(odir, "rna_candidates_by_gene.tsv"), sep = "\t")
fwrite(cand[stars >= 3, .(variant = sprintf("chr%s:%d:%s>%s", chrom, pos, ref, alt), gene, class, avi_phred = s, level = ivM)],
       file.path(out, "rna_candidates_expert_panel.tsv"), sep = "\t")   # input of script 64 (Atlas splicing tracks)

## 4. Exploratory: hybrid strategy on criteria (a) and (b) of A-P2-3 (2024 thresholds) -------------------------------
rs <- readRDS("C:/cadd/calibration/results/paper2/refsets/refsets.rds")[set == "test_2026_new"][, s := phred(avi_cdf)]
ex <- rbindlist(lapply(c("M", "R", "G", "H", "H_fb"), function(st) { a <- apply_strategy(rs, st, THho)
  x <- data.table(y = rs$y, group = a$group, iv = a$iv); tot <- x[, .(nP = sum(y == 1), nB = sum(y == 0)), by = group]
  z <- x[!iv %chin% c("Indeterminate", "None"), .(p = sum(y == 1), b = sum(y == 0)), by = .(group, iv)][tot, on = "group", nomatch = NULL]
  z[, LR := (pmax(p, 0.5) / nP) / (pmax(b, 0.5) / nB)][, level := sub("^(PP3|BP4)_", "", iv)][, dir := sub("_.*", "", iv)]
  z[, evaluable := nP >= 100 & nB >= 100 & (p + b) >= 10][, meets := fifelse(dir == "PP3", LR >= req[level], LR <= 1 / req[level])]
  z[evaluable == TRUE, .(strategy = st, n_eval = .N, n_meet = sum(meets), prop = mean(meets), n_var_evid = sum(p + b), n_var_meet = sum((p + b) * meets))] }))
fwrite(ex, file.path(odir, "exploratory_hybrid_criterion_a.tsv"), sep = "\t"); cat("\n== Exploratory hybrid, criterion (a) ==\n"); print(ex)
log("Resolved VUS for the exploratory criterion (b) are rebuilt from the 2024 and 2026 releases (rules of script 61)")
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"); blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")
rdv <- function(f) { cv <- fread(f, sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5, 8), col.names = c("chrom", "pos", "ref", "alt", "info"), colClasses = list(character = 1))
  cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
  cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                      grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                      grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
  cv }
cv24 <- rdv("C:/cadd/annot/clinvar_archive/clinvar_20240127.vcf.gz")
dg24 <- cv24[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo), unique(sub(":.*", "", unlist(strsplit(geneinfo, "|", fixed = TRUE))))]
rv <- unique(cv24[grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt) & stars >= 1 & clnsig == "Uncertain_significance", ..k]); rm(cv24); invisible(gc())
cv26 <- unique(rdv("C:/cadd/annot/clinvar.vcf.gz")[, .(chrom, pos, ref, alt, clnsig_2026 = clnsig)], by = k)
rv <- merge(rv, cv26, by = k)[, y := fcase(clnsig_2026 %chin% plp, 1L, clnsig_2026 %chin% blb, 0L, default = NA_integer_)][!is.na(y)]; rm(cv26); invisible(gc())
coding_pos <- unique(readRDS(file.path(CDIR, "mane_cds_positions.rds"))[, .(chrom, pos)]); edge <- exon_edge_positions()
uc <- rv[coding_pos, on = .(chrom, pos), nomatch = NULL][, ..k]; un <- rv[!coding_pos, on = .(chrom, pos)][, ..k]
cc <- rbindlist(lapply(split(uc, uc$chrom), function(x) classify_coding_chr(x, edge)))
nc <- merge(un, classify_noncoding(un, mane_structure()), by = c("chrom", "pos"))[!class %chin% c("exonic_other", "intergenic")]
rv <- merge(rv, rbind(cc, nc[, .(chrom, pos, ref, alt, gene, class)]), by = k)[gene %chin% dg24]
af <- unique(rbind(readRDS(file.path(CDIR, "gnomad_af_mane_coding.rds"))[, .(chrom, pos, ref, alt, AF_max)],
                   rbindlist(lapply(list.files("C:/cadd/calibration/results/paper2/refsets/tmp", "^af_.*\\.rds$", full.names = TRUE), readRDS))[, .(chrom, pos, ref, alt, AF_max)],
                   readRDS(file.path(out, "vus_current.rds"))[, .(chrom, pos, ref, alt, AF_max)]))[, .(AF_max = max(AF_max)), by = k]
rv <- merge(rv, af, by = k, all.x = TRUE)
miss <- rv[is.na(AF_max), ..k]
log("Resolved VUS rebuilt: %d; AF looked up in the VEP cache: %d", nrow(rv), nrow(miss))
vep <- "C:/cadd/vep/homo_sapiens/116_GRCh38"   # gnomAD AF from the VEP 116 cache, as in scripts 61 and 62 (absent = 0)
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
if (nrow(miss)) {
  cl <- parallel::makeCluster(6L)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) }))
  parallel::clusterExport(cl, c("vep", "vcols", "af_chr"))
  afm <- rbindlist(parallel::parLapplyLB(cl, split(miss, miss$chrom), af_chr)); parallel::stopCluster(cl)
  if (nrow(afm)) rv[afm, on = k, AF_max := i.AF_max]
}
rv[is.na(AF_max), AF_max := 0]
rv <- attach_fi(rv[AF_max < 0.01], "Resolved rare VUS (exploratory)")
log("Check against script 61 (14,513 resolved VUS): %d", nrow(rv))
net_eval(rv, THho, "exploratory_hybrid_criterion_b", c("M", "R", "G", "X", "H", "H_fb"))
log("Evaluation finished")
