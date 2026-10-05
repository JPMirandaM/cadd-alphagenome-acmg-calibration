# Amendment A15 (Aim 7b): prospective test of calibrated PP3/BP4 evidence on ClinVar VUS of the
# 2024-01-27 release that were resolved by the 2026-09-23 release.
# VUS: CLNSIG == "Uncertain_significance", >= 1 star, SNV, MANE Select missense (same exclusions as
# the reference sets), disease gene (>= 1 P/LP record with >= 1 star in the 2024 release).
# Outcome (2026): resolved to P/LP, resolved to B/LB, or not resolved.
# Primary thresholds: 2024 release series (A6, set r2024, five predictors). Secondary: 2026 main
# (CADD, AVI) and 2026 intersection (five predictors). Primary subset: gnomAD v4.1 AF_max < 0.01.
# Output: C:/cadd/calibration/results/v17_avi/vus_resolved/
suppressPackageStartupMessages({ library(data.table); library(rtracklayer) })
set.seed(20261004)
cdir <- "C:/cadd/data/coding"; cd <- "C:/cadd/comparators"; out <- "C:/cadd/calibration/results/v17_avi/vus_resolved"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")

## ClinVar parsing (rules of scripts 13 and 23) ----------------------------------------------------
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
read_cv <- function(f) {
  cv <- fread(f, sep = "\t", skip = "#CHROM", select = c(1, 2, 3, 4, 5, 8),
              col.names = c("chrom", "pos", "clinvar_id", "ref", "alt", "info"), colClasses = list(character = c(1, 3)))
  cv <- cv[nchar(ref) == 1 & nchar(alt) == 1]
  cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
  cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                      grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                      grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
  cv
}
cv24 <- read_cv("C:/cadd/annot/clinvar_archive/clinvar_20240127.vcf.gz")
dg <- cv24[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo),
           unique(unlist(lapply(strsplit(geneinfo, "|", fixed = TRUE), function(z) sub(":.*", "", z))))]
vus <- cv24[clnsig == "Uncertain_significance" & stars >= 1, .(chrom, pos, ref, alt, clinvar_id_2024 = clinvar_id, stars_2024 = stars)]
rm(cv24); invisible(gc())
cat("ClinVar 2024 SNV VUS (>=1 star):", nrow(vus), "\n")
cv26 <- read_cv("C:/cadd/annot/clinvar.vcf.gz")[, .(chrom, pos, ref, alt, clnsig_2026 = clnsig, stars_2026 = stars, clinvar_id_2026 = clinvar_id)]
cv26 <- unique(cv26, by = k)
vus <- merge(vus, cv26, by = k, all.x = TRUE)
rm(cv26); invisible(gc())
vus[, outcome := fcase(clnsig_2026 %chin% plp, "P/LP", clnsig_2026 %chin% blb, "B/LB", default = "Not resolved")]

## MANE missense in disease genes; AF, CADD, AVI ---------------------------------------------------
af <- readRDS(file.path(cdir, "gnomad_af_mane_coding.rds"))[, .(chrom, pos, ref, alt, AF_max)]
v <- rbindlist(lapply(list.files(cdir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f) {
  d <- readRDS(f)[, .(chrom, pos, ref, alt, gene_name, transcript_id, hgvsp, consequence)]
  d <- d[vus, on = k, nomatch = NULL]
  if (!nrow(d)) return(NULL)
  bad <- unique(d[consequence %chin% c("stop_gained", "stop_lost", "start_lost", "incomplete_terminal_codon_variant"), ..k])
  d <- d[consequence == "missense"][!bad, on = k]
  x <- d[, .(gene = paste(sort(unique(gene_name)), collapse = "&"), transcript = paste(unique(transcript_id), collapse = "&"),
             hgvsp = paste(unique(hgvsp), collapse = "&"), dgene = any(gene_name %chin% dg)), by = k][dgene == TRUE][, dgene := NULL]
  if (!nrow(x)) return(NULL)
  ch <- x$chrom[1]
  s1 <- readRDS(file.path(cdir, sprintf("cadd_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, CADD_PHRED = PHRED)]
  s2 <- readRDS(file.path(cdir, sprintf("avi_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, AVI_PHRED = round(AVI_PHRED, 2))]
  merge(merge(x, s1, by = k), s2, by = k)
}))
v <- merge(v, vus, by = k)
v[af, on = k, AF_max := i.AF_max][is.na(AF_max), AF_max := 0]
rm(af); invisible(gc())
cat("MANE missense VUS (2024) in disease genes:", nrow(v), "| AF < 0.01:", v[AF_max < 0.01, .N], "\n")
print(v[, .N, by = .(outcome, rare = AF_max < 0.01)][order(rare, outcome)])

## Comparator scores (rules of scripts 20 and 23) --------------------------------------------------
keys <- v[, ..k]
am <- fread(file.path(cd, "AlphaMissense_hg38.tsv.gz"), skip = "#CHROM", select = c(1, 2, 3, 4, 9),
            col.names = c("chrom", "pos", "ref", "alt", "AM"), colClasses = list(character = 1))
am[, chrom := sub("^chr", "", chrom)]; am <- am[keys, on = k, nomatch = NULL][, .(AM = max(AM)), by = k]; invisible(gc())
rv <- fread(cmd = sprintf('C:\\Windows\\System32\\tar.exe -xOf "%s" revel_with_transcript_ids',
                          normalizePath(file.path(cd, "revel-v1.3_all_chromosomes.zip"), winslash = "\\")),
            select = c("chr", "grch38_pos", "ref", "alt", "REVEL"), colClasses = list(character = c("chr", "grch38_pos")))
rv <- rv[grch38_pos != "."][, .(chrom = chr, pos = as.integer(grch38_pos), ref, alt, REVEL)]
rv <- rv[keys, on = k, nomatch = NULL][, .(REVEL = max(REVEL)), by = k]; invisible(gc())
chain <- import.chain(file.path(cd, "bayesdel", "hg38ToHg19.over.chain"))
lo <- liftOver(GRanges(paste0("chr", keys$chrom), IRanges(keys$pos, width = 1), strand = "+"), chain)
ok <- which(lengths(lo) == 1); l1 <- unlist(lo[ok]); comp <- c(A = "T", C = "G", G = "C", T = "A")
bd <- data.table(keys[ok], h19_chr = as.character(seqnames(l1)), h19_pos = start(l1), h19_strand = as.character(strand(l1)))
bd[, `:=`(ref19 = ifelse(h19_strand == "-", comp[ref], ref), alt19 = ifelse(h19_strand == "-", comp[alt], alt))]
tb <- TwoBitFile(file.path(cd, "bayesdel", "hg19.2bit")); bd <- bd[h19_chr %in% seqnames(seqinfo(tb))]
bd[, base19 := as.character(getSeq(tb, GRanges(h19_chr, IRanges(h19_pos, width = 1))))]
bd <- bd[base19 == ref19][, BayesDel := NA_real_]
for (b in c("A", "C", "G", "T")) {
  idx <- which(bd$alt19 == b); if (!length(idx)) next
  q <- GRanges(bd$h19_chr[idx], IRanges(bd$h19_pos[idx], width = 1))
  x <- import(BigWigFile(file.path(cd, "bayesdel", sprintf("NoAFBayesDelMut%s.bw", b))), which = reduce(q))
  ov <- findOverlaps(q, x); bd$BayesDel[idx[queryHits(ov)]] <- score(x)[subjectHits(ov)]
}
v <- Reduce(function(a, b) merge(a, b, by = k, all.x = TRUE), list(v, am, rv, bd[, c(k, "BayesDel"), with = FALSE]))
v[, `:=`(REVEL = round(REVEL, 3), AM = round(AM, 4), BayesDel = round(BayesDel, 2))]
saveRDS(v, file.path(out, "vus2024_scores.rds"))

## Thresholds -------------------------------------------------------------------------------------
thr <- list(
  r2024 = fread("C:/cadd/calibration/results/v17_avi/release_series/thresholds_all.tsv")[set == "r2024",
            .(tool, level, PP3 = PP3_discounted, BP4 = BP4_discounted)],
  main2026 = fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main",
            .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)],
  inter2026 = fread("C:/cadd/calibration/results/v17_avi/comparators/summary/thresholds_all_runs.tsv")[,
            .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)])
fwrite(rbindlist(thr, idcol = "thresholds"), file.path(out, "thresholds_used.tsv"), sep = "\t")
pts <- c(PP3_VeryStrong = 8, PP3_Strong = 4, PP3_Moderate = 2, PP3_Supporting = 1, Indeterminate = 0,
         BP4_Supporting = -1, BP4_Moderate = -2, BP4_Strong = -4, BP4_VeryStrong = -8)
lrt <- c(PP3_VeryStrong = 1124, PP3_Strong = 33.53, PP3_Moderate = 5.79, PP3_Supporting = 2.406, Indeterminate = NA,
         BP4_Supporting = 1 / 2.406, BP4_Moderate = 1 / 5.79, BP4_Strong = 1 / 33.53, BP4_VeryStrong = 1 / 1124)
assign_level <- function(s, m) {
  lv <- rep("Indeterminate", length(s))
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { t <- m[level == L, PP3]; if (length(t) && !is.na(t)) lv[!is.na(s) & s >= t] <- paste0("PP3_", L) }
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { t <- m[level == L, BP4]; if (length(t) && !is.na(t)) lv[!is.na(s) & s <= t] <- paste0("BP4_", L) }
  lv[is.na(s)] <- NA; lv
}
cols <- c(CADD = "CADD_PHRED", AVI = "AVI_PHRED", REVEL = "REVEL", AlphaMissense = "AM", BayesDel = "BayesDel")
B <- 2000L
direction <- list(); lr <- list()
for (ts in names(thr)) {
  tl <- if (ts == "main2026") c("CADD", "AVI") else names(cols)
  for (tn in tl) {
    lv <- assign_level(v[[cols[[tn]]]], thr[[ts]][tool == tn])
    v[, (paste0("ev_", ts, "_", tn)) := lv]
    for (sub in c("AF<0.01", "all")) {
      keep <- !is.na(lv) & (sub == "all" | v$AF_max < 0.01)
      d <- data.table(outcome = v$outcome[keep], ev = lv[keep])
      d[, p := pts[ev]]
      direction[[paste(ts, tn, sub)]] <- d[, .(N = .N, PP3_any = round(100 * mean(p > 0), 1), BP4_any = round(100 * mean(p < 0), 1),
                                              indeterminate = round(100 * mean(p == 0), 1), mean_points = round(mean(p), 2)),
                                           by = outcome][, `:=`(thresholds = ts, tool = tn, subset = sub)]
      # interval LRs: resolved P/LP vs resolved B/LB
      pl <- d[outcome == "P/LP", ev]; bl <- d[outcome == "B/LB", ev]
      lev <- names(pts)
      fr <- function(x) tabulate(factor(x, levels = lev), length(lev)) / length(x)
      obs <- fr(pl) / fr(bl)
      bs <- replicate(B, fr(sample(pl, replace = TRUE)) / fr(sample(bl, replace = TRUE)))
      lr[[paste(ts, tn, sub)]] <- data.table(thresholds = ts, tool = tn, subset = sub, evidence = lev,
        n_PLP = tabulate(factor(pl, levels = lev), length(lev)), n_BLB = tabulate(factor(bl, levels = lev), length(lev)),
        LR = signif(obs, 3), LR_lo = signif(apply(bs, 1, quantile, 0.025, na.rm = TRUE), 3),
        LR_hi = signif(apply(bs, 1, quantile, 0.975, na.rm = TRUE), 3), LR_required = signif(lrt[lev], 4))
    }
  }
}
dir_t <- rbindlist(direction); lr_t <- rbindlist(lr)
lr_t[, meets := fifelse(is.na(LR_required) | n_PLP + n_BLB == 0, NA,
                        fifelse(LR_required >= 1, LR_lo >= LR_required, LR_hi <= LR_required))]
fwrite(dir_t, file.path(out, "table_direction.tsv"), sep = "\t")
fwrite(lr_t, file.path(out, "table_interval_LR.tsv"), sep = "\t")
cat("\n== Direction of evidence by 2026 outcome (r2024 thresholds, AF < 0.01) ==\n")
print(dcast(dir_t[thresholds == "r2024" & subset == "AF<0.01"], tool ~ outcome, value.var = c("PP3_any", "BP4_any")))
cat("\n== Interval LRs (r2024 thresholds, AF < 0.01) ==\n")
print(lr_t[thresholds == "r2024" & subset == "AF<0.01" & n_PLP + n_BLB > 0])

## Concordance: evidence pointing in the final direction ------------------------------------------
conc <- rbindlist(lapply(names(cols), function(tn) {
  e <- v[[paste0("ev_r2024_", tn)]]; r <- v$AF_max < 0.01 & !is.na(e) & v$outcome != "Not resolved"
  d <- data.table(outcome = v$outcome[r], p = pts[e[r]])
  d[, .(tool = tn, N = .N, concordant = sum((outcome == "P/LP" & p > 0) | (outcome == "B/LB" & p < 0)),
        discordant = sum((outcome == "P/LP" & p < 0) | (outcome == "B/LB" & p > 0)), none = sum(p == 0))]
}))
conc[, `:=`(pct_concordant = round(100 * concordant / N, 1), pct_discordant = round(100 * discordant / N, 1))]
cat("\n== Concordance with the 2026 classification (r2024 thresholds, AF < 0.01, resolved VUS) ==\n"); print(conc)
fwrite(conc, file.path(out, "table_concordance.tsv"), sep = "\t")

## Illustrative examples (rule fixed in A15) -------------------------------------------------------
v[, AVI_points := pts[ev_r2024_AVI]]
ex <- v[outcome != "Not resolved" & stars_2026 >= 3 & !is.na(AVI_points)]
ex <- ex[order(outcome, -abs(AVI_points), chrom, pos)][, head(.SD, 3), by = outcome]
ex <- ex[, .(outcome, gene, transcript, hgvsp, chrom, pos, ref, alt, clinvar_id_2026, clnsig_2026, stars_2026, AF_max,
             CADD_PHRED, AVI_PHRED, REVEL, AM, BayesDel, ev_r2024_CADD, ev_r2024_AVI, ev_r2024_REVEL, ev_r2024_AlphaMissense, ev_r2024_BayesDel)]
cat("\n== Illustrative examples ==\n"); print(ex)
fwrite(ex, file.path(out, "table_examples.tsv"), sep = "\t")
fwrite(v[, .N, by = .(outcome, rare = AF_max < 0.01)], file.path(out, "set_sizes.tsv"), sep = "\t")
saveRDS(v, file.path(out, "vus2024_evidence.rds"))
