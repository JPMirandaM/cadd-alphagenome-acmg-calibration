# Amendment A16 (Aim 7c): current ClinVar missense VUS whose classification could change with calibrated
# computational evidence, and the evidence they receive under thresholds transferred between CADD versions.
# VUS: the 1,208,461 ClinVar 2026-09-23 VUS of script 23 (>= 1 star, MANE Select missense, disease gene).
# Evidence: CADD v1.7 and AVI at the 2026 main thresholds; REVEL, BayesDel, AlphaMissense at the published
# ClinGen thresholds (as in script 25b). Current practice: CADD v1.7 read with published CADD v1.6 thresholds.
# Concordant evidence: REVEL, BayesDel, AlphaMissense and AVI in the same direction; points = smallest |points|
# among the three missense predictors. No variant is removed; AF_max is reported.
# Output: C:/cadd/calibration/results/v17_avi/vus_candidates/
suppressPackageStartupMessages(library(data.table))
cdir <- "C:/cadd/data/coding"; out <- "C:/cadd/calibration/results/v17_avi/vus_candidates"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")

## VUS scores (script 23) and ClinVar 2026 record ---------------------------------------------------
cache <- file.path(out, "vus2026_annotated.rds")
if (file.exists(cache)) v <- readRDS(cache) else {
v <- readRDS("C:/cadd/calibration/results/v17_avi/vus/vus_scores.rds")
getk <- function(x, kk) { r <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); r[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; r }
cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 3, 4, 5, 8),
            col.names = c("chrom", "pos", "clinvar_id", "ref", "alt", "info"), colClasses = list(character = c(1, 3)))
cv <- cv[nchar(ref) == 1 & nchar(alt) == 1]
cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"))][, info := NULL]
cv <- cv[clnsig == "Uncertain_significance"]
cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                    grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                    grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
cv <- unique(cv[stars >= 1, .(chrom, pos, ref, alt, clinvar_id, stars)], by = k)
v <- merge(v, cv, by = k, all.x = TRUE)
rm(cv); invisible(gc())
hg <- rbindlist(lapply(list.files(cdir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f) {
  d <- readRDS(f)[, .(chrom, pos, ref, alt, transcript_id, hgvsp, consequence)][v[, ..k], on = k, nomatch = NULL]
  d[consequence == "missense", .(transcript = paste(unique(transcript_id), collapse = "&"),
                                 hgvsp = paste(unique(hgvsp), collapse = "&")), by = k]
}))
v <- merge(v, hg, by = k, all.x = TRUE)
saveRDS(v, cache) }
cat("VUS:", nrow(v), "| with ClinVar ID:", v[!is.na(clinvar_id), .N], "| AF < 0.01:", v[AF_max < 0.01, .N], "\n")

## Thresholds ---------------------------------------------------------------------------------------
main <- fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main",
          .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)]
pub <- data.table(tool = c(rep("REVEL", 7), rep("BayesDel", 5), rep("AlphaMissense", 7)),
  crit = c("PP3", "PP3", "PP3", "BP4", "BP4", "BP4", "BP4", "PP3", "PP3", "PP3", "BP4", "BP4", "PP3", "PP3", "PP3", "BP4", "BP4", "BP4", "PP3"),
  level = c("Supporting", "Moderate", "Strong", "Supporting", "Moderate", "Strong", "VeryStrong",
            "Supporting", "Moderate", "Strong", "Supporting", "Moderate",
            "Supporting", "Moderate", "Strong", "Supporting", "Moderate", "Strong", "VeryStrong"),
  value = c(0.644, 0.773, 0.932, 0.290, 0.183, 0.016, 0.003, 0.13, 0.27, 0.50, -0.18, -0.36,
            0.792, 0.906, 0.972, 0.169, 0.099, 0.070, 0.990))   # same values as script 25b (Table S5)
pub <- dcast(pub, tool + level ~ crit, value.var = "value")
v16 <- data.table(tool = "CADD16", level = c("Supporting", "Moderate", "Strong"), PP3 = c(25.3, 28.1, NA), BP4 = c(22.7, 17.3, 0.15))
thr <- rbind(main[tool %chin% c("CADD", "AVI"), .(tool, level, PP3, BP4)], pub[, .(tool, level, PP3, BP4)], v16, fill = TRUE)
fwrite(thr, file.path(out, "thresholds_used.tsv"), sep = "\t")
pts <- c(PP3_VeryStrong = 8, PP3_Strong = 4, PP3_Moderate = 2, PP3_Supporting = 1, Indeterminate = 0,
         BP4_Supporting = -1, BP4_Moderate = -2, BP4_Strong = -4, BP4_VeryStrong = -8)
assign_level <- function(s, m) {
  lv <- rep("Indeterminate", length(s))
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { t <- m[level == L, PP3]; if (length(t) && !is.na(t)) lv[!is.na(s) & s >= t] <- paste0("PP3_", L) }
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { t <- m[level == L, BP4]; if (length(t) && !is.na(t)) lv[!is.na(s) & s <= t] <- paste0("BP4_", L) }
  lv[is.na(s)] <- NA; lv
}
cols <- c(CADD = "CADD_PHRED", AVI = "AVI_PHRED", REVEL = "REVEL", BayesDel = "BayesDel", AlphaMissense = "AM", CADD16 = "CADD_PHRED")
for (tn in names(cols)) v[, (paste0("ev_", tn)) := assign_level(get(cols[[tn]]), thr[tool == tn])]

## Concordant evidence and point scenario ------------------------------------------------------------
P <- sapply(c("REVEL", "BayesDel", "AlphaMissense", "AVI"), function(tn) pts[v[[paste0("ev_", tn)]]])
miss <- P[, 1:3]
v[, concordant := complete.cases(P) & (apply(P > 0, 1, all) | apply(P < 0, 1, all))]
v[, comp_points := fifelse(concordant, sign(miss[, 1]) * apply(abs(miss), 1, min), 0)]
v[, scenario := fcase(comp_points > 0, sprintf("Likely pathogenic if other evidence >= +%d", 6 - comp_points),
                      comp_points < 0, sprintf("Likely benign if other evidence <= +%d", -1 - comp_points),
                      default = "No concordant computational evidence")]

## Transferred CADD v1.6 thresholds versus calibrated evidence ---------------------------------------
dirn <- function(e) sign(pts[e])
v[, `:=`(change_CADD = ev_CADD16 != ev_CADD, change_AVI = ev_CADD16 != ev_AVI,
         flip_CADD = dirn(ev_CADD16) * dirn(ev_CADD) == -1, flip_AVI = dirn(ev_CADD16) * dirn(ev_AVI) == -1)]

## Summaries -------------------------------------------------------------------------------------------
r <- v[AF_max < 0.01]
summ <- rbind(
  r[, .N, keyby = comp_points][, .(item = "concordant points", value = as.character(comp_points), N)],
  r[stars >= 3, .N, keyby = comp_points][, .(item = "concordant points (3-4 stars)", value = as.character(comp_points), N)],
  data.table(item = c("change level: CADD v1.6 transferred vs CADD v1.7 calibrated", "change level: CADD v1.6 transferred vs AVI calibrated",
                      "flip direction: CADD v1.6 transferred vs CADD v1.7 calibrated", "flip direction: CADD v1.6 transferred vs AVI calibrated"),
             value = "TRUE", N = c(r[change_CADD == TRUE, .N], r[change_AVI == TRUE, .N], r[flip_CADD == TRUE, .N], r[flip_AVI == TRUE, .N])),
  data.table(item = "rare VUS", value = "all", N = nrow(r)))
summ[, pct := round(100 * N / nrow(r), 2)]
fwrite(summ, file.path(out, "summary.tsv"), sep = "\t")
print(summ)
xt <- r[, .N, by = .(CADD16 = ev_CADD16, CADD17 = ev_CADD)][order(-N)]
fwrite(xt, file.path(out, "crosstab_cadd16_vs_cadd17.tsv"), sep = "\t")
xa <- r[, .N, by = .(CADD16 = ev_CADD16, AVI = ev_AVI)][order(-N)]
fwrite(xa, file.path(out, "crosstab_cadd16_vs_avi.tsv"), sep = "\t")
print(r[, .N, by = .(dir16 = dirn(ev_CADD16), dirAVI = dirn(ev_AVI))][order(dir16, dirAVI)])

## Lists ---------------------------------------------------------------------------------------------
keep <- c("gene", "transcript", "hgvsp", "chrom", "pos", "ref", "alt", "clinvar_id", "stars", "AF_max",
          "CADD_PHRED", "AVI_PHRED", "REVEL", "BayesDel", "AM", "ev_CADD16", "ev_CADD", "ev_AVI", "ev_REVEL", "ev_BayesDel",
          "ev_AlphaMissense", "comp_points", "scenario", "change_AVI", "flip_AVI")
lst <- v[concordant == TRUE][order(-abs(comp_points), chrom, pos), ..keep]
fwrite(lst, file.path(out, "vus_concordant_evidence.tsv.gz"), sep = "\t")
ep <- lst[stars >= 3]
fwrite(ep, file.path(out, "vus_concordant_expert_panel.tsv"), sep = "\t")
cat("\nConcordant VUS:", nrow(lst), "| rare:", lst[AF_max < 0.01, .N], "| expert panel/guideline:", nrow(ep), "\n")
print(ep[, .N, by = .(gene, comp_points)][order(-abs(comp_points), -N)])
saveRDS(v, file.path(out, "vus2026_evidence.rds"))
