# Main Tables 1 and 2 of the manuscript (TSV + Word), built only from source tables copied into
# figures/data/Table1_thresholds/ and figures/data/Table2_mave/.
# Table 1: calibrated PP3/BP4 thresholds of CADD v1.7 and AVI (main analysis; median [range] over
#          five seeds), with the evidence levels assigned to the reference variants and the published
#          CADD v1.6 thresholds; if the A14 results exist, CADD v1.6 recalibrated on the same 2026 data.
# Table 2: likelihood ratios of the calibrated evidence intervals on MAVE functional classes (A13),
#          five predictors, with pooled AUROC.
# Usage: Rscript 39_make_tables_1_2.R <figures_dir>
suppressPackageStartupMessages({ library(data.table); library(officer) })
fig <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(fig)) fig <- "C:/cadd/figures"
d1 <- file.path(fig, "data", "Table1_thresholds"); d2 <- file.path(fig, "data", "Table2_mave")
fmt <- function(x, d) formatC(x, format = "f", digits = d)
dig <- function(tool) if (tool == "CADD") 2 else 2

## Table 1 -----------------------------------------------------------------------------------------
runs <- fread(file.path(d1, "thresholds_all_runs.tsv"))[analysis == "main"]
ev <- fread(file.path(d1, "table_evidence_levels_reference.tsv"))
lv <- c("VeryStrong", "Strong", "Moderate", "Supporting")
req <- c(VeryStrong = 1124, Strong = 33.53, Moderate = 5.790, Supporting = 2.406)
cell <- function(v, na_frac) {
  v <- v[!is.na(v)]
  if (!length(v) || median(na_frac) > 0.05) return(if (median(na_frac) < 1) sprintf("not reached (reached in %.0f%% of replicates)", 100 * (1 - median(na_frac))) else "not reached")
  m <- median(v); r <- range(v); nd <- if (all(abs(v - round(v, 1)) < 1e-9)) 1 else 2
  if (r[1] == r[2]) fmt(m, nd) else sprintf("%s [%s–%s]", fmt(m, nd), fmt(r[1], nd), fmt(r[2], nd))
}
pub16 <- c(PP3_VeryStrong = "not reached", PP3_Strong = "not reached", PP3_Moderate = "≥ 28.1", PP3_Supporting = "≥ 25.3",
           BP4_Supporting = "≤ 22.7", BP4_Moderate = "≤ 17.3", BP4_Strong = "≤ 0.15", BP4_VeryStrong = "not reached")
rows <- list()
for (crit in c("PP3", "BP4")) {
  lvls <- if (crit == "PP3") lv else rev(lv)
  for (L in lvls) {
    key <- paste0(crit, "_", L)
    get_cell <- function(tn) { x <- runs[tool == tn & level == L]
      v <- if (crit == "PP3") x$PP3_discounted else x$BP4_discounted
      na <- if (crit == "PP3") x$PP3_boot_NA_frac else x$BP4_boot_NA_frac
      s <- cell(v, na); if (grepl("^not", s)) s else paste(if (crit == "PP3") "≥" else "≤", s) }
    evp <- function(tn, set) { e <- ev[tool == tn & evidence == key]; if (!nrow(e)) "—" else fmt(e[[set]], 1) }
    rows[[key]] <- data.table(Evidence = paste(crit, sub("VeryStrong", "Very strong", L)),
      `Required LR` = if (crit == "PP3") as.character(req[L]) else paste0("1/", req[L]),
      `CADD v1.7 PHRED` = get_cell("CADD"), `CADD P/LP %` = evp("CADD", "P/LP"), `CADD B/LB %` = evp("CADD", "B/LB"),
      `AVI Phred` = get_cell("AVI"), `AVI P/LP %` = evp("AVI", "P/LP"), `AVI B/LB %` = evp("AVI", "B/LB"),
      `CADD v1.6 PHRED (Pejaver et al. 2022)` = pub16[key])
  }
}
t1 <- rbindlist(rows)
f16 <- file.path(d1, "cadd16_2026_thresholds_median.tsv")
if (file.exists(f16)) {
  m16 <- fread(f16)[tool == "CADD16"]
  t1[, `CADD v1.6 PHRED (2026 data)` := vapply(seq_len(.N), function(i) {
    crit <- sub(" .*", "", Evidence[i]); L <- sub("Very strong", "VeryStrong", sub("^\\S+ ", "", Evidence[i]))
    v <- m16[level == L][[crit]]; if (!length(v) || is.na(v)) "not reached" else paste(if (crit == "PP3") "≥" else "≤", fmt(v, if (abs(v - round(v, 1)) < 1e-9) 1 else 2)) }, "")]
}
fwrite(t1, file.path(fig, "Table1_thresholds.tsv"), sep = "\t")

## Table 2 -----------------------------------------------------------------------------------------
lr <- fread(file.path(d2, "interval_LR.tsv"))[set == "five-tool intersection"]
au <- fread(file.path(d2, "auroc_pooled.tsv"))[set == "five-tool intersection"]
tools <- c("CADD", "AVI", "REVEL", "BayesDel", "AM"); tl <- c("CADD v1.7", "AVI", "REVEL", "BayesDel", "AlphaMissense")
ord <- c("PP3_Strong", "PP3_Moderate", "PP3_Supporting", "Indeterminate", "BP4_Supporting", "BP4_Moderate", "BP4_Strong")
t2 <- rbindlist(lapply(ord, function(o) {
  r <- data.table(Interval = sub("_", " ", o))
  for (j in seq_along(tools)) { x <- lr[tool == tools[j] & interval == o]
    r[[tl[j]]] <- if (!nrow(x) || x$n_abnormal + x$n_normal == 0) "—" else
      sprintf("%.2f [%.2f–%.2f] (%s / %s)", x$LR, x$LR_lo, x$LR_hi, format(x$n_abnormal, big.mark = ","), format(x$n_normal, big.mark = ",")) }
  r }))
aur <- data.table(Interval = "AUROC [95% CI]")
for (j in seq_along(tools)) { x <- au[tool == tools[j]]; aur[[tl[j]]] <- sprintf("%.3f [%.3f–%.3f]", x$AUROC, x$lo, x$hi) }
t2 <- rbind(t2, aur)
fwrite(t2, file.path(fig, "Table2_mave.tsv"), sep = "\t")

## Word -------------------------------------------------------------------------------------------
add_table <- function(doc, title, dt, notes) {
  doc <- body_add_par(doc, title, style = "heading 2")
  doc <- body_add_table(doc, as.data.frame(dt), style = "table_template")
  for (n in notes) doc <- body_add_par(doc, n, style = "Normal")
  doc
}
doc <- read_docx()
doc <- add_table(doc, "Table 1. Calibrated PP3/BP4 thresholds for CADD v1.7 and AVI (rare missense variants, ClinVar 2026)", t1, c(
  "Thresholds are the one-sided 95% bootstrap bounds (10,000 replicates) of the ClinGen local posterior probability method (prior 0.0441), reported as the median [range] over five random seeds; a single value means identical results in all seeds. PP3 applies when the score is at or above the threshold and BP4 when it is at or below it.",
  "Reference set: 59,979 P/LP and 77,128 B/LB rare MANE Select missense ClinVar variants (4,826 genes; ≥ 1 star) and 1,955,995 rare gnomAD v4.1 missense variants. P/LP % and B/LB %: percentage of reference variants assigned to each evidence level with the median thresholds (in-sample).",
  "LR, likelihood ratio (odds of pathogenicity) required for each strength under the Bayesian adaptation of the ACMG/AMP guidelines. The AVI PP3 strong level was reached in 59% of bootstrap replicates, fewer than the 95% required.",
  "CADD v1.6 thresholds from Pejaver et al. (2022; ClinVar 2019 with HGMD- and training-variant exclusions); these were reproduced from the original data with our implementation (Table S1). The last column gives CADD v1.6 recalibrated on the same 2026 reference sets as CADD v1.7 (median of five seeds)."))
doc <- body_add_break(doc)
doc <- add_table(doc, "Table 2. Calibrated evidence intervals on MAVE functional classes (93,154 MANE missense SNVs, 39 genes)", t2, c(
  "Cells: likelihood ratio P(interval | functionally abnormal) / P(interval | functionally normal) [95% bootstrap interval, 2,000 replicates], and, in parentheses, numbers of abnormal / normal variants in the interval. Intervals are defined by each predictor's thresholds calibrated on ClinVar (CADD v1.7 and AVI: main analysis; REVEL, BayesDel and AlphaMissense: five-predictor intersection analysis). —: level not reached by the predictor.",
  "Functional classes (18,936 abnormal, 74,218 normal) from 89 MaveDB score sets with investigator-provided or other normal/abnormal calibrations; SNVs classified discordantly by different score sets were excluded. Functional abnormality is not equivalent to pathogenicity, so the ratios are not directly comparable with the ClinGen likelihood ratio targets; they test the direction and ordering of the evidence intervals on data independent of ClinVar."))
doc <- body_end_section_landscape(doc)
print(doc, target = file.path(fig, "Tables_1_2.docx"))
print(t1); print(t2)
cat("Written:", file.path(fig, c("Table1_thresholds.tsv", "Table2_mave.tsv", "Tables_1_2.docx")), sep = "\n")
