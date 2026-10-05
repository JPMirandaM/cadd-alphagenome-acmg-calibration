# Table 2 (manuscript v5): calibrated evidence for ClinVar VUS of January 2024 resolved by September
# 2026 (amendment A15), with thresholds calibrated on the January 2024 release only (r2024, five-predictor
# intersection) and gnomAD AF < 0.01. Built only from the outputs of script 42.
# Output: figures/Table2_vus_resolved.tsv and figures/Table2_vus_resolved.md (Markdown pipe table).
# Usage: Rscript 45_make_table2_vus_resolved.R
suppressPackageStartupMessages(library(data.table))
d <- "C:/cadd/calibration/results/v17_avi/vus_resolved"; fig <- "C:/cadd/figures"
dd <- file.path(fig, "data", "Table2_vus_resolved"); dir.create(dd, recursive = TRUE, showWarnings = FALSE)
for (f in c("table_interval_LR.tsv", "table_direction.tsv", "table_concordance.tsv")) file.copy(file.path(d, f), dd, overwrite = TRUE)
tools <- c(CADD = "CADD v1.7", AVI = "AVI", REVEL = "REVEL", BayesDel = "BayesDel", AlphaMissense = "AlphaMissense")
lv <- c("PP3_Strong", "PP3_Moderate", "PP3_Supporting", "Indeterminate", "BP4_Supporting", "BP4_Moderate", "BP4_Strong")
f2 <- function(x) ifelse(x < 0.1, formatC(x, format = "f", digits = 3), formatC(x, format = "f", digits = 2))
lr <- fread(file.path(dd, "table_interval_LR.tsv"))[thresholds == "r2024" & subset == "AF<0.01" & evidence %in% lv]
lr[, cell := fifelse(n_PLP + n_BLB == 0, "—",
               sprintf("%s [%s–%s] (%s / %s)%s", f2(LR), f2(LR_lo), f2(LR_hi), prettyNum(n_PLP, big.mark = ","), prettyNum(n_BLB, big.mark = ","),
                       fifelse(!is.na(meets) & !meets, "*", "")))]
t <- dcast(lr, evidence ~ tool, value.var = "cell")[match(lv, evidence)]
dr <- fread(file.path(dd, "table_direction.tsv"))[thresholds == "r2024" & subset == "AF<0.01"]
cc <- fread(file.path(dd, "table_concordance.tsv"))
extra <- rbind(
  dcast(dr[outcome == "P/LP", .(evidence = "Resolved P/LP with PP3, % (mean points)", tool, cell = sprintf("%.1f (%+.2f)", PP3_any, mean_points))], evidence ~ tool, value.var = "cell"),
  dcast(dr[outcome == "B/LB", .(evidence = "Resolved B/LB with BP4, % (mean points)", tool, cell = sprintf("%.1f (%+.2f)", BP4_any, mean_points))], evidence ~ tool, value.var = "cell"),
  dcast(cc[, .(evidence = "Evidence in the final direction / opposite, %", tool, cell = sprintf("%.1f / %.1f", pct_concordant, pct_discordant))], evidence ~ tool, value.var = "cell"))
t <- rbind(t, extra)[, c("evidence", names(tools)), with = FALSE]
t[, evidence := sub("_", " ", evidence)]
setnames(t, c("Interval", unname(tools)))
fwrite(t, file.path(fig, "Table2_vus_resolved.tsv"), sep = "\t")
md <- c(paste0("| ", paste(names(t), collapse = " | "), " |"), paste0("|", strrep("---|", ncol(t))),
        apply(t, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, file.path(fig, "Table2_vus_resolved.md"), useBytes = TRUE)
print(t)
