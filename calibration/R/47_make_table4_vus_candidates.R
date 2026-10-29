# Table 4 of manuscript v7 (amendment A16): genes with the most ClinVar missense VUS reviewed by an expert
# panel or practice guideline that receive concordant calibrated computational evidence, with the number of
# variants at each level and up to three examples at PP3 strong. Source: script 46 outputs (no new computation).
# Output: C:/cadd/figures/Table4_vus_candidates.tsv
suppressPackageStartupMessages(library(data.table))
vc <- "C:/cadd/calibration/results/v17_avi/vus_candidates"
e <- fread(file.path(vc, "vus_concordant_expert_panel.tsv"))
v <- readRDS(file.path(vc, "vus2026_evidence.rds"))[stars >= 3 & AF_max < 0.01]
tot <- v[, .(expert_panel_VUS = .N), by = gene]
lv <- c("4" = "PP3_strong", "2" = "PP3_moderate", "1" = "PP3_supporting", "-1" = "BP4_supporting", "-2" = "BP4_moderate")
e[, level := factor(lv[as.character(comp_points)], levels = lv)]
t4 <- dcast(e, gene ~ level, fun.aggregate = length, value.var = "pos", drop = FALSE)
t4 <- merge(tot, t4, by = "gene")
t4[, concordant := rowSums(.SD), .SDcols = lv]
ex <- e[comp_points == 4][order(chrom, pos), .(examples_PP3_strong = paste(head(hgvsp, 3), collapse = ", ")), by = gene]
t4 <- merge(t4, ex, by = "gene", all.x = TRUE)[concordant >= 10][order(-concordant, gene)]
setcolorder(t4, c("gene", "expert_panel_VUS", "concordant", lv, "examples_PP3_strong"))
fwrite(t4, "C:/cadd/figures/Table4_vus_candidates.tsv", sep = "\t", na = "")
print(t4)
cat("All expert-panel VUS:", nrow(v), "| concordant:", nrow(e), "| genes:", uniqueN(e$gene), "| genes in Table 4:", nrow(t4),
    "| variants in Table 4:", sum(t4$concordant), "\n")
