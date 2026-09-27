# Correction to the Aim 4 per-gene summary (script 29): A9 specifies "common (group AF > 1%)",
# but script 29 counted common variants with the overall gnomAD AF (max of exome and genome AF).
# This script recomputes the per-gene counts with the gnomAD group maximum (max over genetic
# ancestry groups of max(exome, genome) AF) and keeps the overall-AF count for transparency.
suppressMessages(library(data.table))
out <- "C:/cadd/data/catalogue"
grp <- paste0("AF_", c("AFR", "AMR", "ASJ", "EAS", "FIN", "MID", "NFE", "SAS", "AMI"))
chroms <- c(as.character(1:22), "X", "Y")
res <- rbindlist(lapply(chroms, function(ch) {
  d <- fread(file.path(out, sprintf("catalogue_chr%s.tsv.gz", ch)),
             select = c("PP3_CADD", "PP3_AVI", "calibrated_context", "mane_gene", "AF_max", "grpmax", grp),
             colClasses = list(character = c("mane_gene", "grpmax")))
  d <- d[calibrated_context == TRUE & !grepl("&", mane_gene, fixed = TRUE)]
  d[, grp_AF := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = grp]
  d[, .(common_overall_gt_1pct = sum(AF_max > 0.01, na.rm = TRUE),
        common_group_gt_1pct = sum(grp_AF > 0.01, na.rm = TRUE),
        common_group_gt_5pct = sum(grp_AF > 0.05, na.rm = TRUE),
        common_group_gt_1pct_PP3_CADD = sum(grp_AF > 0.01 & PP3_CADD != "none", na.rm = TRUE),
        common_group_gt_1pct_PP3_AVI = sum(grp_AF > 0.01 & PP3_AVI != "none", na.rm = TRUE)), by = .(gene = mane_gene)]
}))
res <- res[, lapply(.SD, sum), by = gene]
gs <- fread(file.path(out, "gene_summary.tsv"))
gs[, common_gt_1pct := NULL]
gs <- merge(gs, res, by = "gene", all.x = TRUE)
fwrite(gs[order(-missense_catalogued)], file.path(out, "gene_summary.tsv"), sep = "\t")
cat("Genes:", nrow(gs), "\n")
cat("Totals: common overall > 1%:", sum(gs$common_overall_gt_1pct, na.rm = TRUE),
    "| common in some group > 1%:", sum(gs$common_group_gt_1pct, na.rm = TRUE),
    "| > 5%:", sum(gs$common_group_gt_5pct, na.rm = TRUE), "\n")
cat("\nTop genes by overall AF > 1% (as reported in R12):\n"); print(head(gs[order(-common_overall_gt_1pct), .(gene, common_overall_gt_1pct, common_group_gt_1pct)], 10))
cat("\nTop genes by group AF > 1% (corrected definition):\n"); print(head(gs[order(-common_group_gt_1pct), .(gene, missense_catalogued, common_group_gt_1pct, common_group_gt_5pct, common_overall_gt_1pct)], 15))
