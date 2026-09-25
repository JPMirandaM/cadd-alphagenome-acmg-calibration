# Inspects the Aim 1 gene annotation files: columns, key fields and category counts.
# Usage: Rscript A8_inspect_aim1_gene_annotations.R
suppressPackageStartupMessages(library(data.table))
d <- "C:/cadd/annot/genes"
f <- function(p) list.files(d, p, full.names = TRUE)[1]

cat("== gnomAD v4.1 constraint ==\n")
con <- fread(f("constraint_metrics"))
cat("rows:", nrow(con), "\n"); print(grep("gene|transcript|mane|canonical|lof.oe_ci|lof.pLI|mis.z", names(con), value = TRUE))
if ("mane_select" %in% names(con)) cat("MANE Select rows:", con[mane_select == TRUE, .N], "\n")

cat("\n== GenCC submissions ==\n")
gc <- fread(f("gencc_submissions"))
cat("rows:", nrow(gc), "\n"); print(names(gc))
print(gc[, .N, by = moi_title][order(-N)][1:10])
print(gc[, .N, by = classification_title][order(-N)])

cat("\n== ClinGen gene-disease validity ==\n")
cv <- fread(f("clingen_gene_validity"), skip = "GENE SYMBOL", fill = TRUE)
cat("rows:", nrow(cv), "\n"); print(names(cv))
if ("MOI" %in% names(cv)) print(cv[, .N, by = MOI][order(-N)])

cat("\n== ClinGen dosage sensitivity ==\n")
ds <- fread(f("clingen_dosage"), skip = "#Gene Symbol", sep = "\t", fill = TRUE)
cat("rows:", nrow(ds), "\n"); print(names(ds))
hi <- grep("Haploinsufficiency Score", names(ds), value = TRUE)[1]
if (!is.na(hi)) print(ds[, .N, by = hi][order(-N)])

cat("\n== Gene2Phenotype ==\n")
g2 <- fread(f("gene2phenotype"))
cat("rows:", nrow(g2), "\n"); print(names(g2))
for (k in grep("allelic|mechanism|confidence|panel", names(g2), ignore.case = TRUE, value = TRUE))
  { cat("--", k, "\n"); print(g2[, .N, by = k][order(-N)][1:12]) }
