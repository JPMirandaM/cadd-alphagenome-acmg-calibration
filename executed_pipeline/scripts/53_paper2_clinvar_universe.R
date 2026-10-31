# Paper 2 (analysis-plan-paper2-v1): ClinVar SNVs whose AVI attributions are retrieved from the Atlas API.
# All SNVs with >= 1 star classified P/LP, B/LB or VUS in the 2026-09-23 or the 2024-01-27 release, in a
# disease gene of either release (>= 1 P/LP record with >= 1 star; GENEINFO), any variant class. This is a
# superset of the reference sets, the temporal hold-out, the resolved VUS and the current VUS; region class,
# MANE scope and rarity are applied in the analysis scripts. Only variant keys are written; labels are joined
# after the strata are assigned (plan, "Mechanism strata").
# Disease genes are defined from P/LP records of any variant type (as in scripts 30 and 54). A first run that
# defined them from SNV records only (6,337 genes) was already being queried when the difference was found, so
# its list is kept as variants_clinvar_batch1.tsv and the SNVs it missed are written to variants_clinvar_batch2.tsv.
# Output: C:/cadd/calibration/results/paper2/universe/variants_clinvar.tsv (full list), variants_clinvar_batch2.tsv
suppressPackageStartupMessages(library(data.table))
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
read_cv <- function(f) {
  cv <- fread(f, sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5, 8), col.names = c("chrom", "pos", "ref", "alt", "info"),
              colClasses = list(character = 1))
  cv[, snv := grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt)]
  cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
  cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                      grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                      grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
  cv[stars >= 1 & !is.na(geneinfo)]
}
genes_of <- function(g) unique(sub(":.*", "", unlist(strsplit(g, "|", fixed = TRUE))))
cv <- rbind(read_cv("C:/cadd/annot/clinvar.vcf.gz")[, release := "2026"],
            read_cv("C:/cadd/annot/clinvar_archive/clinvar_20240127.vcf.gz")[, release := "2024"])
dg <- cv[clnsig %chin% plp, genes_of(geneinfo)]
cv <- cv[snv == TRUE & clnsig %chin% c(plp, blb, "Uncertain_significance")]
s <- strsplit(cv$geneinfo, "|", fixed = TRUE); g <- sub(":.*", "", unlist(s)); r <- rep(seq_along(s), lengths(s))
cv[, in_dg := seq_len(.N) %in% r[g %chin% dg]]; rm(s, g, r)
cv <- cv[in_dg == TRUE & chrom %chin% c(1:22, "X", "Y", "MT")]
cat("Disease genes:", length(dg), "\n"); print(cv[, .N, by = .(release, group = fcase(clnsig %chin% plp, "P/LP", clnsig %chin% blb, "B/LB", default = "VUS"))])
u <- unique(cv[, .(chrom, pos, ref, alt)])[, chrom := fifelse(chrom == "MT", "M", chrom)]
setorder(u, chrom, pos)
u[, variant := sprintf("chr%s:%d:%s>%s", chrom, pos, ref, alt)]
fwrite(u[, .(variant)], "C:/cadd/calibration/results/paper2/universe/variants_clinvar.tsv", sep = "\t")
cat("Unique SNVs to query:", nrow(u), "\n")
b1 <- fread("C:/cadd/calibration/results/paper2/universe/variants_clinvar_batch1.tsv")
stopifnot(all(b1$variant %chin% u$variant))
b2 <- u[!b1, on = "variant", .(variant)]
fwrite(b2, "C:/cadd/calibration/results/paper2/universe/variants_clinvar_batch2.tsv", sep = "\t")
cat("Already queried (batch 1):", nrow(b1), "| added (batch 2):", nrow(b2), "\n")
