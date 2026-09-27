# Amendment A12(b): overall allele frequencies of PP3-level MANE Select missense SNVs in TOPMed
# BRAVO (Freeze 10, GRCh38; public per-chromosome VCFs; overall AC/AF only, no ancestry).
# The downloaded files are plain gzip; they are recompressed with bgzip and tabix-indexed (the
# originals are kept). Only FILTER = PASS sites are used. Descriptive; no recalibration.
# BRAVO overlaps with gnomAD v4.1 genomes (TOPMed samples): not an independent replication.
suppressMessages({ library(data.table); library(Rsamtools); library(GenomicRanges); library(parallel) })
bdir <- "C:/cadd/data/bravo"; cat_dir <- "C:/cadd/data/catalogue"
out_dir <- "C:/cadd/calibration/results/v17_avi/bravo"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
chroms <- c(1:22, "X")
gnomad_groups <- paste0("AF_", c("AFR", "AMR", "ASJ", "EAS", "FIN", "MID", "NFE", "SAS", "AMI"))

run_chr <- function(chr) {
  suppressMessages({ library(data.table); library(Rsamtools); library(GenomicRanges) })
  src <- file.path(bdir, sprintf("chr%s.bravo.pub.vcf.gz", chr))
  bgz <- file.path(bdir, sprintf("chr%s.bravo.pub.bgz.vcf.gz", chr))
  if (!file.exists(paste0(bgz, ".tbi"))) { bgzip(src, bgz, overwrite = TRUE); indexTabix(bgz, format = "vcf") }
  k <- fread(file.path(cat_dir, paste0("catalogue_chr", chr, ".tsv.gz")),
             select = c("chrom", "pos", "ref", "alt", "PP3_CADD", "PP3_AVI", "calibrated_context", "mane_gene",
                        "clnsig", "clinvar_stars", "AF_max", gnomad_groups))
  k <- k[calibrated_context == TRUE]
  k[, gnomad_grp_AF := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = gnomad_groups][is.na(gnomad_grp_AF), gnomad_grp_AF := 0]
  k[, (gnomad_groups) := NULL]
  gr <- reduce(GRanges(paste0("chr", chr), IRanges(unique(k$pos), unique(k$pos))), min.gapwidth = 1000L)
  tf <- TabixFile(bgz); res <- list()
  for (i in seq(1L, length(gr), by = 2000L)) {
    l <- unlist(scanTabix(tf, param = gr[i:min(i + 1999L, length(gr))]), use.names = FALSE)
    if (!length(l)) next
    x <- tstrsplit(l, "\t", fixed = TRUE, keep = c(2L, 4L, 5L, 7L, 8L))
    res[[length(res) + 1L]] <- data.table(pos = as.integer(x[[1]]), ref = x[[2]], alt = x[[3]], filter = x[[4]], info = x[[5]])
  }
  b <- unique(rbindlist(res), by = c("pos", "ref", "alt"))
  b <- b[nchar(ref) == 1L & nchar(alt) == 1L]
  b[, `:=`(bravo_AC = as.integer(sub(".*(^|;)AC=([0-9]+).*", "\\2", info)),
           bravo_AF = as.numeric(sub(".*(^|;)AF=([^;]+).*", "\\2", info)),
           bravo_hom = as.integer(sub(".*(^|;)Hom=([0-9]+).*", "\\2", info)))][, info := NULL]
  k <- merge(k, b, by = c("pos", "ref", "alt"), all.x = TRUE)
  k[, bravo_pass := !is.na(filter) & filter == "PASS"]
  k[]
}

cl <- makeCluster(3L); clusterExport(cl, c("bdir", "cat_dir", "gnomad_groups"))
t0 <- Sys.time(); k <- rbindlist(parLapplyLB(cl, chroms, run_chr), fill = TRUE); stopCluster(cl)
saveRDS(k, file.path(bdir, "catalogue_missense_bravo.rds"))
cat("MANE missense PP3-level catalogue SNVs:", nrow(k), "| at a BRAVO site:", sum(!is.na(k$filter)),
    "| PASS:", sum(k$bravo_pass), "| elapsed min:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
print(k[!is.na(filter), .N, by = .(pass = bravo_pass)])

summ <- rbindlist(lapply(c("CADD", "AVI"), function(p) rbindlist(lapply(c("Supporting", "Moderate", "Strong"), function(l) {
  s <- k[get(paste0("PP3_", p)) == l]; if (!nrow(s)) return(NULL); ps <- s[bravo_pass == TRUE]
  data.table(predictor = p, level = l, n_snvs = nrow(s), n_observed_bravo_pass = ps[bravo_AC > 0, .N],
             n_bravo_af_gt_0.01 = ps[bravo_AF > 0.01, .N], n_bravo_af_gt_0.05 = ps[bravo_AF > 0.05, .N],
             n_bravo_gt_0.01_gnomad_grpmax_le_0.01 = ps[bravo_AF > 0.01 & gnomad_grp_AF <= 0.01, .N],
             n_bravo_homozygotes_gt0 = ps[bravo_hom > 0, .N])
}))))
fwrite(summ, file.path(out_dir, "pp3_levels_bravo_summary.tsv"), sep = "\t"); print(summ)
plp <- k[grepl("^Pathogenic|^Likely_pathogenic", clnsig) & clinvar_stars >= 1 & bravo_pass == TRUE]
cv <- data.table(set = "ClinVar P/LP >=1 star (catalogue), BRAVO PASS", n_observed = plp[bravo_AC > 0, .N],
                 n_af_gt_0.01 = plp[bravo_AF > 0.01, .N], n_af_gt_0.05 = plp[bravo_AF > 0.05, .N])
fwrite(cv, file.path(out_dir, "clinvar_plp_bravo.tsv"), sep = "\t"); print(cv)
both <- k[bravo_pass == TRUE & bravo_AF > 0 & AF_max > 0]
cc <- data.table(n_both = nrow(both), spearman_bravo_vs_gnomad_overall = cor(both$bravo_AF, both$AF_max, method = "spearman"))
fwrite(cc, file.path(out_dir, "concordance_bravo_gnomad.tsv"), sep = "\t"); print(cc)
conf <- k[(PP3_CADD != "none" | PP3_AVI != "none") & bravo_pass == TRUE & bravo_AF > 0.01 & gnomad_grp_AF <= 0.01,
          .(chrom, pos, ref, alt, gene = mane_gene, PP3_CADD, PP3_AVI, bravo_AF = signif(bravo_AF, 3),
            gnomAD_group_max = signif(gnomad_grp_AF, 3), gnomAD_overall = signif(AF_max, 3), clnsig)][order(-bravo_AF)]
fwrite(conf, file.path(out_dir, "bravo_only_frequency_conflicts.tsv"), sep = "\t")
cat("BRAVO-only frequency conflicts (PP3 either predictor, BRAVO PASS AF > 0.01, gnomAD group max <= 0.01):", nrow(conf), "\n")
print(head(conf, 15))
