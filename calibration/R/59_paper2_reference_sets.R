# Paper 2 (plan v1 "Reference sets"; amendment A-P2-1 "First step and stopping rule"): reference sets of the
# 2026-09-23 (calibration) and 2024-01-27 (hold-out calibration) ClinVar releases for every region class, with the
# mechanism stratum of each SNV from the Atlas attributions, and their sizes. No score is summarised here: the log
# reports only counts of P/LP and B/LB variants per mechanism stratum and per region class.
# Rules (as in scripts 13, 18 and 30): SNVs with >= 1 star and an exact P/LP or B/LB label; gene (MANE Select class
# gene) among the disease genes of the same release (>= 1 P/LP record of any variant type with >= 1 star);
# rare (gnomAD v4.1 exomes/genomes group-maximum AF < 0.01, VEP 116 cache; absent = 0); intergenic excluded.
# Test set (temporal hold-out, as in script 18): 2026 set variants absent from the 2024 set; label changes are
# excluded and counted. SNVs without Atlas attributions (chrM) are excluded and counted.
# Output: C:/cadd/calibration/results/paper2/refsets/ (refsets.rds, sizes_*.tsv)
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
source("C:/cadd/calibration/R/paper2_regions.R")
out <- "C:/cadd/calibration/results/paper2/refsets"; tmpd <- file.path(out, "tmp"); dir.create(tmpd, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
read_cv <- function(f) {
  cv <- fread(f, sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5, 8), col.names = c("chrom", "pos", "ref", "alt", "info"),
              colClasses = list(character = 1))
  cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
  cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                      grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                      grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
  dg <- cv[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo), unique(sub(":.*", "", unlist(strsplit(geneinfo, "|", fixed = TRUE))))]
  list(dg = dg, ref = cv[grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt) & stars >= 1 & clnsig %chin% c(plp, blb),
                        .(chrom, pos, ref, alt, y = as.integer(clnsig %chin% plp))])
}
rel <- list(r2026 = read_cv("C:/cadd/annot/clinvar.vcf.gz"),
            r2024 = read_cv("C:/cadd/annot/clinvar_archive/clinvar_20240127.vcf.gz"))
for (r in names(rel)) cat(sprintf("%s: %d disease genes, %d P/LP or B/LB SNVs (>= 1 star)\n", r, length(rel[[r]]$dg), nrow(rel[[r]]$ref)))

## Region classes of the union of SNVs -------------------------------------------------------------
u <- unique(rbind(rel$r2026$ref[, ..k], rel$r2024$ref[, ..k]))
coding_pos <- unique(readRDS(file.path(CDIR, "mane_cds_positions.rds"))[, .(chrom, pos)])
uc <- u[coding_pos, on = .(chrom, pos), nomatch = NULL]; un <- u[!coding_pos, on = .(chrom, pos)]
edge <- exon_edge_positions()
cc <- rbindlist(lapply(split(uc, uc$chrom), function(x) classify_coding_chr(x, edge)))
m <- mane_structure()
nc <- merge(un, classify_noncoding(un, m), by = c("chrom", "pos"))[!class %chin% c("exonic_other", "intergenic")]
reg <- rbind(cc, nc[, .(chrom, pos, ref, alt, gene, class)])
cat("SNVs with a region class:", nrow(reg), "of", nrow(u), "\n")

## gnomAD AF: coding table (13a) and VEP cache for noncoding SNVs, cached per chromosome ---------------
af_c <- readRDS(file.path(CDIR, "gnomad_af_mane_coding.rds"))[, .(chrom, pos, ref, alt, AF_max)]
vep <- "C:/cadd/vep/homo_sapiens/116_GRCh38"
vcols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", readLines(file.path(vep, "info.txt")), value = TRUE)), ",")[[1]]
af_chr <- function(x) {   # as in scripts 30 and 54
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
xn <- unique(nc[, ..k]); chs <- unique(xn$chrom); todo <- chs[!file.exists(file.path(tmpd, sprintf("af_%s.rds", chs)))]
if (length(todo)) {
  cl <- parallel::makeCluster(6L)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) }))
  parallel::clusterExport(cl, c("vep", "vcols", "af_chr", "tmpd"))
  parallel::parLapplyLB(cl, lapply(todo, function(ch) xn[chrom == ch]), function(x) { saveRDS(af_chr(x), file.path(tmpd, sprintf("af_%s.rds", x$chrom[1]))); NULL })
  parallel::stopCluster(cl); cat("AF lookups:", length(todo), "chromosomes", format(Sys.time(), "%H:%M"), "\n")
}
af_n <- rbindlist(lapply(chs, function(ch) readRDS(file.path(tmpd, sprintf("af_%s.rds", ch)))))
af <- unique(rbind(af_c, af_n[, .(chrom, pos, ref, alt, AF_max)])[, .(AF_max = max(AF_max)), by = k])
reg <- merge(reg, af, by = k, all.x = TRUE)[is.na(AF_max), AF_max := 0]

## Mechanism strata from the Atlas attributions (plan v1 rule) ---------------------------------------
fi <- rbindlist(lapply(list.files(c("C:/cadd/calibration/results/paper2/atlas_clinvar", "C:/cadd/calibration/results/paper2/atlas_clinvar_batch2"),
                                  "^chunk_.*\\.tsv\\.gz$", full.names = TRUE), fread, sep = "\t", colClasses = list(character = "error")))
fi <- fi[error == "" | is.na(error)]
pos_sum <- function(cols) Reduce(`+`, lapply(cols, function(cn) pmax(fi[[cn]], 0)))
reg_f <- grep("^fi_MAX_ABS_", names(fi), value = TRUE)
fi[, `:=`(g_spl = pos_sum("fi_MERGED_SPLICING"), g_pro = pos_sum(c("fi_ALPHAMISSENSE", "fi_PROTEIN_TERMINATION", "fi_START_LOST", "fi_STOP_LOST")),
          g_reg = pos_sum(reg_f), g_con = pos_sum(c("fi_CACTUS_241_WAY", "fi_PHASTCONS_470_WAY")))]
fi[, mech := g_spl + g_pro + g_reg]
fi[, stratum := fifelse(mech >= 0.05,
                        fifelse(g_spl >= 0.5 * mech & g_spl >= g_pro & g_spl >= g_reg, "splicing",
                        fifelse(g_pro >= 0.5 * mech & g_pro >= g_reg, "protein",
                        fifelse(g_reg >= 0.5 * mech, "regulatory", "mixed"))),
                        fifelse(g_con >= 0.05, "conservation_only", "no_predicted_effect"))]
v <- tstrsplit(sub("^chr", "", fi$variant), "[:>]"); fi[, `:=`(chrom = fifelse(v[[1]] == "M", "MT", v[[1]]), pos = as.integer(v[[2]]), ref = v[[3]], alt = v[[4]])]
fi <- fi[, c(k, "avi_raw", "avi_cdf", grep("^fi_", names(fi), value = TRUE), "stratum"), with = FALSE]
reg <- merge(reg, fi, by = k, all.x = TRUE)
cat("SNVs with a region class but no Atlas attributions (before the gene and rarity filters):", reg[is.na(stratum), .N],
    "| of which chrM:", reg[is.na(stratum) & chrom == "MT", .N], "\n")
# The query list (script 53) took disease genes from GENEINFO; here the gene is the MANE Select class gene, so
# SNVs whose GENEINFO gene is not a disease gene can enter a reference set through an overlapping MANE gene.
# They have no attributions; they are counted after the reference-set filters and excluded.

## Reference sets ------------------------------------------------------------------------------------
mk <- function(r) {
  s <- merge(rel[[r]]$ref, reg, by = k)[gene %chin% rel[[r]]$dg & AF_max < 0.01]
  cat(sprintf("%s reference SNVs without Atlas attributions (excluded): %d (P/LP %d, B/LB %d; chrM %d)\n", r,
              s[is.na(stratum), .N], s[is.na(stratum) & y == 1, .N], s[is.na(stratum) & y == 0, .N], s[is.na(stratum) & chrom == "MT", .N]))
  print(head(s[is.na(stratum), .N, by = class][order(-N)], 8))
  s[!is.na(stratum)]
}
s26 <- mk("r2026"); s24 <- mk("r2024")
chg <- merge(s26[, c(k, "y"), with = FALSE], s24[, c(k, "y"), with = FALSE], by = k)[y.x != y.y]
test <- s26[!s24, on = k]
cat(sprintf("Reference SNVs: 2026 %d, 2024 %d; test (in 2026, absent from 2024) %d; label changes between releases (excluded from test) %d\n",
            nrow(s26), nrow(s24), nrow(test), nrow(chg)))
refsets <- rbind(s26[, set := "cal_2026"], s24[, set := "cal_2024"], test[, set := "test_2026_new"])
saveRDS(refsets, file.path(out, "refsets.rds"))

## Sizes (counts only) ---------------------------------------------------------------------------------
lev_s <- c("splicing", "protein", "regulatory", "mixed", "conservation_only", "no_predicted_effect")
size_tab <- function(by) {
  z <- dcast(refsets[, .N, by = c("set", by, "y")][, lab := fifelse(y == 1, "P_LP", "B_LB")], as.formula(paste(by, "~ set + lab")), value.var = "N", fill = 0)
  z[, calibratable_2026 := cal_2026_P_LP >= 1000 & cal_2026_B_LB >= 1000][, calibratable_2024 := cal_2024_P_LP >= 1000 & cal_2024_B_LB >= 1000]
  z
}
ss <- size_tab("stratum")[order(match(stratum, lev_s))]; sc <- size_tab("class")[order(match(class, REGION_CLASSES))]
options(width = 250)
cat("\n== Sizes by mechanism stratum ==\n"); print(ss); fwrite(ss, file.path(out, "sizes_by_stratum.tsv"), sep = "\t")
cat("\n== Sizes by region class ==\n"); print(sc); fwrite(sc, file.path(out, "sizes_by_class.tsv"), sep = "\t")
x <- dcast(refsets[set == "cal_2026", .N, by = .(class, stratum, lab = fifelse(y == 1, "P", "B"))], class ~ stratum + lab, value.var = "N", fill = 0)
x <- x[order(match(class, REGION_CLASSES))]
cat("\n== 2026 calibration set: region class x mechanism stratum (P = P/LP, B = B/LB) ==\n"); print(x)
fwrite(x, file.path(out, "class_by_stratum_2026.tsv"), sep = "\t")
