# Amendment A12: population frequencies of PP3-level MANE Select missense SNVs in the Mexico City
# Prospective Study (MCPS; RGC-MCPS release 20230130; ~140,000 adults; local-ancestry-specific
# frequencies: AFR, EUR, IMX = Indigenous Mexican). Descriptive; no recalibration.
# Data: RGC-MCPS allele-frequency VCFs (licence: research use, citation/acknowledgement required,
# no redistribution). Only aggregate counts are written to the repository.
suppressMessages({ library(data.table); library(Rsamtools); library(GenomicRanges); library(parallel) })

cat_dir <- "C:/cadd/data/catalogue"
mcps_dir <- "C:/cadd/data/mcps"
out_dir <- "C:/cadd/calibration/results/v17_avi/mcps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
chroms <- c(1:22, "X")
gnomad_groups <- paste0("AF_", c("AFR", "AMR", "ASJ", "EAS", "FIN", "MID", "NFE", "SAS", "AMI"))
fields <- c("AC_RAW", "AN_RAW", "AF_RAW", "AF_AFR", "AF_EUR", "AF_IMX", "AN_IMX", "AF_RAW_UNREL")

get_info <- function(info, key) {
  v <- regmatches(info, regexpr(paste0("(^|;)", key, "=[^;]*"), info))
  out <- rep(NA_character_, length(info))
  hit <- regexpr(paste0("(^|;)", key, "=[^;]*"), info) > 0
  out[hit] <- sub(paste0("^;?", key, "="), "", v)
  suppressWarnings(as.numeric(out))
}

query_mcps <- function(chr, qpos) {
  # positions -> reduced ranges (gap 1 kb) -> tabix -> exact position filter
  gr <- reduce(GRanges(chr, IRanges(qpos, qpos)), min.gapwidth = 1000L)
  tf <- TabixFile(file.path(mcps_dir, paste0("chr", chr, ".freq.vcf.gz")))
  res <- list(); step <- 2000L
  for (i in seq(1L, length(gr), by = step)) {
    lines <- unlist(scanTabix(tf, param = gr[i:min(i + step - 1L, length(gr))]), use.names = FALSE)
    if (!length(lines)) next
    x <- tstrsplit(lines, "\t", fixed = TRUE, keep = c(1L, 2L, 4L, 5L, 8L))
    d <- data.table(chrom = x[[1]], pos = as.integer(x[[2]]), ref = x[[3]], alt = x[[4]], info = x[[5]])
    d <- d[nchar(ref) == 1L & nchar(alt) == 1L & pos %in% qpos]
    res[[length(res) + 1L]] <- d
  }
  d <- unique(rbindlist(res), by = c("pos", "ref", "alt"))
  if (!nrow(d)) return(NULL)
  d[, source := sub(";.*", "", sub("^SOURCE=", "", info))]
  for (f in fields) set(d, j = paste0("mcps_", f), value = get_info(d$info, f))
  d[, info := NULL]
  d[]
}

run_chr <- function(chr) {
  suppressMessages({ library(data.table); library(Rsamtools); library(GenomicRanges) })
  cat_f <- file.path(cat_dir, paste0("catalogue_chr", chr, ".tsv.gz"))
  k <- fread(cat_f, select = c("chrom", "pos", "ref", "alt", "CADD_PHRED", "AVI_PHRED", "PP3_CADD",
                               "PP3_AVI", "calibrated_context", "mane_gene", "hgvsp", "clnsig",
                               "clinvar_stars", "AF_max", "grpmax", "gnomAD_status", gnomad_groups))
  k <- k[calibrated_context == TRUE]
  # gnomAD group maximum (max over genetic ancestry groups of max(exome, genome) AF); AF_max in the
  # catalogue is the overall AF (max of exome and genome), used for the plan-v1 rarity rule.
  k[, gnomad_grp_AF := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = gnomad_groups]
  k[is.na(gnomad_grp_AF), gnomad_grp_AF := 0]
  k[, (setdiff(gnomad_groups, "AF_AMR")) := NULL]
  k[, chrom := as.character(chrom)]
  m <- query_mcps(chr, unique(k$pos))
  if (is.null(m)) { k[, mcps_observed := FALSE]; return(k) }
  m[, chrom := as.character(chrom)]
  k <- merge(k, m, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE)
  k[, mcps_observed := !is.na(mcps_AN_RAW) & mcps_AC_RAW > 0]
  k[]
}

cl <- makeCluster(4L)
clusterExport(cl, c("cat_dir", "mcps_dir", "fields", "get_info", "query_mcps", "gnomad_groups"))
t0 <- Sys.time()
res <- parLapply(cl, chroms, run_chr)
stopCluster(cl)
k <- rbindlist(res, fill = TRUE)
cat("MANE missense PP3-level catalogue SNVs:", nrow(k), "| observed (AC>0) in MCPS:", sum(k$mcps_observed),
    "| elapsed min:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
saveRDS(k, file.path(mcps_dir, "catalogue_missense_mcps.rds"))

# Reference-allele sanity: MCPS rows present at a catalogue position with a different REF
cat("Source of matched variants:\n"); print(k[mcps_observed == TRUE, .N, by = source][order(-N)])

lvl_order <- c("Supporting", "Moderate", "Strong")
summ <- rbindlist(lapply(c("CADD", "AVI"), function(p) {
  col <- paste0("PP3_", p)
  rbindlist(lapply(lvl_order, function(l) {
    s <- k[get(col) == l]
    if (!nrow(s)) return(NULL)
    gmax <- s$gnomad_grp_AF
    data.table(
      predictor = p, level = l, n_snvs = nrow(s),
      n_observed_mcps = sum(s$mcps_observed),
      n_mcps_af_gt_0.01 = sum(s$mcps_AF_RAW > 0.01, na.rm = TRUE),
      n_mcps_af_gt_0.05 = sum(s$mcps_AF_RAW > 0.05, na.rm = TRUE),
      n_mcps_gt_0.01_gnomad_grpmax_le_0.01 = sum(s$mcps_AF_RAW > 0.01 & gmax <= 0.01, na.rm = TRUE),
      n_mcps_gt_0.05_gnomad_grpmax_le_0.05 = sum(s$mcps_AF_RAW > 0.05 & gmax <= 0.05, na.rm = TRUE),
      n_imx_af_gt_0.01 = sum(s$mcps_AF_IMX > 0.01, na.rm = TRUE),
      n_imx_gt_0.01_gnomad_grpmax_le_0.01 = sum(s$mcps_AF_IMX > 0.01 & gmax <= 0.01, na.rm = TRUE),
      n_unrel_af_gt_0.01 = sum(s$mcps_AF_RAW_UNREL > 0.01, na.rm = TRUE),
      n_observed_mcps_not_gnomad = sum(s$mcps_observed & s$gnomAD_status == "not observed")
    )
  }))
}))
fwrite(summ, file.path(out_dir, "pp3_levels_mcps_summary.tsv"), sep = "\t")
print(summ)

# ClinVar P/LP (>= 1 star) among catalogue SNVs: quality check
plp <- k[grepl("^Pathogenic|^Likely_pathogenic", clnsig) & clinvar_stars >= 1]
cv <- data.table(set = "ClinVar P/LP >=1 star (catalogue)", n = nrow(plp),
                 n_observed_mcps = sum(plp$mcps_observed),
                 n_mcps_af_gt_0.01 = sum(plp$mcps_AF_RAW > 0.01, na.rm = TRUE),
                 n_mcps_af_gt_0.05 = sum(plp$mcps_AF_RAW > 0.05, na.rm = TRUE),
                 n_imx_af_gt_0.01 = sum(plp$mcps_AF_IMX > 0.01, na.rm = TRUE))
fwrite(cv, file.path(out_dir, "clinvar_plp_mcps.tsv"), sep = "\t"); print(cv)

# Concordance MCPS AF vs gnomAD v4.1 AMR (variants with AF > 0 in both)
both <- k[mcps_AF_RAW > 0 & AF_AMR > 0]
cc <- data.table(n_both = nrow(both),
                 spearman_mcps_vs_gnomad_amr = if (nrow(both) > 2) cor(both$mcps_AF_RAW, both$AF_AMR, method = "spearman") else NA_real_,
                 spearman_imx_vs_gnomad_amr = if (nrow(both) > 2) cor(both$mcps_AF_IMX, both$AF_AMR, method = "spearman", use = "complete.obs") else NA_real_)
fwrite(cc, file.path(out_dir, "concordance_mcps_gnomad_amr.tsv"), sep = "\t"); print(cc)

# MCPS-only frequency conflicts: PP3 (either predictor) with MCPS AF > 0.01 and gnomAD group max <= 0.01
conf <- k[(PP3_CADD != "none" | PP3_AVI != "none") & mcps_AF_RAW > 0.01 & gnomad_grp_AF <= 0.01,
          .(chrom, pos, ref, alt, gene = mane_gene, hgvsp, CADD_PHRED, AVI_PHRED, PP3_CADD, PP3_AVI,
            mcps_AF_RAW = signif(mcps_AF_RAW, 3), mcps_AF_IMX = signif(mcps_AF_IMX, 3),
            mcps_AF_EUR = signif(mcps_AF_EUR, 3), mcps_AF_AFR = signif(mcps_AF_AFR, 3),
            gnomAD_AF_AMR = signif(AF_AMR, 3), gnomAD_AF_group_max = signif(gnomad_grp_AF, 3), gnomAD_AF_overall = signif(AF_max, 3), grpmax, clnsig)][order(-mcps_AF_RAW)]
fwrite(conf, file.path(out_dir, "mcps_only_frequency_conflicts.tsv"), sep = "\t")
cat("MCPS-only frequency conflicts (PP3 by either predictor, MCPS AF > 0.01, gnomAD group max <= 0.01):", nrow(conf), "\n")
print(head(conf, 20))

# Per-gene counts of MCPS-common (AF > 0.01) PP3-level SNVs
gs <- k[(PP3_CADD != "none" | PP3_AVI != "none") & mcps_AF_RAW > 0.01, .N, by = .(gene = mane_gene)][order(-N)]
fwrite(gs, file.path(out_dir, "genes_mcps_common_pp3.tsv"), sep = "\t")
cat("Done.\n")
