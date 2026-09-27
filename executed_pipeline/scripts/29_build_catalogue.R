# Aim 4 (amendment A9): genome-wide catalogue of SNVs with CADD v1.7 PHRED >= 25.6 or AVI Phred >= 25.18,
# with both scores, PP3 levels, MANE consequence / GENCODE gene overlap, ClinVar and gnomAD v4.1.
# Processed per chromosome (4 parallel workers); output per chromosome plus a per-gene summary.
# Output: C:/cadd/data/catalogue/catalogue_chr<chrom>.tsv.gz, gene_summary.tsv, catalogue_summary.tsv
suppressPackageStartupMessages({ library(data.table); library(parallel) })
cdir <- "C:/cadd/data/coding"; out <- "C:/cadd/data/catalogue"; tmp <- file.path(out, "tmp")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
thr <- fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main", .(PP3 = median(PP3_discounted)), by = .(tool, level)]
TH <- list(CADD = setNames(thr[tool == "CADD", PP3], thr[tool == "CADD", level]), AVI = setNames(thr[tool == "AVI", PP3], thr[tool == "AVI", level]))
print(TH)

## Split inputs per chromosome (so that workers load only their chromosome) ----------------------
if (!file.exists(file.path(tmp, "split_done"))) {
  ca <- fread("C:/cadd/data/PHRED_ge25.6_CADDv1.7_GRCh38.tsv.gz", select = c("Chrom", "Pos", "Ref", "Alt", "PHRED"),
              col.names = c("chrom", "pos", "ref", "alt", "CADD_PHRED"), colClasses = list(character = "Chrom"))
  for (ch in unique(ca$chrom)) saveRDS(ca[chrom == ch], file.path(tmp, sprintf("cadd_%s.rds", ch)))
  rm(ca); invisible(gc())
  av <- fread("C:/cadd/avi/AVI_PHRED_ge25.18_GRCh38.tsv.gz", select = c("Chrom", "Pos", "Ref", "Alt", "AVI_PHRED"),
              col.names = c("chrom", "pos", "ref", "alt", "AVI_PHRED"), colClasses = list(character = "Chrom"))
  for (ch in unique(av$chrom)) saveRDS(av[chrom == ch], file.path(tmp, sprintf("avi_%s.rds", ch)))
  rm(av); invisible(gc())
  cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 3, 4, 5, 8),
              col.names = c("chrom", "pos", "clinvar_id", "ref", "alt", "info"), colClasses = list(character = c(1, 3)))
  cv <- cv[nchar(ref) == 1 & nchar(alt) == 1]
  getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
  cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"))][, info := NULL]
  cv[, clinvar_stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                              grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                              grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)][, revstat := NULL]
  for (ch in unique(cv$chrom)) saveRDS(cv[chrom == ch], file.path(tmp, sprintf("clinvar_%s.rds", ch)))
  rm(cv); invisible(gc())
  gtf <- fread("C:/cadd/annot/gencode.v50.basic.annotation.gtf.gz", sep = "\t", header = FALSE, skip = "chr1\t",
               col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))[type == "gene"]
  gtf[, `:=`(chrom = sub("^chr", "", chr), gene_name = sub('.*gene_name "([^"]+)".*', "\\1", attr),
             gene_type = sub('.*gene_type "([^"]+)".*', "\\1", attr))]
  gtf[chrom == "M", chrom := "MT"]
  saveRDS(gtf[, .(chrom, start, end, gene_name, gene_type)], file.path(tmp, "genes_gencode.rds"))
  af <- readRDS(file.path(cdir, "gnomad_af_mane_coding.rds"))
  for (ch in unique(af$chrom)) saveRDS(af[chrom == ch], file.path(tmp, sprintf("afcoding_%s.rds", ch)))
  rm(af, gtf); invisible(gc())
  file.create(file.path(tmp, "split_done"))
}

## Per-chromosome worker -----------------------------------------------------------------------------
build_chrom <- function(ch) {
  suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
  of <- file.path(out, sprintf("catalogue_chr%s.tsv.gz", ch))
  if (file.exists(of)) return(sprintf("%s: exists", ch))
  rd <- function(p) { f <- file.path(tmp, sprintf(p, ch)); if (file.exists(f)) readRDS(f) else NULL }
  ca <- rd("cadd_%s.rds"); av <- rd("avi_%s.rds")
  keys <- unique(rbindlist(list(if (!is.null(ca)) ca[, ..k], if (!is.null(av)) av[, ..k])))
  if (!nrow(keys)) return(sprintf("%s: empty", ch))
  d <- merge(keys, ca, by = k, all.x = TRUE); d <- merge(d, av, by = k, all.x = TRUE)
  # Missing score of the other predictor: tabix lookup over merged position regions
  lookup <- function(file, seqname, need, col, round2 = FALSE) {
    if (!nrow(need)) return(NULL)
    p <- sort(unique(need$pos)); grp <- cumsum(c(1L, diff(p) > 1L))
    rg <- data.table(s = p, g = grp)[, .(start = min(s), end = max(s)), by = g]
    tf <- open(TabixFile(file)); on.exit(close(tf))
    res <- rbindlist(lapply(split(rg, ceiling(seq_len(nrow(rg)) / 5000)), function(b) {
      l <- unlist(scanTabix(tf, param = GenomicRanges::GRanges(seqname, IRanges::IRanges(b$start, b$end))), use.names = FALSE)
      if (!length(l)) return(NULL)
      fread(text = l, sep = "\t", header = FALSE, select = c(2, 3, 4, 6), col.names = c("pos", "ref", "alt", "score"), nThread = 1)
    }))
    if (round2) res[, score := round(score, 2)]
    res[, chrom := ch]; setnames(res, "score", col); res[need, on = k, nomatch = NULL]
  }
  x <- lookup("C:/cadd/data/whole_genome_SNVs.tsv.gz", ch, d[is.na(CADD_PHRED), ..k], "CADD_PHRED")
  if (!is.null(x)) d[x, on = k, CADD_PHRED := i.CADD_PHRED]
  x <- lookup("C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz", paste0("chr", ch), d[is.na(AVI_PHRED), ..k], "AVI_PHRED", TRUE)
  if (!is.null(x)) d[x, on = k, AVI_PHRED := i.AVI_PHRED]
  d[, AVI_PHRED := round(AVI_PHRED, 2)]
  lvl <- function(s, t) fifelse(is.na(s), NA_character_, fifelse(s >= t["Moderate"], "Moderate", fifelse(s >= t["Supporting"], "Supporting", "none")))
  d[, `:=`(PP3_CADD = lvl(CADD_PHRED, TH$CADD), PP3_AVI = lvl(AVI_PHRED, TH$AVI))]
  # MANE consequence (coding) and GENCODE gene overlap (all)
  cf <- file.path(cdir, sprintf("coding_snv_%s.rds", ch))
  if (file.exists(cf)) {
    cs <- readRDS(cf)[, .(chrom, pos, ref, alt, gene_name, transcript_id, consequence, hgvsp)][d[, ..k], on = k, nomatch = NULL]
    cs <- cs[, .(mane_gene = paste(unique(gene_name), collapse = "&"), mane_transcript = paste(unique(transcript_id), collapse = "&"),
                 mane_consequence = paste(unique(consequence), collapse = "&"), hgvsp = paste(unique(hgvsp), collapse = "&")), by = k]
    d <- merge(d, cs, by = k, all.x = TRUE)
  } else d[, `:=`(mane_gene = NA_character_, mane_transcript = NA_character_, mane_consequence = NA_character_, hgvsp = NA_character_)]
  gn <- readRDS(file.path(tmp, "genes_gencode.rds"))[chrom == ch]
  setkey(gn, start, end)
  q <- unique(d[, .(pos)])[, `:=`(start = pos, end = pos)]; setkey(q, start, end)
  ov <- foverlaps(q, gn, type = "within", nomatch = NULL)
  ov <- ov[, .(genes = paste(sort(unique(gene_name)), collapse = "&"),
               gene_types = paste(sort(unique(gene_type)), collapse = "&")), by = pos]
  d <- merge(d, ov, by = "pos", all.x = TRUE)
  d[, region := fcase(!is.na(mane_consequence), "MANE coding", !is.na(genes), "gene body (non-MANE-coding)", default = "intergenic")]
  d[, calibrated_context := !is.na(mane_consequence) & mane_consequence == "missense"]
  # ClinVar
  cv <- rd("clinvar_%s.rds")
  if (!is.null(cv)) d <- merge(d, cv[, .(chrom, pos, ref, alt, clinvar_id, clnsig, clinvar_stars)], by = k, all.x = TRUE)
  # gnomAD: coding from the MANE table; all positions from a stream of the VEP cache chromosome file
  grp <- c("AFR", "AMR", "ASJ", "EAS", "FIN", "MID", "NFE", "SAS", "AMI")
  vep_dir <- "C:/cadd/vep/homo_sapiens/116_GRCh38"
  cols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", readLines(file.path(vep_dir, "info.txt")), value = TRUE)), ",")[[1]]
  afcols <- grep("^gnomAD[eg]", cols, value = TRUE)
  f <- file.path(vep_dir, ch, "all_vars.gz"); posset <- unique(d$pos); acc <- list()
  if (file.exists(f)) {
    con <- gzfile(f, "r")
    repeat {
      l <- readLines(con, n = 2e6); if (!length(l)) break
      st <- as.integer(sub("^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]*)\t.*", "\\1", l))
      l <- l[st %in% posset]; if (!length(l)) next
      x <- fread(text = l, sep = "\t", header = FALSE, colClasses = "character", nThread = 1,
                 select = which(cols %in% c("start", "allele_string", afcols)), col.names = cols[cols %in% c("start", "allele_string", afcols)])
      x <- x[grepl("^[ACGT](/[ACGT])+$", allele_string)]
      if (!nrow(x)) next
      lg <- rbindlist(lapply(afcols, function(kk) { y <- x[get(kk) != ".", .(pos = as.integer(start), allele_string, v = get(kk))]
        if (!nrow(y)) return(NULL); y <- y[, .(e = unlist(strsplit(v, ",", fixed = TRUE))), by = .(pos, allele_string)]
        y[, .(pos, ref = sub("/.*", "", allele_string), alt = sub(":.*", "", e), field = kk, af = as.numeric(sub(".*:", "", e)))] }))
      acc[[length(acc) + 1]] <- lg
    }
    close(con)
  }
  if (length(acc)) {
    lg <- rbindlist(acc)[nchar(alt) == 1]
    w <- dcast(lg, pos + ref + alt ~ field, value.var = "af", fun.aggregate = max, fill = NA)
    for (kk in afcols) if (kk %in% names(w)) set(w, which(is.infinite(w[[kk]])), kk, NA_real_)
    w[, AF_exomes := if ("gnomADe" %in% names(w)) gnomADe else NA_real_][, AF_genomes := if ("gnomADg" %in% names(w)) gnomADg else NA_real_]
    for (g in grp) { e <- paste0("gnomADe_", g); gg <- paste0("gnomADg_", g)
      set(w, j = paste0("AF_", g), value = pmax(if (e %in% names(w)) w[[e]] else NA_real_, if (gg %in% names(w)) w[[gg]] else NA_real_, na.rm = TRUE)) }
    w <- w[, c("pos", "ref", "alt", "AF_exomes", "AF_genomes", paste0("AF_", grp)), with = FALSE]
    d <- merge(d, w, by = c("pos", "ref", "alt"), all.x = TRUE)
  } else for (cc in c("AF_exomes", "AF_genomes", paste0("AF_", grp))) set(d, j = cc, value = NA_real_)
  d[, AF_max := pmax(AF_exomes, AF_genomes, na.rm = TRUE)]
  afc <- paste0("AF_", grp)
  d[, grpmax := { m <- as.matrix(.SD); m[is.na(m)] <- -1; ifelse(apply(m, 1, max) > 0, grp[max.col(m, ties.method = "first")], NA_character_) }, .SDcols = afc]
  d[, gnomAD_status := fifelse(!is.na(AF_max) & AF_max > 0, "observed", "not observed")]
  setcolorder(d, c("chrom", "pos", "ref", "alt", "CADD_PHRED", "AVI_PHRED", "PP3_CADD", "PP3_AVI", "calibrated_context", "region",
                   "mane_gene", "mane_transcript", "mane_consequence", "hgvsp", "genes", "gene_types"))
  setorder(d, pos, alt)
  fwrite(d, of, sep = "\t")
  sprintf("%s: %d SNVs", ch, nrow(d))
}

all_chroms <- c(as.character(1:22), "X", "Y")
a <- commandArgs(trailingOnly = TRUE)
chroms <- if (length(a) >= 1) strsplit(a[1], ",")[[1]] else all_chroms       # e.g. "21" for a test run
cl <- makeCluster(min(3, length(chroms)))    # 3 workers: large chromosomes have ~10x more SNVs than chr21
clusterExport(cl, c("cdir", "out", "tmp", "k", "TH"))
t0 <- Sys.time()
msg <- parLapplyLB(cl, chroms, build_chrom)
stopCluster(cl)
print(unlist(msg)); cat(sprintf("catalogue built in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## Summaries (only when every chromosome is done) ----------------------------------------------------
if (!all(file.exists(file.path(out, sprintf("catalogue_chr%s.tsv.gz", all_chroms))))) { cat("summaries skipped: not all chromosomes built\n"); quit(save = "no") }
chroms <- all_chroms
cat_all <- rbindlist(lapply(chroms, function(ch) fread(file.path(out, sprintf("catalogue_chr%s.tsv.gz", ch)),
  select = c("chrom", "PP3_CADD", "PP3_AVI", "calibrated_context", "region", "mane_gene", "clnsig", "gnomAD_status", "AF_max"),
  colClasses = list(character = c("chrom", "mane_gene", "clnsig")))))
summ <- cat_all[, .N, by = .(region, calibrated_context, PP3_CADD, PP3_AVI)][order(region, -N)]
print(summ); fwrite(summ, file.path(out, "catalogue_summary.tsv"), sep = "\t")
cat("total SNVs in catalogue:", nrow(cat_all), "\n")
gs <- cat_all[calibrated_context == TRUE & !grepl("&", mane_gene, fixed = TRUE), .(
  missense_catalogued = .N,
  PP3_CADD_Supporting = sum(PP3_CADD == "Supporting", na.rm = TRUE), PP3_CADD_Moderate = sum(PP3_CADD == "Moderate", na.rm = TRUE),
  PP3_AVI_Supporting = sum(PP3_AVI == "Supporting", na.rm = TRUE), PP3_AVI_Moderate = sum(PP3_AVI == "Moderate", na.rm = TRUE),
  observed_gnomAD = sum(gnomAD_status == "observed"), common_gt_1pct = sum(AF_max > 0.01, na.rm = TRUE),
  ClinVar_PLP = sum(clnsig %chin% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"), na.rm = TRUE),
  ClinVar_BLB = sum(clnsig %chin% c("Benign", "Likely_benign", "Benign/Likely_benign"), na.rm = TRUE)), by = .(gene = mane_gene)]
fwrite(gs[order(-missense_catalogued)], file.path(out, "gene_summary.tsv"), sep = "\t")
cat("genes in summary:", nrow(gs), "\n")
