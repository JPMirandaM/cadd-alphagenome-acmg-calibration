# Aim 9 (amendment A10, component 1): calibration of the noncoding region classes that reach
# >= 1,000 P/LP and >= 1,000 B/LB ClinVar reference variants (from script 30: only the intronic
# splice-region class, +/-3-8 nt from a MANE Select exon boundary). As registered in A10: CADD v1.7 and
# AVI, one seed (20260925), B = 10,000, increment 0.1; the gnomAD set is the rare (0 < AF < 0.01) SNVs
# observed in gnomAD v4.1 in the same region class and the same genes, absent from ClinVar (any record).
# Region definition as in script 30: intronic position at distance 3-8 from the nearest exon of a MANE
# Select transcript, not in MANE CDS and not at a canonical splice position (+/-1-2) of any MANE transcript.
# Output: C:/cadd/calibration/results/v17_avi/noncoding/splice_region_intronic/
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/pejaver_port_fast.R")
out <- "C:/cadd/calibration/results/v17_avi/noncoding/splice_region_intronic"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
tmpd <- file.path(out, "tmp"); dir.create(tmpd, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
ref <- readRDS("C:/cadd/calibration/results/v17_avi/noncoding/clinvar_noncoding_reference.rds")[class == "splice_region_intronic"]
genes <- unique(ref$gene)
cat(sprintf("ClinVar splice-region intronic reference: %d P/LP, %d B/LB, %d genes\n", sum(ref$y == 1), sum(ref$y == 0), length(genes)))

## Candidate positions ----------------------------------------------------------------------------
gtf <- fread("C:/cadd/annot/gencode.v50.basic.annotation.gtf.gz", sep = "\t", header = FALSE, skip = "chr1\t",
             col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
gtf <- gtf[grepl('tag "MANE_Select"', attr, fixed = TRUE) & type == "exon"]
gtf[, `:=`(chrom = sub("^chr", "", chr), tx = sub('.*transcript_id "([^"]+)".*', "\\1", attr), gene = sub('.*gene_name "([^"]+)".*', "\\1", attr))]
gtf <- gtf[!grepl("_PAR_Y$", tx)]
introns <- gtf[order(tx, start)][, .(start = end[-.N] + 1L, end = start[-1] - 1L, gene = gene[1]), by = .(chrom, tx)][start <= end]
mk <- function(d) introns[, .(chrom, gene, pos = c(start + d - 1L, end - d + 1L))]
canon <- unique(rbind(mk(1L), mk(2L))[, .(chrom, pos)])
sr <- unique(rbindlist(lapply(3:8, mk))[gene %chin% genes])
coding_pos <- unique(readRDS("C:/cadd/data/coding/mane_cds_positions.rds")[, .(chrom, pos)])
sr <- sr[!canon, on = .(chrom, pos)][!coding_pos, on = .(chrom, pos)]
sr <- unique(sr, by = c("chrom", "pos"))
cat("Candidate splice-region intronic positions in the reference genes:", nrow(sr), "\n")

## gnomAD AF (VEP cache), ClinVar exclusion, scores ------------------------------------------------
vep <- "C:/cadd/vep/homo_sapiens/116_GRCh38"
vcols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", readLines(file.path(vep, "info.txt")), value = TRUE)), ",")[[1]]
af_chr <- function(x) {
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
  r[, chrom := x$chrom[1]][x[, .(chrom, pos, gene)], on = .(chrom, pos), nomatch = NULL]
}
lookup_chr <- function(file, pref, x, col, round2 = FALSE) {
  p <- sort(unique(x$pos)); grp <- cumsum(c(1L, diff(p) > 1L)); rg <- data.table(s = p, g = grp)[, .(start = min(s), end = max(s)), by = g]
  tf <- open(TabixFile(file)); on.exit(close(tf))
  r <- rbindlist(lapply(split(rg, ceiling(seq_len(nrow(rg)) / 5000)), function(b) {
    l <- unlist(scanTabix(tf, param = GenomicRanges::GRanges(paste0(pref, x$chrom[1]), IRanges::IRanges(b$start, b$end))), use.names = FALSE)
    if (!length(l)) return(NULL); fread(text = l, sep = "\t", header = FALSE, select = c(2, 3, 4, 6), col.names = c("pos", "ref", "alt", "sc"), nThread = 1) }))
  if (!nrow(r)) return(NULL)
  if (round2) r[, sc := round(sc, 2)]; r[, chrom := x$chrom[1]]; setnames(r, "sc", col); r[x[, ..k], on = k, nomatch = NULL] }
run_cached <- function(tag, d, fun) {
  chs <- unique(d$chrom); todo <- chs[!file.exists(file.path(tmpd, sprintf("%s_%s.rds", tag, chs)))]
  if (length(todo)) {
    cl <- parallel::makeCluster(4L); on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) }))
    parallel::clusterExport(cl, c("k", "vep", "vcols", "lookup_chr", "af_chr", "tmpd"), envir = globalenv())
    parallel::clusterExport(cl, c("tag", "fun"), envir = environment())
    parallel::parLapplyLB(cl, lapply(todo, function(ch) d[chrom == ch]), function(x) { saveRDS(fun(x), file.path(tmpd, sprintf("%s_%s.rds", tag, x$chrom[1]))); NULL })
  }
  rbindlist(lapply(chs, function(ch) readRDS(file.path(tmpd, sprintf("%s_%s.rds", tag, ch)))), fill = TRUE)
}
gn <- run_cached("sr_af", sr, af_chr)
gn <- gn[AF_max > 0 & AF_max < 0.01]
cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5), col.names = k, colClasses = list(character = 1))
gn <- gn[!cv, on = k]; rm(cv); invisible(gc())
cat("Rare gnomAD SNVs at splice-region intronic positions, absent from ClinVar:", nrow(gn), "\n")
sn_cadd <- seqnamesTabix(TabixFile("C:/cadd/data/whole_genome_SNVs.tsv.gz"))
sn_avi <- sub("^chr", "", seqnamesTabix(TabixFile("C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz")))
s1 <- run_cached("sr_cadd", gn[chrom %chin% sn_cadd], function(x) lookup_chr("C:/cadd/data/whole_genome_SNVs.tsv.gz", "", x, "CADD_PHRED"))
s2 <- run_cached("sr_avi", gn[chrom %chin% sn_avi], function(x) lookup_chr("C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz", "chr", x, "AVI_PHRED", TRUE))
gn <- merge(merge(gn, s1, by = k), s2, by = k)
saveRDS(gn, file.path(out, "gnomad_set.rds"))
cat("gnomAD set with both scores:", nrow(gn), "\n")

## Calibration (A10: one seed, B = 10,000, increment 0.1) -----------------------------------------
res_all <- list()
for (tn in c("CADD", "AVI")) {
  col <- if (tn == "CADD") "CADD_PHRED" else "AVI_PHRED"
  o <- file.path(out, tn); dir.create(o, showWarnings = FALSE)
  if (!file.exists(file.path(o, "thresholds.tsv"))) {
    t0 <- Sys.time()
    r <- calibrate_fast(ref[[col]], ref$y, gn[[col]], B = 10000L, workers = 6L, seed = 20260925L)
    fwrite(r$thresholds, file.path(o, "thresholds.tsv"), sep = "\t"); fwrite(r$curves, file.path(o, "posterior_curves.tsv"), sep = "\t")
    fwrite(r$boot_thresholds, file.path(o, "bootstrap_thresholds.tsv"), sep = "\t"); writeLines(capture.output(str(r$params)), file.path(o, "parameters.txt"))
    cat(sprintf("%s calibrated in %.1f min\n", tn, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  res_all[[tn]] <- data.table(tool = tn, fread(file.path(o, "thresholds.tsv"), colClasses = c(level = "character")))
}
thr <- rbindlist(res_all); fwrite(thr, file.path(out, "thresholds.tsv"), sep = "\t")
cat("\n== Splice-region intronic thresholds (one seed) ==\n"); print(thr)
main <- fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main", .(PP3_missense = median(PP3_discounted), BP4_missense = median(BP4_discounted)), by = .(tool, level)]
cmp <- merge(thr[, .(tool, level, PP3_splice_region = PP3_discounted, BP4_splice_region = BP4_discounted)], main, by = c("tool", "level"))
fwrite(cmp, file.path(out, "comparison_with_missense.tsv"), sep = "\t"); cat("\n== Comparison with missense thresholds ==\n"); print(cmp)
