# Aim 9 (amendments A10/A11), components (1) ClinVar by region and (2) MPRA.
# (1) Region class of every ClinVar 2026-09-23 SNV relative to MANE Select transcripts; reference
#     variants per class (plan-v1 rules); class sizes; evaluation of the missense-calibrated CADD v1.7
#     and AVI thresholds in each class (interval likelihood ratios with 1,000 bootstrap replicates,
#     AUROC). Class priority when a SNV falls in several transcripts: canonical splice > splice region
#     intronic > 5'UTR > 3'UTR > near-intronic > deep intronic > upstream > downstream.
# (2) MPRA saturation mutagenesis (Kircher et al. 2019): functional = pValue < 1e-5 and |Coefficient| >= 0.1.
# Output: C:/cadd/calibration/results/v17_avi/noncoding/
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
out <- "C:/cadd/calibration/results/v17_avi/noncoding"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
cdir <- "C:/cadd/data/coding"; k <- c("chrom", "pos", "ref", "alt")
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"); blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")

## MANE transcript structure -----------------------------------------------------------------------
gtf <- fread("C:/cadd/annot/gencode.v50.basic.annotation.gtf.gz", sep = "\t", header = FALSE, skip = "chr1\t",
             col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
gtf <- gtf[grepl('tag "MANE_Select"', attr, fixed = TRUE) & type %chin% c("transcript", "exon", "CDS")]
gtf[, `:=`(chrom = sub("^chr", "", chr), tx = sub('.*transcript_id "([^"]+)".*', "\\1", attr), gene = sub('.*gene_name "([^"]+)".*', "\\1", attr))]
gtf <- gtf[!grepl("_PAR_Y$", tx)][chrom == "M", chrom := "MT"][, c("chr", "src", "score", "frame", "attr") := NULL]
txs <- gtf[type == "transcript", .(chrom, tx, gene, strand, tstart = start, tend = end)]
exons <- gtf[type == "exon", .(chrom, tx, start, end)]
cds <- gtf[type == "CDS", .(cds_start = min(start), cds_end = max(end)), by = tx]
introns <- exons[order(tx, start)][, .(start = end[-.N] + 1L, end = start[-1] - 1L), by = .(chrom, tx)][start <= end]

## ClinVar SNVs outside MANE CDS ---------------------------------------------------------------------
cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 3, 4, 5, 8),
            col.names = c("chrom", "pos", "clinvar_id", "ref", "alt", "info"), colClasses = list(character = c(1, 3)))
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                    grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                    grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
dg <- cv[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo), unique(unlist(lapply(strsplit(geneinfo, "|", fixed = TRUE), function(z) sub(":.*", "", z))))]
ref <- cv[nchar(ref) == 1 & nchar(alt) == 1 & stars >= 1 & (clnsig %chin% plp | clnsig %chin% blb), .(chrom, pos, ref, alt, y = as.integer(clnsig %chin% plp))]
rm(cv); invisible(gc())
coding_pos <- unique(readRDS(file.path(cdir, "mane_cds_positions.rds"))[, .(chrom, pos)])
nc <- ref[!coding_pos, on = .(chrom, pos)]
cat("ClinVar P/LP+B/LB SNVs (>=1 star):", nrow(ref), "| outside MANE CDS:", nrow(nc), "\n")

## Region classes ---------------------------------------------------------------------------------
q <- unique(nc[, .(chrom, pos)])[, `:=`(s = pos, e = pos)]
setkey(q, chrom, s, e)
cls <- list()
ex <- merge(merge(exons, cds, by = "tx"), txs[, .(tx, strand, gene)], by = "tx")   # inner joins: transcripts with exons and CDS only
setnames(ex, c("start", "end"), c("s", "e")); ex <- ex[!is.na(s) & !is.na(e)]; setkey(ex, chrom, s, e)
o <- foverlaps(q, ex, nomatch = NULL)[, .(chrom, pos, tx, gene, class = fifelse((strand == "+" & pos < cds_start) | (strand == "-" & pos > cds_end), "5UTR",
                                        fifelse((strand == "+" & pos > cds_end) | (strand == "-" & pos < cds_start), "3UTR", "exonic_other")))]
cls[[1]] <- o
it <- introns[txs[, .(tx, gene)], on = "tx", nomatch = NULL]; setnames(it, c("start", "end"), c("s", "e")); it <- it[!is.na(s) & !is.na(e)]; setkey(it, chrom, s, e)
o <- foverlaps(q, it, nomatch = NULL)[, d := pmin(pos - s + 1L, e - pos + 1L)]
cls[[2]] <- o[, .(chrom, pos, tx, gene, class = fcase(d <= 2, "canonical_splice", d <= 8, "splice_region_intronic", d <= 50, "near_intronic", default = "deep_intronic"))]
fl <- txs[, .(chrom, tx, gene, strand, s = pmax(1L, tstart - 5000L), e = tend + 5000L, tstart, tend)]; setkey(fl, chrom, s, e)
o <- foverlaps(q, fl, nomatch = NULL)[pos < tstart | pos > tend]
cls[[3]] <- o[, .(chrom, pos, tx, gene, class = fifelse((strand == "+" & pos < tstart) | (strand == "-" & pos > tend), "upstream_5kb", "downstream_5kb"))]
pri <- c("canonical_splice", "splice_region_intronic", "5UTR", "3UTR", "near_intronic", "deep_intronic", "upstream_5kb", "downstream_5kb", "exonic_other")
reg <- rbindlist(cls)[, p := match(class, pri)][order(chrom, pos, p)][, .SD[1], by = .(chrom, pos)][, p := NULL]
nc <- merge(nc, reg, by = c("chrom", "pos"), all.x = TRUE)[is.na(class), `:=`(class = "intergenic", gene = NA_character_)]
# Synonymous (coding) from the consequence tables, split by exon edge
pos_all <- readRDS(file.path(cdir, "mane_cds_positions.rds")); setorder(pos_all, transcript_id, cds_pos)
pos_all[, jump := c(abs(diff(pos)) != 1L, FALSE), by = transcript_id][, idx := seq_len(.N), by = transcript_id]
jn <- pos_all[jump == TRUE, .(transcript_id, kk = idx)]
fl2 <- rbindlist(lapply(c(-2:0, 1:3), function(o) jn[, .(transcript_id, idx = kk + o)]))
pos_all[, near := FALSE][fl2, on = .(transcript_id, idx), near := TRUE]
edge <- unique(pos_all[near == TRUE, .(chrom, pos)]); rm(pos_all); invisible(gc())
syn <- rbindlist(lapply(list.files(cdir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f) {
  d <- readRDS(f)[, .(chrom, pos, ref, alt, gene_name, consequence)][ref, on = k, nomatch = NULL]
  d[, .(gene = gene_name[1], only_syn = all(consequence == "synonymous"), y = y[1]), by = k][only_syn == TRUE][, only_syn := NULL]
}))
syn[, class := "synonymous_other"][edge, on = .(chrom, pos), class := "synonymous_exon_edge"]
allc <- rbind(nc[class != "exonic_other", .(chrom, pos, ref, alt, y, gene, class)], syn[, .(chrom, pos, ref, alt, y, gene, class)])
allc <- allc[!is.na(gene) & gene %chin% dg | class == "intergenic"]
cat("ClinVar noncoding/synonymous reference SNVs (disease genes; intergenic kept for description):", nrow(allc), "\n")

## Scores (tabix) and gnomAD AF (VEP cache) -------------------------------------------------------
# Lookups run per chromosome in 4 parallel workers; each chromosome's result is cached in out/tmp
# so that an interrupted run resumes where it stopped (added after a 4-hour run was lost when the
# session ended). Results are identical to the sequential version.
tmpd <- file.path(out, "tmp"); dir.create(tmpd, showWarnings = FALSE)
vep <- "C:/cadd/vep/homo_sapiens/116_GRCh38"
vcols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", readLines(file.path(vep, "info.txt")), value = TRUE)), ",")[[1]]
lookup_chr <- function(file, pref, x, col, round2 = FALSE) {
  p <- sort(unique(x$pos)); grp <- cumsum(c(1L, diff(p) > 1L)); rg <- data.table(s = p, g = grp)[, .(start = min(s), end = max(s)), by = g]
  tf <- open(TabixFile(file)); on.exit(close(tf))
  r <- rbindlist(lapply(split(rg, ceiling(seq_len(nrow(rg)) / 5000)), function(b) {
    l <- unlist(scanTabix(tf, param = GenomicRanges::GRanges(paste0(pref, x$chrom[1]), IRanges::IRanges(b$start, b$end))), use.names = FALSE)
    if (!length(l)) return(NULL); fread(text = l, sep = "\t", header = FALSE, select = c(2, 3, 4, 6), col.names = c("pos", "ref", "alt", "sc"), nThread = 1) }))
  if (!nrow(r)) return(NULL)
  if (round2) r[, sc := round(sc, 2)]; r[, chrom := x$chrom[1]]; setnames(r, "sc", col); r[x[, ..k], on = k, nomatch = NULL] }
run_cached <- function(tag, d, fun) {
  chs <- unique(d$chrom)
  todo <- chs[!file.exists(file.path(tmpd, sprintf("%s_%s.rds", tag, chs)))]
  if (length(todo)) {
    cl <- parallel::makeCluster(4L); on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) }))
    parallel::clusterExport(cl, c("k", "vep", "vcols", "lookup_chr", "af_chr", "tmpd"), envir = globalenv())
    parallel::clusterExport(cl, c("tag", "fun"), envir = environment())
    parallel::parLapplyLB(cl, lapply(todo, function(ch) d[chrom == ch]), function(x) {
      r <- fun(x); saveRDS(r, file.path(tmpd, sprintf("%s_%s.rds", tag, x$chrom[1]))); NULL })
    cat(sprintf("%s: %d chromosomes computed (%s)\n", tag, length(todo), format(Sys.time(), "%H:%M")))
  }
  rbindlist(lapply(chs, function(ch) readRDS(file.path(tmpd, sprintf("%s_%s.rds", tag, ch)))), fill = TRUE)
}
af_chr <- function(x) {
  f <- file.path(vep, x$chrom[1], "all_vars.gz"); if (!file.exists(f)) return(NULL)
  p <- sort(unique(x$pos)); grp <- cumsum(c(1L, diff(p) > 1L))
  rg <- data.table(s = p, g = grp)[, .(start = min(s), end = max(s)), by = g]; tf <- open(TabixFile(f, index = paste0(f, ".csi"))); on.exit(close(tf))
  rbindlist(lapply(split(rg, ceiling(seq_len(nrow(rg)) / 5000)), function(b) {
    l <- unlist(scanTabix(tf, param = GenomicRanges::GRanges(x$chrom[1], IRanges::IRanges(b$start, b$end))), use.names = FALSE)
    if (!length(l)) return(NULL); sel <- which(vcols %in% c("start", "allele_string", "gnomADe", "gnomADg"))
    d <- fread(text = l, sep = "\t", header = FALSE, colClasses = "character", select = sel, col.names = vcols[sel], nThread = 1)
    d <- d[grepl("^[ACGT](/[ACGT])+$", allele_string) & (gnomADe != "." | gnomADg != ".")]; if (!nrow(d)) return(NULL)
    lg <- rbindlist(lapply(c("gnomADe", "gnomADg"), function(kk) { y <- d[get(kk) != ".", .(pos = as.integer(start), allele_string, v = get(kk))]
      if (!nrow(y)) return(NULL); y <- y[, .(e = unlist(strsplit(v, ",", fixed = TRUE))), by = .(pos, allele_string)]
      y[, .(pos, ref = sub("/.*", "", allele_string), alt = sub(":.*", "", e), af = as.numeric(sub(".*:", "", e)))] }))
    lg[, .(AF_max = max(af)), by = .(pos, ref, alt)][, chrom := x$chrom[1]] })) }
sn_cadd <- seqnamesTabix(TabixFile("C:/cadd/data/whole_genome_SNVs.tsv.gz"))                   # e.g. no MT in the CADD file
sn_avi <- sub("^chr", "", seqnamesTabix(TabixFile("C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz")))
sc1 <- run_cached("cv_cadd", allc[chrom %chin% sn_cadd], function(x) lookup_chr("C:/cadd/data/whole_genome_SNVs.tsv.gz", "", x, "CADD_PHRED"))
sc2 <- run_cached("cv_avi", allc[chrom %chin% sn_avi], function(x) lookup_chr("C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz", "chr", x, "AVI_PHRED", TRUE))
allc <- merge(merge(allc, sc1, by = k, all.x = TRUE), sc2, by = k, all.x = TRUE)
af <- run_cached("cv_af", allc, af_chr)
allc <- merge(allc, af[, .(AF_max = max(AF_max)), by = k], by = k, all.x = TRUE)[is.na(AF_max), AF_max := 0]
allc <- allc[AF_max < 0.01 & !is.na(CADD_PHRED) & !is.na(AVI_PHRED)]
saveRDS(allc, file.path(out, "clinvar_noncoding_reference.rds"))
sizes <- allc[, .(P_LP = sum(y == 1), B_LB = sum(y == 0), genes = uniqueN(gene)), by = class][order(match(class, c(pri, "intergenic", "synonymous_exon_edge", "synonymous_other")))]
sizes[, calibratable := P_LP >= 1000 & B_LB >= 1000]
cat("\n== Reference SNVs per class (rare, >= 1 star) ==\n"); print(sizes); fwrite(sizes, file.path(out, "class_sizes.tsv"), sep = "\t")

## Evaluation of missense-calibrated thresholds per class ----------------------------------------
thr <- fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main",
          .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)]
req <- c(Supporting = 2.406, Moderate = 5.790, Strong = 33.53)
auroc <- function(s, y) { r <- rank(s); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0)); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
ivl <- function(s, t) { lv <- rep("Indeterminate", length(s))
  for (L in c("Supporting", "Moderate", "Strong")) { v <- t[level == L, PP3]; if (length(v) && !is.na(v)) lv[s >= v] <- paste0("PP3_", L) }
  for (L in c("Supporting", "Moderate", "Strong")) { v <- t[level == L, BP4]; if (length(v) && !is.na(v)) lv[s <= v] <- paste0("BP4_", L) }; lv }
evalc <- rbindlist(lapply(sizes[P_LP >= 20 & B_LB >= 20, class], function(cl) rbindlist(lapply(c(CADD = "CADD_PHRED", AVI = "AVI_PHRED"), function(col) {
  x <- allc[class == cl]; tn <- if (col == "CADD_PHRED") "CADD" else "AVI"; iv <- ivl(x[[col]], thr[tool == tn]); yy <- x$y
  lr <- function(iv, yy) sapply(sort(unique(iv)), function(L) mean(iv[yy == 1] == L) / mean(iv[yy == 0] == L))
  e <- lr(iv, yy); set.seed(20260925)
  bs <- replicate(1000, { i <- sample.int(length(yy), replace = TRUE); v <- lr(iv[i], yy[i]); v[names(e)] })
  data.table(class = cl, tool = tn, AUROC = auroc(x[[col]], yy), interval = names(e), n_PLP = as.integer(table(factor(iv[yy == 1], levels = names(e)))),
             n_BLB = as.integer(table(factor(iv[yy == 0], levels = names(e)))), LR = e,
             lo = apply(bs, 1, quantile, 0.025, na.rm = TRUE), hi = apply(bs, 1, quantile, 0.975, na.rm = TRUE))
}))))
evalc[, lvl := sub("^(PP3|BP4)_", "", interval)]
evalc[grepl("^PP3", interval), meets := lo >= req[lvl]][grepl("^BP4", interval), meets := hi <= 1 / req[lvl]]
cat("\n== AUROC per class ==\n"); print(unique(evalc[, .(class, tool, AUROC)])[order(class, tool)], digits = 3)
cat("\n== Interval LRs of missense thresholds per class (meets = conservative 95% bound meets required LR) ==\n")
print(evalc[interval != "Indeterminate", .(class, tool, interval, n_PLP, n_BLB, LR = signif(LR, 3), lo = signif(lo, 3), hi = signif(hi, 3), meets)][order(class, tool, interval)])
fwrite(evalc, file.path(out, "evaluation_missense_thresholds_by_class.tsv"), sep = "\t")

## MPRA -------------------------------------------------------------------------------------------
m <- fread("C:/cadd/annot/mpra/elements.tsv.gz")[Release == "GRCh38" & Ref %chin% c("A", "C", "G", "T") & Alt %chin% c("A", "C", "G", "T")]
m <- m[, .(chrom = as.character(Chrom), pos = Pos, ref = Ref, alt = Alt, Coefficient, pValue, Element)]
m[, functional := as.integer(pValue < 1e-5 & abs(Coefficient) >= 0.1)]
s1 <- run_cached("mpra_cadd", m[chrom %chin% sn_cadd], function(x) lookup_chr("C:/cadd/data/whole_genome_SNVs.tsv.gz", "", x, "CADD_PHRED"))
s2 <- run_cached("mpra_avi", m[chrom %chin% sn_avi], function(x) lookup_chr("C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz", "chr", x, "AVI_PHRED", TRUE))
m <- merge(merge(m, s1, by = k, all.x = TRUE), s2, by = k, all.x = TRUE)[!is.na(CADD_PHRED) & !is.na(AVI_PHRED)]
saveRDS(m, file.path(out, "mpra_scored.rds"))
mres <- rbind(data.table(Element = "ALL", n = nrow(m), functional = sum(m$functional), CADD_AUROC = auroc(m$CADD_PHRED, m$functional), AVI_AUROC = auroc(m$AVI_PHRED, m$functional)),
              m[, .(n = .N, functional = sum(functional), CADD_AUROC = if (sum(functional) >= 10 & sum(functional == 0) >= 10) auroc(CADD_PHRED, functional) else NA_real_,
                    AVI_AUROC = if (sum(functional) >= 10 & sum(functional == 0) >= 10) auroc(AVI_PHRED, functional) else NA_real_), by = Element])
set.seed(20260925)
bsd <- replicate(1000, { i <- sample.int(nrow(m), replace = TRUE); auroc(m$AVI_PHRED[i], m$functional[i]) - auroc(m$CADD_PHRED[i], m$functional[i]) })
cat("\n== MPRA: AUROC for functional (p < 1e-5, |coef| >= 0.1) vs non-functional SNVs ==\n"); print(mres, digits = 3)
cat(sprintf("AVI - CADD AUROC (all elements): %.3f [95%% CI %.3f to %.3f]\n", mres[Element == "ALL", AVI_AUROC - CADD_AUROC], quantile(bsd, 0.025), quantile(bsd, 0.975)))
fwrite(mres, file.path(out, "mpra_auroc.tsv"), sep = "\t")
