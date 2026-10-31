# Paper 2 (analysis-plan-paper2-v1, "Reference sets"): unlabeled gnomAD v4.1 set for window sizing, per region
# class. Rare SNVs (0 < AF_max < 0.01, gnomAD exomes or genomes, VEP 116 cache) observed in gnomAD, absent from
# ClinVar 2026-09-23 (any record, as for the gnomAD sets of the main article), in disease genes (>= 1 P/LP
# record with >= 1 star in ClinVar 2026-09-23). Up to 100,000 SNVs per region class (all when fewer), seed
# 20261030. Region classes from paper2_regions.R. Noncoding candidate positions: every position for the
# canonical splice and intronic splice-region classes; for the other noncoding classes, 30,000 blocks of 20
# consecutive positions drawn at random (weighted by length) from the class intervals, then reclassified with
# the priority rules (blocks reduce the number of tabix queries). Coding classes (missense, synonymous): all
# MANE coding SNVs with gnomAD AF (data/coding/gnomad_af_mane_coding.rds).
# Output: C:/cadd/calibration/results/paper2/universe/gnomad_sample.rds and variants_gnomad.tsv
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
source("C:/cadd/calibration/R/paper2_regions.R")
set.seed(20261030)
out <- "C:/cadd/calibration/results/paper2/universe"; tmpd <- file.path(out, "tmp"); dir.create(tmpd, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
N_CLASS <- 100000L; N_BLOCKS <- 30000L; BLOCK <- 20L

## ClinVar 2026: disease genes and all record keys --------------------------------------------------
cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5, 8),
            col.names = c("chrom", "pos", "ref", "alt", "info"), colClasses = list(character = 1))
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
one_star <- grepl("practice_guideline|reviewed_by_expert_panel|criteria_provided,_multiple_submitters,_no_conflicts|criteria_provided,_single_submitter|criteria_provided,_conflicting", cv$revstat)   # >= 1 star (rules of scripts 13, 42 and 53)
dg <- unique(sub(":.*", "", unlist(strsplit(cv[clnsig %chin% plp & one_star & !is.na(geneinfo), geneinfo], "|", fixed = TRUE))))
cvk <- unique(cv[, ..k]); rm(cv, one_star); invisible(gc())
cat("Disease genes:", length(dg), "\n")

## Candidate noncoding positions ---------------------------------------------------------------------
m <- mane_structure()
txd <- m$txs[gene %chin% dg]
coding_pos <- unique(readRDS(file.path(CDIR, "mane_cds_positions.rds"))[, .(chrom, pos)])
intr <- m$introns[tx %chin% txd$tx]
all_pos <- function(iv) iv[, .(pos = seq.int(start, end)), by = .(chrom, i = seq_len(nrow(iv)))][, .(chrom, pos)]
splice_iv <- rbind(intr[, .(chrom, start, end = pmin(end, start + 7L))], intr[, .(chrom, start = pmax(start, end - 7L), end)])
utr_iv <- merge(m$exons[tx %chin% txd$tx], m$cds, by = "tx")[, .(chrom, start, end)]   # exons; CDS positions removed below
flank_iv <- txd[, .(chrom, start = c(pmax(1L, tstart - 5000L), tend + 1L), end = c(tstart - 1L, tend + 5000L))]
near_iv <- rbind(intr[end - start >= 16L, .(chrom, start = start + 8L, end = pmin(end, start + 49L))],
                 intr[end - start >= 16L, .(chrom, start = pmax(start, end - 49L), end = end - 8L)])
deep_iv <- intr[end - start > 100L, .(chrom, start = start + 50L, end = end - 50L)]
blocks <- function(iv, n) {
  iv <- iv[end - start + 1L >= BLOCK]; w <- as.numeric(iv$end - iv$start - BLOCK + 2L)
  i <- sample.int(nrow(iv), n, replace = TRUE, prob = w)
  st <- iv$start[i] + floor(runif(n) * w[i])
  data.table(chrom = iv$chrom[i], s = as.integer(st))[, .(pos = s + 0:(BLOCK - 1L)), by = .(chrom, b = seq_len(n))][, .(chrom, pos)]
}
cand <- unique(rbind(all_pos(splice_iv), blocks(utr_iv, N_BLOCKS), blocks(flank_iv, N_BLOCKS),
                     blocks(near_iv, N_BLOCKS), blocks(deep_iv, N_BLOCKS)))
cand <- cand[!coding_pos, on = .(chrom, pos)]
cls <- classify_noncoding(cand, m)[gene %chin% dg & class != "exonic_other" & class != "intergenic"]
cat("Candidate noncoding positions by class\n"); print(cls[, .N, by = class])

## gnomAD AF at candidate positions (VEP 116 cache), cached per chromosome ----------------------------
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
  r[, chrom := x$chrom[1]][x[, .(chrom, pos)], on = .(chrom, pos), nomatch = NULL]
}
chs <- unique(cls$chrom); todo <- chs[!file.exists(file.path(tmpd, sprintf("gn_af_%s.rds", chs)))]
if (length(todo)) {
  cl <- parallel::makeCluster(4L)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) }))
  parallel::clusterExport(cl, c("vep", "vcols", "af_chr", "tmpd"))
  parallel::parLapplyLB(cl, lapply(todo, function(ch) cls[chrom == ch]), function(x) { saveRDS(af_chr(x), file.path(tmpd, sprintf("gn_af_%s.rds", x$chrom[1]))); NULL })
  parallel::stopCluster(cl)
}
gn_nc <- rbindlist(lapply(chs, function(ch) readRDS(file.path(tmpd, sprintf("gn_af_%s.rds", ch)))))
gn_nc <- unique(gn_nc[AF_max > 0 & AF_max < 0.01], by = k)
gn_nc <- merge(gn_nc, cls, by = c("chrom", "pos"))[, .(chrom, pos, ref, alt, AF_max, gene, class)]

## Coding SNVs with gnomAD AF ----------------------------------------------------------------------
gc_ <- readRDS(file.path(CDIR, "gnomad_af_mane_coding.rds"))[AF_max > 0 & AF_max < 0.01, .(chrom, pos, ref, alt, AF_max)]
edge <- exon_edge_positions()
gn_c <- rbindlist(lapply(split(gc_, gc_$chrom), function(x) {
  r <- classify_coding_chr(x, edge); if (is.null(r)) return(NULL); merge(x, r, by = k) }))
gn_c <- gn_c[gene %chin% dg & class %chin% c("missense", "synonymous_exon_edge", "synonymous_other")]

## ClinVar exclusion and sampling ---------------------------------------------------------------------
gn <- rbind(gn_nc, gn_c)[!cvk, on = k]
cat("\nRare gnomAD SNVs absent from ClinVar, by class (before sampling)\n"); print(gn[, .N, by = class])
gn <- gn[gn[, .I[if (.N > N_CLASS) sample(.N, N_CLASS) else seq_len(.N)], by = class]$V1]
cat("\nSampled set by class\n"); print(gn[, .N, by = class][order(match(class, REGION_CLASSES))])
saveRDS(gn, file.path(out, "gnomad_sample.rds"))
u <- unique(gn[, ..k])[, chrom := fifelse(chrom == "MT", "M", chrom)]; setorder(u, chrom, pos)
fwrite(u[, .(variant = sprintf("chr%s:%d:%s>%s", chrom, pos, ref, alt))], file.path(out, "variants_gnomad.tsv"), sep = "\t")
cat("Unique SNVs to query:", nrow(u), "\n")
