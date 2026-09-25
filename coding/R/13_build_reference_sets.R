# Builds the calibration reference sets exactly as pre-specified in STUDY_PLAN.md (analysis plan v1):
#   ClinVar set: MANE missense SNVs, P/LP vs B/LB (exact CLNSIG), >= 1 star, rare (max gnomAD v4.1
#                exome/genome AF < 0.01; absent = rare), in disease genes.
#   gnomAD set : MANE missense SNVs observed in gnomAD v4.1 with AF < 0.01, absent from ClinVar,
#                in the same disease genes (window sizing only).
# Also flags SNVs within 3 nt of an internal exon boundary (sensitivity S4) and attaches
# CADD v1.7 PHRED and AVI Phred (rounded to 2 decimals).
# Output: C:/cadd/data/refsets/{clinvar_ref.rds, gnomad_ref.rds, summary.tsv}
# Usage: Rscript 13_build_reference_sets.R [workers]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 6L
cdir <- "C:/cadd/data/coding"; odir <- "C:/cadd/data/refsets"
dir.create(odir, recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat(sprintf(...), "\n")

## 1. ClinVar: all records (for disease genes) and SNVs ------------------------------------------
cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 3, 4, 5, 8),
            col.names = c("chrom", "pos", "clinvar_id", "ref", "alt", "info"), colClasses = list(character = c(1, 3)))
getk <- function(x, k) { v <- sub(paste0(".*(^|;)", k, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", k, "="), x)] <- NA; v }
cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
cv[, stars := fcase(grepl("practice_guideline", revstat), 4L,
                    grepl("reviewed_by_expert_panel", revstat), 3L,
                    grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                    grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L,
                    default = 0L)]
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")
dg <- cv[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo), unique(unlist(lapply(strsplit(geneinfo, "|", fixed = TRUE),
                                                                            function(z) sub(":.*", "", z))))]
log("ClinVar records: %d | disease genes (>=1 P/LP any type, >=1 star): %d", nrow(cv), length(dg))
cv_snv <- cv[nchar(ref) == 1 & nchar(alt) == 1, .(chrom, pos, ref, alt, clinvar_id, clnsig, stars)]
rm(cv); invisible(gc())

## 2. MANE missense SNVs (variant level) + exon-boundary flag ----------------------------------
pos <- readRDS(file.path(cdir, "mane_cds_positions.rds"))
setorder(pos, transcript_id, cds_pos)
pos[, jump := c(abs(diff(pos)) != 1L, FALSE), by = transcript_id]           # junction after this nt
pos[, near_junction := FALSE]
pos[, idx := seq_len(.N), by = transcript_id]
jn <- pos[jump == TRUE, .(transcript_id, k = idx)]
flag <- rbindlist(lapply(c(-2:0, 1:3), function(o) jn[, .(transcript_id, idx = k + o)]))
pos[flag, on = .(transcript_id, idx), near_junction := TRUE]
edge <- unique(pos[near_junction == TRUE, .(chrom, pos)])[, exon_edge := TRUE]
log("coding positions within 3 nt of an internal exon boundary: %d", nrow(edge))
rm(pos, flag, jn); invisible(gc())

# Memory-lean (v2, 2026-09-25): one chromosome at a time; rows restricted to disease-gene
# transcripts before any string aggregation; only the needed columns are kept.
n_mis_total <- 0L
mis <- rbindlist(lapply(list.files(cdir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f) {
  d <- readRDS(f)[, .(chrom, pos, ref, alt, gene_name, consequence, hgvsp)]
  bad <- unique(d[consequence %chin% c("stop_gained", "stop_lost", "start_lost", "incomplete_terminal_codon_variant"),
                  .(chrom, pos, alt)])
  d <- d[consequence == "missense"][!bad, on = .(chrom, pos, alt)]
  n_mis_total <<- n_mis_total + uniqueN(d, by = c("chrom", "pos", "alt"))
  keys <- unique(d[gene_name %chin% dg, .(chrom, pos, alt)])       # SNV missense in >= 1 disease gene
  d <- d[keys, on = .(chrom, pos, alt)]
  r <- d[, .(gene = paste(sort(unique(gene_name)), collapse = "&"), hgvsp = hgvsp[1]), by = .(chrom, pos, ref, alt)]
  rm(d); invisible(gc()); r
}))
mis <- merge(mis, edge, by = c("chrom", "pos"), all.x = TRUE)[is.na(exon_edge), exon_edge := FALSE]
mis[, in_disease_gene := TRUE]
rm(edge); invisible(gc())
log("MANE missense SNVs: %s | in disease genes: %s", format(n_mis_total, big.mark = ","),
    format(nrow(mis), big.mark = ","))

## 3. gnomAD v4.1 AF from the VEP 116 cache, over CDS regions of disease genes ------------------
reg <- unique(mis[in_disease_gene == TRUE, .(chrom, pos)])
setorder(reg, chrom, pos)
reg[, grp := cumsum(c(1L, diff(pos) != 1L) | c(TRUE, chrom[-1] != chrom[-.N]))]
reg <- reg[, .(start = min(pos), end = max(pos)), by = .(chrom, grp)]
reg[, batch := paste(chrom, (seq_len(.N) - 1L) %/% 2000L), by = chrom]
vep_dir <- "C:/cadd/vep/homo_sapiens/116_GRCh38"
cols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", readLines(file.path(vep_dir, "info.txt")), value = TRUE)), ",")[[1]]
cl <- makeCluster(workers)
clusterExport(cl, c("vep_dir", "cols"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })))
t0 <- Sys.time()
af <- rbindlist(parLapplyLB(cl, split(reg, reg$batch), function(b) {
  f <- file.path(vep_dir, b$chrom[1], "all_vars.gz")
  l <- unlist(scanTabix(TabixFile(f, index = paste0(f, ".csi")),
                        param = GenomicRanges::GRanges(b$chrom, IRanges::IRanges(b$start, b$end))), use.names = FALSE)
  if (!length(l)) return(NULL)
  sel <- which(cols %in% c("chr", "variation_name", "start", "end", "allele_string", "gnomADe", "gnomADg"))
  d <- fread(text = l, sep = "\t", header = FALSE, colClasses = "character", nThread = 1,
             select = sel, col.names = cols[sel])
  # SNVs identified by single-base alleles (the cache stores end = "." when end == start)
  d <- d[grepl("^[ACGT](/[ACGT])+$", allele_string) & (gnomADe != "." | gnomADg != ".")]
  if (!nrow(d)) return(NULL)
  expand <- function(field) {
    x <- d[get(field) != ".", .(chr, pos = as.integer(start), allele_string, v = get(field))]
    if (!nrow(x)) return(data.table(chrom = character(), pos = integer(), ref = character(), alt = character(), af = numeric()))
    x <- x[, .(e = unlist(strsplit(v, ",", fixed = TRUE))), by = .(chr, pos, allele_string)]
    x[, .(chrom = chr, pos, ref = sub("/.*", "", allele_string), alt = sub(":.*", "", e), af = as.numeric(sub(".*:", "", e)))]
  }
  e <- expand("gnomADe"); g <- expand("gnomADg")
  setnames(e, "af", "AF_exomes"); setnames(g, "af", "AF_genomes")
  merge(e[, .(AF_exomes = max(AF_exomes)), by = .(chrom, pos, ref, alt)],
        g[, .(AF_genomes = max(AF_genomes)), by = .(chrom, pos, ref, alt)], by = c("chrom", "pos", "ref", "alt"), all = TRUE)
}))
stopCluster(cl)
af <- af[nchar(alt) == 1 & nchar(ref) == 1]
af <- af[, .(AF_exomes = max(AF_exomes, na.rm = TRUE), AF_genomes = max(AF_genomes, na.rm = TRUE)), by = .(chrom, pos, ref, alt)]
af[is.infinite(AF_exomes), AF_exomes := NA][is.infinite(AF_genomes), AF_genomes := NA]
af[, AF_max := pmax(AF_exomes, AF_genomes, na.rm = TRUE)]
log("gnomAD SNV allele frequencies from VEP cache: %s rows in %.1f min", format(nrow(af), big.mark = ","),
    as.numeric(difftime(Sys.time(), t0, units = "mins")))

## 4. Reference sets ----------------------------------------------------------------------------
m <- merge(mis[in_disease_gene == TRUE], af, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE)
m[is.na(AF_max), AF_max := 0]
cvm <- merge(m, cv_snv[, .(chrom, pos, ref, alt, clinvar_id, clnsig, stars)], by = c("chrom", "pos", "ref", "alt"), all.x = TRUE)

clin <- cvm[!is.na(clnsig) & (clnsig %chin% plp | clnsig %chin% blb) & stars >= 1 & AF_max < 0.01]
clin[, y := as.integer(clnsig %chin% plp)]
gno  <- cvm[is.na(clinvar_id) & AF_max > 0 & AF_max < 0.01]
gno  <- gno[!cv_snv, on = .(chrom, pos, ref, alt)]

attach_scores <- function(d) {
  rbindlist(lapply(split(d, d$chrom), function(x) {
    ch <- x$chrom[1]
    s1 <- readRDS(file.path(cdir, sprintf("cadd_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, CADD_PHRED = PHRED)]
    s2 <- readRDS(file.path(cdir, sprintf("avi_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, AVI_PHRED = round(AVI_PHRED, 2))]
    merge(merge(x, s1, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE), s2, by = c("chrom", "pos", "ref", "alt"), all.x = TRUE)
  }))
}
clin <- attach_scores(clin); gno <- attach_scores(gno)
stopifnot(!anyNA(clin$CADD_PHRED), !anyNA(clin$AVI_PHRED), !anyNA(gno$CADD_PHRED), !anyNA(gno$AVI_PHRED))

saveRDS(clin, file.path(odir, "clinvar_ref.rds")); saveRDS(gno, file.path(odir, "gnomad_ref.rds"))
summ <- rbind(
  data.table(set = "ClinVar main (>=1 star, AF<0.01)", P_LP = clin[y == 1, .N], B_LB = clin[y == 0, .N], genes = uniqueN(clin$gene)),
  data.table(set = "S1: >=2 stars", P_LP = clin[y == 1 & stars >= 2, .N], B_LB = clin[y == 0 & stars >= 2, .N], genes = uniqueN(clin[stars >= 2, gene])),
  data.table(set = "S2: AF<0.001", P_LP = clin[y == 1 & AF_max < 0.001, .N], B_LB = clin[y == 0 & AF_max < 0.001, .N], genes = uniqueN(clin[AF_max < 0.001, gene])),
  data.table(set = "S4: no exon-edge", P_LP = clin[y == 1 & !exon_edge, .N], B_LB = clin[y == 0 & !exon_edge, .N], genes = uniqueN(clin[exon_edge == FALSE, gene])),
  data.table(set = "gnomAD set", P_LP = NA_integer_, B_LB = NA_integer_, genes = uniqueN(gno$gene))[, gnomAD := nrow(gno)],
  fill = TRUE)
print(summ)
fwrite(summ, file.path(odir, "summary.tsv"), sep = "\t")
log("P/LP exon-edge fraction: %.2f%% | B/LB exon-edge fraction: %.2f%%",
    100 * clin[y == 1, mean(exon_edge)], 100 * clin[y == 0, mean(exon_edge)])
