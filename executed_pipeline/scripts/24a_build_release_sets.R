# Amendment A6 (a): reference sets from ClinVar GRCh38 releases 2019-12-02 and 2022-11-19 with the
# plan-v1 criteria (same code as 18_temporal_holdout.R v2). The 2024-01-27 and 2026-09-23 sets
# already exist (temporal_holdout/refsets_2024.rds; data/refsets/*).
# Output: C:/cadd/calibration/results/v17_avi/release_series/refsets_<release>.rds
suppressPackageStartupMessages(library(data.table))
cdir <- "C:/cadd/data/coding"; out <- "C:/cadd/calibration/results/v17_avi/release_series"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
blb <- c("Benign", "Likely_benign", "Benign/Likely_benign")
read_clinvar <- function(f) {
  cv <- fread(f, sep = "\t", skip = "#CHROM", select = c(1, 2, 3, 4, 5, 8),
              col.names = c("chrom", "pos", "clinvar_id", "ref", "alt", "info"), colClasses = list(character = c(1, 3)))
  getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
  cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
  cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                      grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                      grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
  dg <- cv[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo),
           unique(unlist(lapply(strsplit(geneinfo, "|", fixed = TRUE), function(z) sub(":.*", "", z))))]
  list(dg = dg, snv = cv[nchar(ref) == 1 & nchar(alt) == 1, .(chrom, pos, ref, alt, clinvar_id, clnsig, stars)], n = nrow(cv),
       revstat_values = cv[, .N, by = revstat][order(-N)][1:8])
}
rel <- c(r2019 = "C:/cadd/annot/clinvar_archive/clinvar_20191202.vcf.gz", r2022 = "C:/cadd/annot/clinvar_archive/clinvar_20221119.vcf.gz")
cvs <- lapply(rel, read_clinvar)
for (r in names(cvs)) { cat(sprintf("%s: %d records, %d disease genes\n", r, cvs[[r]]$n, length(cvs[[r]]$dg))); print(cvs[[r]]$revstat_values) }
dgu <- Reduce(union, lapply(cvs, `[[`, "dg"))

pos <- readRDS(file.path(cdir, "mane_cds_positions.rds"))
setorder(pos, transcript_id, cds_pos)
pos[, jump := c(abs(diff(pos)) != 1L, FALSE), by = transcript_id][, idx := seq_len(.N), by = transcript_id]
jn <- pos[jump == TRUE, .(transcript_id, kk = idx)]
flag <- rbindlist(lapply(c(-2:0, 1:3), function(o) jn[, .(transcript_id, idx = kk + o)]))
pos[, near := FALSE][flag, on = .(transcript_id, idx), near := TRUE]
edge <- unique(pos[near == TRUE, .(chrom, pos)])[, exon_edge := TRUE]
rm(pos, jn, flag); invisible(gc())
mis <- rbindlist(lapply(list.files(cdir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f) {
  d <- readRDS(f)[, .(chrom, pos, ref, alt, gene_name, consequence)]
  bad <- unique(d[consequence %chin% c("stop_gained", "stop_lost", "start_lost", "incomplete_terminal_codon_variant"), .(chrom, pos, alt)])
  d <- d[consequence == "missense"][!bad, on = .(chrom, pos, alt)][, consequence := NULL]
  d <- d[d[gene_name %chin% dgu, unique(.SD), .SDcols = c("chrom", "pos", "alt")], on = .(chrom, pos, alt)]
  r <- d[, c(list(gene = paste(sort(unique(gene_name)), collapse = "&")),
             lapply(cvs, function(cv) any(gene_name %chin% cv$dg))), by = .(chrom, pos, ref, alt)]
  rm(d); invisible(gc()); r
}))
mis[edge, on = .(chrom, pos), exon_edge := i.exon_edge][is.na(exon_edge), exon_edge := FALSE]
af <- readRDS(file.path(cdir, "gnomad_af_mane_coding.rds")); af <- af[, .(chrom, pos, ref, alt, AF_max)]; invisible(gc())
mis[af, on = .(chrom, pos, ref, alt), AF_max := i.AF_max][is.na(AF_max), AF_max := 0]
rm(af, edge); invisible(gc())

attach_scores <- function(d) rbindlist(lapply(split(d, d$chrom), function(x) {
  ch <- x$chrom[1]
  s1 <- readRDS(file.path(cdir, sprintf("cadd_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, CADD_PHRED = PHRED)]
  s2 <- readRDS(file.path(cdir, sprintf("avi_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, AVI_PHRED = round(AVI_PHRED, 2))]
  merge(merge(x, s1, by = k, all.x = TRUE), s2, by = k, all.x = TRUE)
}))
for (r in names(cvs)) {
  cv <- cvs[[r]]
  m <- mis[get(r) == TRUE]
  cvm <- merge(m, cv$snv, by = k, all.x = TRUE)
  clin <- cvm[!is.na(clnsig) & (clnsig %chin% plp | clnsig %chin% blb) & stars >= 1 & AF_max < 0.01][, y := as.integer(clnsig %chin% plp)]
  gno <- cvm[is.na(clinvar_id) & AF_max > 0 & AF_max < 0.01][!cv$snv, on = k]
  s <- list(clin = attach_scores(clin), gno = attach_scores(gno))
  cat(sprintf("%s sets: %d P/LP, %d B/LB (%d genes); gnomAD %d\n", r, s$clin[y == 1, .N], s$clin[y == 0, .N], uniqueN(s$clin$gene), nrow(s$gno)))
  saveRDS(s, file.path(out, sprintf("refsets_%s.rds", r)))
}
