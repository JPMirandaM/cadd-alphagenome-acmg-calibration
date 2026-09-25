# Aim 7 (amendment A5): computational evidence (PP3/BP4) for ClinVar missense VUS.
# VUS: CLNSIG == "Uncertain_significance", >= 1 star, MANE Select missense (same exclusions as the
# reference sets), disease genes. Scores: CADD v1.7, AVI, and AlphaMissense/REVEL/BayesDel extracted
# with the same rules as script 20. Thresholds: medians of (a) the main analysis (CADD, AVI) and
# (b) the Aim 6 intersection analysis (all five predictors).
# Output: C:/cadd/calibration/results/v17_avi/vus/
suppressPackageStartupMessages({ library(data.table); library(rtracklayer) })
cdir <- "C:/cadd/data/coding"; cd <- "C:/cadd/comparators"; out <- "C:/cadd/calibration/results/v17_avi/vus"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("chrom", "pos", "ref", "alt")
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")

## ClinVar VUS and disease genes (rules of script 13) ---------------------------------------------
cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 3, 4, 5, 8),
            col.names = c("chrom", "pos", "clinvar_id", "ref", "alt", "info"), colClasses = list(character = c(1, 3)))
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                    grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                    grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
dg <- cv[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo),
         unique(unlist(lapply(strsplit(geneinfo, "|", fixed = TRUE), function(z) sub(":.*", "", z))))]
vus <- cv[clnsig == "Uncertain_significance" & stars >= 1 & nchar(ref) == 1 & nchar(alt) == 1, .(chrom, pos, ref, alt, clinvar_id, stars)]
rm(cv); invisible(gc())
cat("ClinVar SNV VUS (>=1 star):", nrow(vus), "\n")

## Restrict to MANE missense in disease genes; attach AF, CADD, AVI --------------------------------
af <- readRDS(file.path(cdir, "gnomad_af_mane_coding.rds"))[, .(chrom, pos, ref, alt, AF_max)]
v <- rbindlist(lapply(list.files(cdir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f) {
  d <- readRDS(f)[, .(chrom, pos, ref, alt, gene_name, consequence)]
  d <- d[vus, on = k, nomatch = NULL]
  if (!nrow(d)) return(NULL)
  bad <- unique(d[consequence %chin% c("stop_gained", "stop_lost", "start_lost", "incomplete_terminal_codon_variant"), ..k])
  d <- d[consequence == "missense"][!bad, on = k]
  x <- d[, .(gene = paste(sort(unique(gene_name)), collapse = "&"), dgene = any(gene_name %chin% dg)), by = k][dgene == TRUE][, dgene := NULL]
  if (!nrow(x)) return(NULL)
  ch <- x$chrom[1]
  s1 <- readRDS(file.path(cdir, sprintf("cadd_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, CADD_PHRED = PHRED)]
  s2 <- readRDS(file.path(cdir, sprintf("avi_coding_%s.rds", ch)))[, .(chrom, pos, ref, alt, AVI_PHRED = round(AVI_PHRED, 2))]
  merge(merge(x, s1, by = k), s2, by = k)
}))
v[af, on = k, AF_max := i.AF_max][is.na(AF_max), AF_max := 0]
rm(af); invisible(gc())
cat("MANE missense VUS in disease genes:", nrow(v), "| AF < 0.01:", v[AF_max < 0.01, .N], "\n")

## Comparator scores (same rules as script 20) -----------------------------------------------------
keys <- v[, ..k]
am <- fread(file.path(cd, "AlphaMissense_hg38.tsv.gz"), skip = "#CHROM", select = c(1, 2, 3, 4, 9),
            col.names = c("chrom", "pos", "ref", "alt", "AM"), colClasses = list(character = 1))
am[, chrom := sub("^chr", "", chrom)]; am <- am[keys, on = k, nomatch = NULL][, .(AM = max(AM)), by = k]; invisible(gc())
rv <- fread(cmd = sprintf('C:\\Windows\\System32\\tar.exe -xOf "%s" revel_with_transcript_ids',
                          normalizePath(file.path(cd, "revel-v1.3_all_chromosomes.zip"), winslash = "\\")),
            select = c("chr", "grch38_pos", "ref", "alt", "REVEL"), colClasses = list(character = c("chr", "grch38_pos")))
rv <- rv[grch38_pos != "."][, .(chrom = chr, pos = as.integer(grch38_pos), ref, alt, REVEL)]
rv <- rv[keys, on = k, nomatch = NULL][, .(REVEL = max(REVEL)), by = k]; invisible(gc())
chain <- import.chain(file.path(cd, "bayesdel", "hg38ToHg19.over.chain"))
lo <- liftOver(GRanges(paste0("chr", keys$chrom), IRanges(keys$pos, width = 1), strand = "+"), chain)
ok <- which(lengths(lo) == 1); l1 <- unlist(lo[ok]); comp <- c(A = "T", C = "G", G = "C", T = "A")
bd <- data.table(keys[ok], h19_chr = as.character(seqnames(l1)), h19_pos = start(l1), h19_strand = as.character(strand(l1)))
bd[, `:=`(ref19 = ifelse(h19_strand == "-", comp[ref], ref), alt19 = ifelse(h19_strand == "-", comp[alt], alt))]
tb <- TwoBitFile(file.path(cd, "bayesdel", "hg19.2bit")); bd <- bd[h19_chr %in% seqnames(seqinfo(tb))]
bd[, base19 := as.character(getSeq(tb, GRanges(h19_chr, IRanges(h19_pos, width = 1))))]
cat("VUS liftOver: unmapped", sum(lengths(lo) != 1), "| reference mismatch", bd[base19 != ref19, .N], "\n")
bd <- bd[base19 == ref19][, BayesDel := NA_real_]
for (b in c("A", "C", "G", "T")) {
  idx <- which(bd$alt19 == b); if (!length(idx)) next
  q <- GRanges(bd$h19_chr[idx], IRanges(bd$h19_pos[idx], width = 1))
  x <- import(BigWigFile(file.path(cd, "bayesdel", sprintf("NoAFBayesDelMut%s.bw", b))), which = reduce(q))
  ov <- findOverlaps(q, x); bd$BayesDel[idx[queryHits(ov)]] <- score(x)[subjectHits(ov)]
}
v <- Reduce(function(a, b) merge(a, b, by = k, all.x = TRUE), list(v, am, rv, bd[, c(k, "BayesDel"), with = FALSE]))
v[, `:=`(REVEL = round(REVEL, 3), AM = round(AM, 4), BayesDel = round(BayesDel, 2))]
saveRDS(v, file.path(out, "vus_scores.rds"))

## Evidence levels and points ---------------------------------------------------------------------
main <- fread("C:/cadd/calibration/results/v17_avi/summary/thresholds_all_runs.tsv")[analysis == "main",
          .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)]
inter <- fread("C:/cadd/calibration/results/v17_avi/comparators/summary/thresholds_all_runs.tsv")[,
          .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)]
pts <- c(PP3_VeryStrong = 8, PP3_Strong = 4, PP3_Moderate = 2, PP3_Supporting = 1, Indeterminate = 0,
         BP4_Supporting = -1, BP4_Moderate = -2, BP4_Strong = -4, BP4_VeryStrong = -8)
assign_level <- function(s, m) {
  lv <- rep("Indeterminate", length(s))
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { t <- m[level == L, PP3]; if (length(t) && !is.na(t)) lv[!is.na(s) & s >= t] <- paste0("PP3_", L) }
  for (L in c("Supporting", "Moderate", "Strong", "VeryStrong")) { t <- m[level == L, BP4]; if (length(t) && !is.na(t)) lv[!is.na(s) & s <= t] <- paste0("BP4_", L) }
  lv[is.na(s)] <- NA; lv
}
cols <- c(CADD = "CADD_PHRED", AVI = "AVI_PHRED", REVEL = "REVEL", AlphaMissense = "AM", BayesDel = "BayesDel")
res <- list()
for (set in c("main", "intersection")) {
  thr <- if (set == "main") main else inter
  tl <- if (set == "main") c("CADD", "AVI") else names(cols)
  for (tn in tl) {
    lv <- assign_level(v[[cols[[tn]]]], thr[tool == tn])
    v[, (paste0("ev_", set, "_", tn)) := lv]
    for (sub in c("all", "AF<0.01")) {
      keep <- !is.na(lv) & (sub == "all" | v$AF_max < 0.01)
      t <- as.data.table(table(factor(lv[keep], levels = names(pts))))
      setnames(t, c("evidence", "N"))
      res[[paste(set, tn, sub)]] <- t[, `:=`(thresholds = set, tool = tn, subset = sub, pct = round(100 * N / sum(N), 1), points = pts[evidence])]
    }
  }
}
tab <- rbindlist(res)
cat("\n== VUS per evidence level (thresholds: main = CADD/AVI main analysis; intersection = Aim 6) ==\n")
print(dcast(tab[subset == "all"], thresholds + evidence ~ tool, value.var = "pct")[order(thresholds, -match(evidence, names(pts)))])
fwrite(tab, file.path(out, "table_vus_evidence.tsv"), sep = "\t")
xt <- v[, .N, by = .(CADD = ev_main_CADD, AVI = ev_main_AVI)]
cat("\n== CADD vs AVI evidence for VUS (main thresholds) ==\n"); print(dcast(xt, CADD ~ AVI, value.var = "N", fill = 0))
fwrite(xt, file.path(out, "table_vus_cadd_vs_avi.tsv"), sep = "\t")
