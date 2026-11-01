# Paper 2, amendment A-P2-4: functional label sets (SGE, MPRA) and current VUS, with region classes, and the list of
# SNVs whose AVI attributions are still to be retrieved from the Atlas API. No AVI score or attribution is read here.
# 1. SGE: the 12 MaveDB score sets of the main article with a genomic mapping (status ok in scoreset_qc.tsv); classes
#    from the calibration selected in A13 (cached files of script 33b); intermediate, unspecified and conflicting
#    classes excluded; SNVs classified discordantly by two sets excluded; exonic_other and intergenic excluded.
# 2. MPRA (Kircher et al. 2019, GRCh38): functional = pValue < 1e-5 and |Coefficient| >= 0.1; all SNVs kept.
# 3. Current VUS: exact Uncertain_significance, >= 1 star, 2026-09-23 release, 2026 disease gene (MANE class gene),
#    AF < 0.01 (rules of script 59), class other than exonic_other or intergenic; expert-panel flag.
# Output: C:/cadd/calibration/results/paper2/functional/ (sge_labels.rds, mpra_labels.rds, vus_current.rds, counts) and
#         C:/cadd/calibration/results/paper2/universe/variants_functional.tsv (Atlas query list)
suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(Rsamtools) })
source("C:/cadd/calibration/R/paper2_regions.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
out <- "C:/cadd/calibration/results/paper2/functional"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
mdir <- "C:/cadd/data/mave"; k <- c("chrom", "pos", "ref", "alt")
log <- function(...) cat(sprintf(...), "\n")
m <- mane_structure(); edge <- exon_edge_positions()
coding_pos <- unique(readRDS(file.path(CDIR, "mane_cds_positions.rds"))[, .(chrom, pos)])
region <- function(x) {   # x: chrom, pos, ref, alt -> + gene, class (all classes, exonic_other and intergenic kept)
  x <- unique(x[, ..k])
  uc <- x[coding_pos, on = .(chrom, pos), nomatch = NULL]; un <- x[!coding_pos, on = .(chrom, pos)]
  cc <- rbindlist(lapply(split(uc, uc$chrom), function(z) classify_coding_chr(z, edge)))
  nc <- merge(un, classify_noncoding(un, m), by = c("chrom", "pos"))
  rbind(cc, nc[, .(chrom, pos, ref, alt, gene, class)])
}

## 1. SGE ------------------------------------------------------------------------------------------------
urns <- c("urn:mavedb:00000662-0-1", "urn:mavedb:00000097-0-2", "urn:mavedb:00001242-a-1", "urn:mavedb:00001226-a-1",
          "urn:mavedb:00001262-a-3", "urn:mavedb:00000658-0-1", "urn:mavedb:00001259-a-2", "urn:mavedb:00000673-0-1",
          "urn:mavedb:00001260-a-2", "urn:mavedb:00001265-a-2", "urn:mavedb:00000675-a-1", "urn:mavedb:00001264-a-2")
qc <- fread("C:/cadd/calibration/results/v17_avi/mave/scoreset_qc.tsv")[urn %chin% urns]
stopifnot(nrow(qc) == 12, all(qc$status == "ok"), all(qc$mapping == "genomic"))
pick_calibration <- function(cals) {   # as in script 33b
  ok <- Filter(function(c) { cl <- vapply(c$functionalClassifications, function(f) f$functionalClassification %||% "NA", "")
                             any(cl == "normal") && any(cl == "abnormal") }, cals)
  if (!length(ok)) return(NULL)
  inv_p <- Filter(function(c) isTRUE(c$investigatorProvided) || grepl("investigator", c$title %||% "", ignore.case = TRUE), ok)
  if (length(inv_p)) return(inv_p[[1]])
  nr <- Filter(function(c) !isTRUE(c$researchUseOnly), ok)
  if (length(nr)) return(nr[[1]]); ok[[1]]
}
classify <- function(score, fcs, acc, key) {   # as in script 33b, from cached class-member files only
  o <- rep(NA_character_, length(score))
  for (f in fcs) {
    cls <- f$functionalClassification %||% "NA"; if (!cls %in% c("normal", "abnormal")) next
    if (!length(f$range)) {
      mf <- file.path(mdir, "classes", sprintf("%s_%s.json", key, f$id)); stopifnot(file.exists(mf))
      inr <- acc %in% vapply(fromJSON(mf, simplifyVector = FALSE)$variants, function(v) v$urn, "")
    } else {
      lo <- f$range[[1]]; hi <- f$range[[2]]
      lo <- if (is.null(lo)) -Inf else as.numeric(lo); hi <- if (is.null(hi)) Inf else as.numeric(hi)
      inc_lo <- f$inclusiveLowerBound %||% TRUE; inc_hi <- f$inclusiveUpperBound %||% FALSE
      inr <- (if (inc_lo) score >= lo else score > lo) & (if (inc_hi) score <= hi else score < hi); inr[is.na(inr)] <- FALSE
    }
    o[inr & is.na(o)] <- cls; o[inr & !is.na(o) & o != cls] <- "conflict"
  }
  o
}
grch38 <- c(setNames(as.character(1:22), sprintf("NC_%06d.%s", 1:22,
  c(11, 12, 12, 12, 10, 12, 14, 11, 12, 11, 10, 12, 11, 9, 10, 10, 11, 10, 10, 11, 9, 11))), "NC_000023.11" = "X", "NC_000024.10" = "Y")
sge <- rbindlist(lapply(seq_len(nrow(qc)), function(i) {
  u <- qc$urn[i]; key <- gsub(":", "_", u)
  s <- fromJSON(file.path(mdir, "scoreset_json", paste0(key, ".json")), simplifyVector = FALSE)
  cal <- pick_calibration(s$scoreCalibrations)
  sc <- fread(file.path(mdir, "scores", paste0(key, ".csv")), select = c("accession", "hgvs_nt", "score"), colClasses = list(character = c("accession", "hgvs_nt")))
  sc[, score := suppressWarnings(as.numeric(score))][, fclass := classify(score, cal$functionalClassifications, accession, key)]
  mv <- fromJSON(file.path(mdir, "mapped", paste0(key, ".json")), simplifyVector = FALSE)
  mvt <- rbindlist(lapply(mv, function(z) { ex <- z$postMapped$expressions
    data.table(accession = z$variantUrn, current = isTRUE(z$current), date = z$mappedDate %||% "",
               expr = if (length(ex)) vapply(ex, function(x) x$value %||% NA_character_, "") else NA_character_) }))[!is.na(expr)]
  mvt[, keep := if (any(current)) current else date == max(date), by = accession]
  g <- merge(sc[fclass %chin% c("normal", "abnormal", "conflict")], unique(mvt[keep == TRUE, .(accession, expr)]), by = "accession")
  g <- g[grepl(":g\\.[0-9]+[ACGT]>[ACGT]$", expr)]
  g[, `:=`(acc = sub(":.*", "", expr), pos = as.integer(sub(".*:g\\.([0-9]+).*", "\\1", expr)),
           ref = sub(".*[0-9]([ACGT])>.*", "\\1", expr), alt = sub(".*>([ACGT])$", "\\1", expr))]
  g <- g[acc %in% names(grch38), .(accession, hgvs_nt, chrom = grch38[acc], pos, ref, alt, fclass, src = "mavedb")]
  # Intronic SNVs (c.N+k, c.N-k) have no MaveDB genomic mapping: anchored on the genomic position of exonic base c.N
  # from the MaveDB mapping of the same set, shifted by k along the transcript strand (deviation of 2026-11-01).
  ct <- function(h) sub("^.*c\\.([-*]?[0-9]+)([ACGT])>([ACGT])$", "\\1", h)
  an <- g[grepl("c\\.[-*]?[0-9]+[ACGT]>[ACGT]$", hgvs_nt), .(tok = ct(hgvs_nt), chrom, pos)]
  an <- an[, .(chrom = chrom[1], pos = pos[1], n = uniqueN(paste(chrom, pos))), by = tok][n == 1]
  cd <- an[grepl("^[0-9]+$", tok)][, cpos := as.integer(tok)]
  strand <- if (nrow(cd) >= 2) sign(cor(cd$cpos, cd$pos, method = "spearman")) else NA
  it <- sc[fclass %chin% c("normal", "abnormal", "conflict") & grepl("c\\.[-*]?[0-9]+[+-][0-9]+[ACGT]>[ACGT]$", hgvs_nt) & !accession %chin% g$accession]
  if (nrow(it) && !is.na(strand)) {
    pr <- regmatches(it$hgvs_nt, regexec("c\\.([-*]?[0-9]+)([+-][0-9]+)([ACGT])>([ACGT])$", it$hgvs_nt))
    it[, `:=`(tok = vapply(pr, `[`, "", 2), off = as.integer(vapply(pr, `[`, "", 3)), tr = vapply(pr, `[`, "", 4), ta = vapply(pr, `[`, "", 5))]
    it <- merge(it, an[, .(tok, chrom, apos = pos)], by = "tok")
    cmp <- c(A = "T", C = "G", G = "C", T = "A")
    it[, `:=`(pos = apos + strand * off, ref = if (strand > 0) tr else cmp[tr], alt = if (strand > 0) ta else cmp[ta])]
    g <- rbind(g, it[, .(accession, hgvs_nt, chrom, pos, ref, alt, fclass, src = "anchored_intronic")])
  }
  log("%s %s: %d MaveDB-mapped SNVs, %d intronic SNVs anchored (strand %s)", u, qc$gene[i], g[src == "mavedb", .N], g[src != "mavedb", .N], strand)
  g[, .(urn = u, set_gene = qc$gene[i], chrom, pos, ref, alt, fclass, src)]
}))
log("SGE classified genomic SNV records: %d (conflicting within a set: %d)", nrow(sge), sge[fclass == "conflict", .N])
# Reference check of the anchored intronic SNVs against GRCh38 (reference column of the CADD v1.7 SNV file)
tf <- TabixFile("C:/cadd/data/whole_genome_SNVs.tsv.gz"); open(tf)
ai <- unique(sge[src != "mavedb", .(chrom, pos)])
gref <- rbindlist(lapply(split(ai, ai$chrom), function(x) {
  l <- unlist(scanTabix(tf, param = GenomicRanges::GRanges(x$chrom[1], IRanges::IRanges(x$pos, x$pos))), use.names = FALSE)
  if (!length(l)) return(NULL); unique(fread(text = l, sep = "\t", header = FALSE, select = 1:3, col.names = c("chrom", "pos", "gref"), colClasses = list(character = 1)))
})); close(tf)
sge <- merge(sge, gref, by = c("chrom", "pos"), all.x = TRUE)
print(sge[src != "mavedb", .(anchored = .N, ref_match = sum(ref == gref, na.rm = TRUE)), by = set_gene])
sge <- sge[src == "mavedb" | (!is.na(gref) & ref == gref)][, gref := NULL]
sge <- sge[fclass != "conflict"]
sge <- sge[, .(set_gene = set_gene[1], urn = paste(unique(urn), collapse = ","), src = src[1], classes = paste(sort(unique(fclass)), collapse = "/")), by = k]
log("SGE unique SNVs: %d; discordant between sets (excluded): %d", nrow(sge), sge[grepl("/", classes), .N])
sge <- sge[!grepl("/", classes)][, y := as.integer(classes == "abnormal")][, classes := NULL]
sge <- merge(sge, region(sge), by = k)
log("SGE SNVs by region class (before exclusion):"); print(sge[, .N, by = class][order(-N)])
sge <- sge[!class %chin% c("exonic_other", "intergenic")]
log("SGE analysis set: %d SNVs (%d abnormal, %d normal)", nrow(sge), sge[y == 1, .N], sge[y == 0, .N])
fwrite(dcast(sge[, .N, by = .(set_gene, class, y)], set_gene + class ~ y, value.var = "N", fill = 0), file.path(out, "sge_counts.tsv"), sep = "\t")
saveRDS(sge, file.path(out, "sge_labels.rds"))

## 2. MPRA -----------------------------------------------------------------------------------------------
mp <- fread("C:/cadd/annot/mpra/elements.tsv.gz")[Release == "GRCh38" & Ref %chin% c("A", "C", "G", "T") & Alt %chin% c("A", "C", "G", "T")]
mp <- mp[, .(chrom = as.character(Chrom), pos = as.integer(Pos), ref = Ref, alt = Alt, element = Element,
             y = as.integer(pValue < 1e-5 & abs(Coefficient) >= 0.1))]
log("MPRA SNV records: %d", nrow(mp)); mp <- unique(mp, by = c(k, "element"))
mp <- mp[, .(element = paste(unique(element), collapse = ","), y = max(y), ny = uniqueN(y)), by = k]
log("MPRA unique SNVs: %d (discordant across elements: %d, kept as functional if any)", nrow(mp), mp[ny > 1, .N]); mp[, ny := NULL]
mp <- merge(mp, region(mp), by = k, all.x = TRUE)
log("MPRA: %d SNVs (%d functional, %d non-functional)", nrow(mp), mp[y == 1, .N], mp[y == 0, .N])
fwrite(dcast(mp[, .N, by = .(element, class, y)], element + class ~ y, value.var = "N", fill = 0), file.path(out, "mpra_counts.tsv"), sep = "\t")
saveRDS(mp, file.path(out, "mpra_labels.rds"))

## 3. Current VUS ------------------------------------------------------------------------------------------
getk <- function(x, kk) { v <- sub(paste0(".*(^|;)", kk, "=([^;]*).*"), "\\2", x); v[!grepl(paste0("(^|;)", kk, "="), x)] <- NA; v }
plp <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5, 8), col.names = c("chrom", "pos", "ref", "alt", "info"), colClasses = list(character = 1))
cv[, `:=`(clnsig = getk(info, "CLNSIG"), revstat = getk(info, "CLNREVSTAT"), geneinfo = getk(info, "GENEINFO"))][, info := NULL]
cv[, stars := fcase(grepl("practice_guideline", revstat), 4L, grepl("reviewed_by_expert_panel", revstat), 3L,
                    grepl("criteria_provided,_multiple_submitters,_no_conflicts", revstat), 2L,
                    grepl("criteria_provided,_single_submitter|criteria_provided,_conflicting", revstat), 1L, default = 0L)]
dg <- cv[clnsig %chin% plp & stars >= 1 & !is.na(geneinfo), unique(sub(":.*", "", unlist(strsplit(geneinfo, "|", fixed = TRUE))))]
vus <- unique(cv[grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt) & stars >= 1 & clnsig == "Uncertain_significance", .(chrom, pos, ref, alt, stars)], by = k)
rm(cv); invisible(gc())
log("2026 SNV VUS (>= 1 star): %d (expert panel or guideline: %d)", nrow(vus), vus[stars >= 3, .N])
vus <- merge(vus, region(vus), by = k)[!class %chin% c("exonic_other", "intergenic") & gene %chin% dg]
af <- unique(rbind(readRDS(file.path(CDIR, "gnomad_af_mane_coding.rds"))[, .(chrom, pos, ref, alt, AF_max)],
                   rbindlist(lapply(list.files("C:/cadd/calibration/results/paper2/refsets/tmp", "^af_.*\\.rds$", full.names = TRUE), readRDS))[, .(chrom, pos, ref, alt, AF_max)]))
af <- af[, .(AF_max = max(AF_max)), by = k]; vus <- merge(vus, af, by = k, all.x = TRUE); rm(af); invisible(gc())
vep <- "C:/cadd/vep/homo_sapiens/116_GRCh38"   # gnomAD AF from the VEP 116 cache, as in scripts 54, 59 and 61 (absent = 0)
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
  r[, chrom := x$chrom[1]]
}
miss <- vus[is.na(AF_max), ..k]
if (nrow(miss)) {
  cl <- parallel::makeCluster(6L)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) }))
  parallel::clusterExport(cl, c("vep", "vcols", "af_chr"))
  afm <- rbindlist(parallel::parLapplyLB(cl, split(miss, miss$chrom), af_chr)); parallel::stopCluster(cl)
  if (nrow(afm)) vus[afm, on = k, AF_max := i.AF_max]
}
log("Current VUS in a 2026 disease gene with a region class: %d; AF looked up in the VEP cache: %d; absent from gnomAD (AF = 0): %d",
    nrow(vus), nrow(miss), vus[is.na(AF_max), .N])
vus[is.na(AF_max), AF_max := 0]; vus <- vus[AF_max < 0.01]
log("Current rare VUS: %d (expert panel or guideline: %d)", nrow(vus), vus[stars >= 3, .N])
saveRDS(vus, file.path(out, "vus_current.rds"))

## 4. Atlas query list --------------------------------------------------------------------------------------
have <- rbindlist(lapply(list.files(c("C:/cadd/calibration/results/paper2/atlas_clinvar", "C:/cadd/calibration/results/paper2/atlas_clinvar_batch2",
                                      "C:/cadd/calibration/results/paper2/atlas_gnomad"), "^chunk_.*\\.tsv\\.gz$", full.names = TRUE),
                         function(f) fread(f, sep = "\t", select = c("variant", "error"), colClasses = "character")))[error == "" | is.na(error), unique(variant)]
q <- unique(rbind(sge[, ..k], mp[, ..k], vus[, ..k]))[, variant := sprintf("chr%s:%d:%s>%s", chrom, pos, ref, alt)]
log("SNVs needed: %d; already retrieved: %d; to query: %d (of which current VUS: %d)", nrow(q), sum(q$variant %chin% have),
    sum(!q$variant %chin% have), sum(!vus[, sprintf("chr%s:%d:%s>%s", chrom, pos, ref, alt)] %chin% have))
fwrite(q[!variant %chin% have, .(variant)], "C:/cadd/calibration/results/paper2/universe/variants_functional.tsv", sep = "\t")
