# Aim 3 (A13), step b: download scores and GRCh38 mappings of the included MaveDB score sets,
# assign the functional class (normal / abnormal) from the selected score calibration, and map
# the classified variants to MANE Select missense SNVs.
suppressMessages({ library(jsonlite); library(data.table); library(curl) })
`%||%` <- function(a, b) if (is.null(a)) b else a
api <- "https://api.mavedb.org/api/v1"
mdir <- "C:/cadd/data/mave"; cdir <- "C:/cadd/data/coding"
res_dir <- "C:/cadd/calibration/results/v17_avi/mave"
for (d in c("scores", "mapped", "classes")) dir.create(file.path(mdir, d), showWarnings = FALSE)

h <- new_handle(); handle_setopt(h, ssl_options = 2L, timeout = 1200L)
dl <- function(url, file) {
  if (!file.exists(file) || file.size(file) == 0) {
    for (i in 1:3) { ok <- tryCatch({ curl_download(url, file, handle = h, quiet = TRUE); TRUE }, error = function(e) FALSE)
      if (ok) break; Sys.sleep(10) }
    Sys.sleep(0.5)
  }
  file.exists(file) && file.size(file) > 0
}

inv <- fread(file.path(res_dir, "mavedb_human_scoreset_inventory.tsv"))
cand <- unique(inv[n_normal_ranges > 0 & n_abnormal_ranges > 0 & mapping_state != "failed", .(urn, genes, licence, mapping_state)])
cat("Candidate score sets (normal + abnormal classes, mapping not failed):", nrow(cand), "\n")

# gene symbol: first HGNC-like token of the target name (e.g. "TP53 (P72R)" -> TP53)
cand[, gene := sub("^([A-Za-z0-9-]+).*$", "\\1", genes)]

pick_calibration <- function(cals) {
  ok <- Filter(function(c) { cl <- vapply(c$functionalClassifications, function(f) f$functionalClassification %||% "NA", "")
                             any(cl == "normal") && any(cl == "abnormal") }, cals)
  if (!length(ok)) return(NULL)
  inv_p <- Filter(function(c) isTRUE(c$investigatorProvided) || grepl("investigator", c$title %||% "", ignore.case = TRUE), ok)
  if (length(inv_p)) return(list(cal = inv_p[[1]], rule = "investigator-provided"))
  nr <- Filter(function(c) !isTRUE(c$researchUseOnly), ok)
  if (length(nr)) return(list(cal = nr[[1]], rule = "first non-research-use-only"))
  list(cal = ok[[1]], rule = "first (research-use-only)")
}

# Range-based classes are applied to the score; class-based (categorical) classes have no range and
# their member variants are downloaded from the calibration endpoint.
classify <- function(score, fcs, acc = NULL, cal_urn = NULL, key = NULL) {
  out <- rep(NA_character_, length(score))
  for (f in fcs) {
    cls <- f$functionalClassification %||% "NA"; if (!cls %in% c("normal", "abnormal")) next
    if (!length(f$range)) {
      mf <- file.path(mdir, "classes", sprintf("%s_%s.json", key, f$id))
      if (!dl(sprintf("%s/score-calibrations/%s/functional-classifications/%s/variants", api, cal_urn, f$id), mf)) stop("class download failed")
      mem <- vapply(fromJSON(mf, simplifyVector = FALSE)$variants, function(v) v$urn, "")
      inr <- acc %in% mem
      out[inr & is.na(out)] <- cls
      out[inr & !is.na(out) & out != cls] <- "conflict"
      next
    }
    lo <- f$range[[1]]; hi <- f$range[[2]]
    lo <- if (is.null(lo)) -Inf else as.numeric(lo); hi <- if (is.null(hi)) Inf else as.numeric(hi)
    inc_lo <- f$inclusiveLowerBound %||% TRUE; inc_hi <- f$inclusiveUpperBound %||% FALSE
    inr <- (if (inc_lo) score >= lo else score > lo) & (if (inc_hi) score <= hi else score < hi)
    inr[is.na(inr)] <- FALSE
    out[inr & is.na(out)] <- cls
    out[inr & !is.na(out) & out != cls] <- "conflict"
  }
  out
}

# GRCh38 RefSeq chromosome accessions
grch38 <- c(setNames(as.character(1:22), sprintf("NC_%06d.%s", 1:22,
  c(11, 12, 12, 12, 10, 12, 14, 11, 12, 11, 10, 12, 11, 9, 10, 10, 11, 10, 10, 11, 9, 11))),
  "NC_000023.11" = "X", "NC_000024.10" = "Y")

aa3 <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
         Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V")

rows <- list(); qc <- list()
for (i in seq_len(nrow(cand))) {
  u <- cand$urn[i]; key <- gsub(":", "_", u)
  s <- fromJSON(file.path(mdir, "scoreset_json", paste0(key, ".json")), simplifyVector = FALSE)
  pc <- pick_calibration(s$scoreCalibrations)
  if (is.null(pc)) { qc[[u]] <- data.table(urn = u, status = "no usable calibration"); next }
  if (!dl(paste0(api, "/score-sets/", u, "/scores"), file.path(mdir, "scores", paste0(key, ".csv"))) ||
      !dl(paste0(api, "/score-sets/", u, "/mapped-variants"), file.path(mdir, "mapped", paste0(key, ".json")))) {
    qc[[u]] <- data.table(urn = u, status = "download failed"); next }
  sc <- fread(file.path(mdir, "scores", paste0(key, ".csv")), select = c("accession", "score"), colClasses = list(character = "accession"))
  sc[, score := suppressWarnings(as.numeric(score))]
  sc[, class := classify(score, pc$cal$functionalClassifications, accession, pc$cal$urn, key)]
  mv <- fromJSON(file.path(mdir, "mapped", paste0(key, ".json")), simplifyVector = FALSE)
  mvt <- rbindlist(lapply(mv, function(m) {
    ex <- m$postMapped$expressions
    e <- if (length(ex)) vapply(ex, function(x) x$value %||% NA_character_, "") else NA_character_
    data.table(accession = m$variantUrn, current = isTRUE(m$current), date = m$mappedDate %||% "", level = m$alignmentLevel %||% NA_character_,
               expr = e)
  }))
  mvt <- mvt[!is.na(expr)]
  # current mapping if any, else the most recent mapping date
  mvt[, keep := if (any(current)) current else date == max(date), by = accession]
  mvt <- unique(mvt[keep == TRUE, .(accession, level, expr)])
  x <- merge(sc[class %in% c("normal", "abnormal", "conflict")], mvt, by = "accession")
  g <- x[grepl(":g\\.[0-9]+[ACGT]>[ACGT]$", expr)]
  p <- x[grepl(":p\\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}$", expr)]
  out <- NULL; mtype <- NA_character_
  if (nrow(g)) {
    mtype <- "genomic"
    g[, `:=`(acc = sub(":.*", "", expr), pos = as.integer(sub(".*:g\\.([0-9]+).*", "\\1", expr)),
             ref = sub(".*[0-9]([ACGT])>.*", "\\1", expr), alt = sub(".*>([ACGT])$", "\\1", expr))]
    g <- g[acc %in% names(grch38)][, chrom := grch38[acc]]
    out <- g[, .(urn = u, gene = cand$gene[i], map = "genomic", chrom, pos, ref, alt, hgvsp = NA_character_, class, score)]
  } else if (nrow(p)) {
    mtype <- "protein"
    p[, hgvsp := sub(".*:(p\\..*)$", "\\1", expr)]
    out <- p[, .(urn = u, gene = cand$gene[i], map = "protein", chrom = NA_character_, pos = NA_integer_, ref = NA_character_,
                 alt = NA_character_, hgvsp, class, score)]
  }
  qc[[u]] <- data.table(urn = u, gene = cand$gene[i], status = if (is.null(out)) "no SNV / missense mapping" else "ok",
                        calibration = pc$cal$title %||% NA_character_, calibration_rule = pc$rule, mapping = mtype,
                        n_scored = nrow(sc), n_normal = sum(sc$class == "normal", na.rm = TRUE),
                        n_abnormal = sum(sc$class == "abnormal", na.rm = TRUE), n_mapped_classified = if (is.null(out)) 0L else nrow(out))
  if (!is.null(out)) rows[[u]] <- out
  cat(sprintf("[%d/%d] %s %s: %s\n", i, nrow(cand), u, cand$gene[i], qc[[u]]$status))
}
lab <- rbindlist(rows); qc <- rbindlist(qc, fill = TRUE)

# Map to MANE Select missense SNVs
genes <- unique(lab$gene)
chroms <- c(1:22, "X", "Y")
mane <- rbindlist(lapply(chroms, function(ch) {
  x <- readRDS(file.path(cdir, sprintf("coding_snv_%s.rds", ch)))
  x[consequence == "missense" & gene_name %chin% genes, .(chrom, pos, ref, alt, gene = gene_name, hgvsp, CADD_PHRED = PHRED)]
}))
cat("MANE missense SNVs in MAVE genes:", nrow(mane), "\n")

gm <- merge(lab[map == "genomic"], mane[, .(chrom, pos, ref, alt, gene_mane = gene, hgvsp_mane = hgvsp)],
            by = c("chrom", "pos", "ref", "alt"))
gm <- gm[gene == gene_mane]
pm_in <- lab[map == "protein"]
pm <- merge(pm_in[, !c("chrom", "pos", "ref", "alt")], mane[, .(chrom, pos, ref, alt, gene, hgvsp_mane = hgvsp)],
            by.x = c("gene", "hgvsp"), by.y = c("gene", "hgvsp_mane"), allow.cartesian = TRUE)
# protein-level reference residue check: share of mapped (position, reference residue) pairs whose
# residue equals the MANE Select protein residue at that position (a DMS covers all 19 substitutions,
# only ~1/3 of which are reachable by one SNV, so matching whole changes would understate the match)
mane_res <- unique(mane[, .(gene, protein_pos = as.integer(sub("^p\\.[A-Za-z]{3}([0-9]+).*", "\\1", hgvsp)),
                            aa_mane = sub("^p\\.([A-Za-z]{3}).*", "\\1", hgvsp))])
pres <- unique(pm_in[, .(urn, gene, protein_pos = as.integer(sub("^p\\.[A-Za-z]{3}([0-9]+).*", "\\1", hgvsp)),
                         aa_mave = sub("^p\\.([A-Za-z]{3}).*", "\\1", hgvsp))])
pres <- merge(pres, mane_res, by = c("gene", "protein_pos"), all.x = TRUE)
pr <- pres[, .(n_positions = .N, match_rate = mean(!is.na(aa_mane) & aa_mane == aa_mave)), by = urn]
gr <- merge(lab[map == "genomic", .(n_snv = .N), by = urn], gm[, .(n_match = .N), by = urn], by = "urn", all.x = TRUE)
qc <- merge(qc, rbind(pr[, .(urn, match_rate)], gr[, .(urn, match_rate = NA_real_)]), by = "urn", all.x = TRUE)
drop_sets <- pr[match_rate < 0.9, urn]
cat("Protein-level score sets dropped (< 90% MANE match):", length(drop_sets), "\n")
qc[urn %chin% drop_sets, status := "dropped: < 90% MANE match"]

all <- rbind(gm[, .(urn, gene, chrom, pos, ref, alt, hgvsp = hgvsp_mane, class, score, map)],
             pm[!urn %chin% drop_sets, .(urn, gene, chrom, pos, ref, alt, hgvsp, class, score, map)])
all <- all[class %in% c("normal", "abnormal")]
# consensus per SNV across score sets
snv <- all[, .(n_sets = uniqueN(urn), classes = paste(sort(unique(class)), collapse = "/"), urns = paste(sort(unique(urn)), collapse = ";")),
           by = .(gene, chrom, pos, ref, alt, hgvsp)]
cat("Classified MANE missense SNVs:", nrow(snv), "| discordant across sets (excluded):", snv[classes == "abnormal/normal", .N], "\n")
snv <- snv[classes != "abnormal/normal"][, abnormal := as.integer(classes == "abnormal")][, classes := NULL]
saveRDS(snv, file.path(mdir, "mave_labels_mane_missense.rds"))
fwrite(qc, file.path(res_dir, "scoreset_qc.tsv"), sep = "\t")
per_gene <- snv[, .(n = .N, abnormal = sum(abnormal), normal = sum(abnormal == 0)), by = gene][order(-n)]
fwrite(per_gene, file.path(res_dir, "labels_per_gene.tsv"), sep = "\t")
print(qc[, .N, by = status]); print(per_gene, nrows = 100)
