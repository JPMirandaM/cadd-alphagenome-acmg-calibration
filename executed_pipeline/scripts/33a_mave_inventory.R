# Aim 3 (MAVE validation), step a: inventory of MaveDB score sets (metadata only; no scores).
# Lists every published human score set (target taxonomy 9606 or a human transcript/genomic
# accession) with its investigator-provided functional classifications ("normal" / "abnormal"),
# licence and variant-mapping state. Used to select score sets before any score is looked at.
suppressMessages({ library(jsonlite); library(data.table); library(curl) })
`%||%` <- function(a, b) if (is.null(a)) b else a
api <- "https://api.mavedb.org/api/v1"
out <- "C:/cadd/data/mave"
res_dir <- "C:/cadd/calibration/results/v17_avi/mave"
dir.create(file.path(out, "scoreset_json"), recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

h <- new_handle(ssl_verifypeer = TRUE)
handle_setopt(h, ssl_options = 2L)  # CURLSSLOPT_NO_REVOKE (corporate network)
get_json <- function(url, file) {
  if (!file.exists(file)) {
    for (i in 1:3) {
      ok <- tryCatch({ curl_download(url, file, handle = h, quiet = TRUE); TRUE }, error = function(e) FALSE)
      if (ok) break; Sys.sleep(5)
    }
    Sys.sleep(0.3)
  }
  if (!file.exists(file)) return(NULL)
  fromJSON(file, simplifyVector = FALSE)
}

lst <- fromJSON("C:/cadd/data/mave/scoreset_search_2026-09-27.json", simplifyVector = FALSE)
is_human <- vapply(lst, function(s) any(vapply(s$targetGenes, function(t) {
  code <- t$targetSequence$taxonomy$code
  if (!is.null(code)) return(code == 9606)
  acc <- t$targetAccession$accession
  !is.null(acc) && grepl("^(ENST|NM_|NP_|NC_0000|ENSP|NG_)", acc)
}, logical(1))), logical(1))
urns <- vapply(lst[is_human], `[[`, "", "urn")
cat("Published score sets:", length(lst), "| human targets:", length(urns), "\n")

rows <- list()
for (u in urns) {
  s <- get_json(paste0(api, "/score-sets/", u), file.path(out, "scoreset_json", paste0(gsub(":", "_", u), ".json")))
  if (is.null(s)) { rows[[u]] <- data.table(urn = u, fetch_ok = FALSE); next }
  tg <- s$targetGenes
  genes <- paste(unique(vapply(tg, function(t) { g <- t$targetAccession$gene; if (is.null(g)) t$name else g }, "")), collapse = ",")
  acc <- paste(unique(unlist(lapply(tg, function(t) t$targetAccession$accession))), collapse = ",")
  ttype <- paste(unique(vapply(tg, function(t) if (!is.null(t$targetAccession)) "accession" else if (!is.null(t$targetSequence)) t$targetSequence$sequenceType else "NA", "")), collapse = ",")
  cal <- s$scoreCalibrations
  base <- data.table(urn = u, fetch_ok = TRUE, title = s$title, genes = genes, target_accession = acc,
                     target_type = ttype, num_variants = s$numVariants %||% NA_integer_,
                     licence = s$license$shortName %||% NA_character_,
                     mapping_state = s$mappingState %||% NA_character_,
                     superseded_by = if (is.null(s$supersedingScoreSet)) NA_character_ else s$supersedingScoreSet$urn,
                     meta_analysis = length(s$metaAnalyzesScoreSetUrns) > 0,
                     pmid = paste(unlist(lapply(s$primaryPublicationIdentifiers, `[[`, "identifier")), collapse = ","),
                     published = s$publishedDate %||% NA_character_)
  if (!length(cal)) { rows[[u]] <- cbind(base, calibration = NA_character_, research_use_only = NA, investigator_provided = NA,
                                           n_normal_ranges = 0L, n_abnormal_ranges = 0L, normal_ranges = NA_character_, abnormal_ranges = NA_character_); next }
  rows[[u]] <- rbindlist(lapply(cal, function(c) {
    fc <- c$functionalClassifications
    cls <- vapply(fc, function(f) f$functionalClassification %||% "NA", "")
    rng <- vapply(fc, function(f) paste0("[", paste(vapply(f$range, function(x) if (is.null(x)) "NA" else format(x), ""), collapse = ","), "]"), "")
    cbind(base, data.table(calibration = c$title %||% NA_character_, research_use_only = isTRUE(c$researchUseOnly),
                           investigator_provided = isTRUE(c$investigatorProvided) || grepl("investigator", c$title %||% "", ignore.case = TRUE),
                           n_normal_ranges = sum(cls == "normal"), n_abnormal_ranges = sum(cls == "abnormal"),
                           normal_ranges = paste(rng[cls == "normal"], collapse = ";"),
                           abnormal_ranges = paste(rng[cls == "abnormal"], collapse = ";")))
  }))
}
inv <- rbindlist(rows, fill = TRUE)
fwrite(inv, file.path(res_dir, "mavedb_human_scoreset_inventory.tsv"), sep = "\t")
cat("Human score sets fetched:", uniqueN(inv[fetch_ok == TRUE]$urn), "\n")
cat("With >= 1 calibration having normal and abnormal classes:",
    uniqueN(inv[n_normal_ranges > 0 & n_abnormal_ranges > 0]$urn), "\n")
print(inv[n_normal_ranges > 0 & n_abnormal_ranges > 0, .(urn, genes, num_variants, target_type, licence, mapping_state, calibration, research_use_only)][order(genes)], nrows = 300)
