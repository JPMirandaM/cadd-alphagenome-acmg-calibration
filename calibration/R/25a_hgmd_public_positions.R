# Amendment A7 (a): HGMD Public (release 2020.4, as shipped in the Ensembl VEP 116 cache) records.
# Every all_vars.gz file is streamed in chunks; records with allele_string "HGMD_MUTATION" are kept.
# Output: C:/cadd/data/refsets/hgmd_public_2020.4.rds (chrom, start, end, accession) + counts.
suppressPackageStartupMessages({ library(data.table); library(parallel) })
vep_dir <- "C:/cadd/vep/homo_sapiens/116_GRCh38"
cols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", readLines(file.path(vep_dir, "info.txt")), value = TRUE)), ",")[[1]]
i_name <- which(cols == "variation_name"); i_start <- which(cols == "start"); i_end <- which(cols == "end"); i_al <- which(cols == "allele_string")
cat("HGMD source in info.txt:", grep("HGMD", readLines(file.path(vep_dir, "info.txt")), value = TRUE), "\n")
chroms <- c(as.character(1:22), "X", "Y", "MT")
cl <- makeCluster(6)
clusterExport(cl, c("vep_dir", "i_name", "i_start", "i_end", "i_al"))
invisible(clusterEvalQ(cl, library(data.table)))
t0 <- Sys.time()
res <- parLapplyLB(cl, chroms, function(ch) {
  f <- file.path(vep_dir, ch, "all_vars.gz")
  if (!file.exists(f)) return(NULL)
  con <- gzfile(f, "r"); on.exit(close(con))
  out <- list(); n_lines <- 0
  repeat {
    l <- readLines(con, n = 2e6)
    if (!length(l)) break
    n_lines <- n_lines + length(l)
    h <- l[grepl("\tHGMD_MUTATION\t", l, fixed = TRUE)]
    if (length(h)) {
      p <- tstrsplit(h, "\t", fixed = TRUE, keep = c(i_name, i_start, i_end, i_al))
      out[[length(out) + 1]] <- data.table(chrom = ch, accession = p[[1]], start = as.integer(p[[2]]), end = p[[3]], allele = p[[4]])
    }
  }
  d <- rbindlist(out)
  if (nrow(d)) d <- d[allele == "HGMD_MUTATION"][, allele := NULL]
  list(d = d, n_lines = n_lines)
})
stopCluster(cl)
h <- rbindlist(lapply(res, `[[`, "d"))
h[, end := fifelse(end == ".", start, as.integer(end))]
cat(sprintf("all_vars records scanned: %s in %.1f min\n", format(sum(vapply(res, function(x) if (is.null(x)) 0 else x$n_lines, 0)), big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("HGMD Public records:", nrow(h), "| unique accessions:", uniqueN(h$accession),
    "| unique positions (start):", uniqueN(h, by = c("chrom", "start")),
    "| single-nucleotide records (start = end):", h[start == end, .N], "\n")
h[, class := substr(accession, 1, 2)]
print(h[, .(records = .N, accessions = uniqueN(accession)), by = class][order(-records)])
saveRDS(h, "C:/cadd/data/refsets/hgmd_public_2020.4.rds")
