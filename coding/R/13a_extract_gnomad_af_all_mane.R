# gnomAD v4.1 allele frequencies (exomes, genomes; overall, and per genetic ancestry group) for SNVs
# at every MANE Select coding position, from the Ensembl VEP 116 indexed cache. Same parsing as
# 13_build_reference_sets.R (step 3), extended to all MANE genes and to per-population columns, so
# that the table can be reused (temporal hold-out, Aim 4 catalogue, Aim 8 ancestry analysis).
# Output: C:/cadd/data/coding/gnomad_af_mane_coding.rds
# Usage: Rscript 13a_extract_gnomad_af_all_mane.R [workers]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 6L
cdir <- "C:/cadd/data/coding"
vep_dir <- "C:/cadd/vep/homo_sapiens/116_GRCh38"
cols <- strsplit(sub("^variation_cols\t", "", grep("^variation_cols", readLines(file.path(vep_dir, "info.txt")), value = TRUE)), ",")[[1]]
afcols <- grep("^gnomAD[eg]", cols, value = TRUE)          # gnomADe, gnomADe_AFR, ..., gnomADg_SAS
cat("frequency columns:", paste(afcols, collapse = ", "), "\n")

pos <- readRDS(file.path(cdir, "mane_cds_positions.rds"))
reg <- unique(pos[, .(chrom, pos)]); rm(pos); invisible(gc())
setorder(reg, chrom, pos)
reg[, grp := cumsum(c(1L, diff(pos) != 1L) | c(TRUE, chrom[-1] != chrom[-.N]))]
reg <- reg[, .(start = min(pos), end = max(pos)), by = .(chrom, grp)]
reg[, batch := paste(chrom, (seq_len(.N) - 1L) %/% 1000L), by = chrom]
cat("regions:", nrow(reg), "| batches:", uniqueN(reg$batch), "\n")

cl <- makeCluster(workers)
clusterExport(cl, c("vep_dir", "cols", "afcols"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })))
t0 <- Sys.time()
af <- rbindlist(parLapplyLB(cl, split(reg, reg$batch), function(b) {
  f <- file.path(vep_dir, b$chrom[1], "all_vars.gz")
  l <- unlist(scanTabix(TabixFile(f, index = paste0(f, ".csi")),
                        param = GenomicRanges::GRanges(b$chrom, IRanges::IRanges(b$start, b$end))), use.names = FALSE)
  if (!length(l)) return(NULL)
  sel <- which(cols %in% c("chr", "start", "allele_string", afcols))
  d <- fread(text = l, sep = "\t", header = FALSE, colClasses = "character", nThread = 1, select = sel, col.names = cols[sel])
  d <- d[grepl("^[ACGT](/[ACGT])+$", allele_string)]
  d <- d[Reduce(`|`, lapply(afcols, function(k) d[[k]] != "."))]
  if (!nrow(d)) return(NULL)
  long <- rbindlist(lapply(afcols, function(k) {
    x <- d[get(k) != ".", .(chr, pos = as.integer(start), allele_string, v = get(k))]
    if (!nrow(x)) return(NULL)
    x <- x[, .(e = unlist(strsplit(v, ",", fixed = TRUE))), by = .(chr, pos, allele_string)]
    x[, .(chrom = chr, pos, ref = sub("/.*", "", allele_string), alt = sub(":.*", "", e),
          field = k, af = as.numeric(sub(".*:", "", e)))]
  }))
  long <- long[nchar(alt) == 1]
  dcast(long, chrom + pos + ref + alt ~ field, value.var = "af", fun.aggregate = max, fill = NA)
}), fill = TRUE)
stopCluster(cl)
for (k in afcols) if (k %in% names(af)) set(af, which(is.infinite(af[[k]])), k, NA_real_)
af[, AF_max := pmax(gnomADe, gnomADg, na.rm = TRUE)]
setnames(af, c("gnomADe", "gnomADg"), c("AF_exomes", "AF_genomes"))
cat(sprintf("SNV alleles with gnomAD frequency fields: %s in %.1f min\n", format(nrow(af), big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("duplicated alleles:", sum(duplicated(af, by = c("chrom", "pos", "ref", "alt"))), "\n")
saveRDS(af, file.path(cdir, "gnomad_af_mane_coding.rds"))
