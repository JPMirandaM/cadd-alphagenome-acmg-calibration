# Extracts AlphaGenome AVI scores (raw and Phred-scaled) for every possible SNV at the MANE Select
# coding positions (same regions as 09_extract_cadd_coding.R), by tabix queries on the local AVI
# file (sequence names carry the "chr" prefix; output uses CADD-style names without it).
# Output: avi_coding_<chrom>.rds  (chrom, pos, ref, alt, AVI_raw, AVI_PHRED)
# Usage: Rscript 12_extract_avi_coding.R [workers]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 10L
avi    <- "C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz"
outdir <- "C:/cadd/data/coding"

pos <- readRDS(file.path(outdir, "mane_cds_positions.rds"))
u <- unique(pos[, .(chrom, pos)])
setorder(u, chrom, pos)
u[, grp := cumsum(c(1L, diff(pos) != 1L) | c(TRUE, chrom[-1] != chrom[-.N]))]
reg <- u[, .(start = min(pos), end = max(pos)), by = .(chrom, grp)]
cat("coding positions:", nrow(u), "| merged regions:", nrow(reg), "\n")
reg[, batch := paste(chrom, (seq_len(.N) - 1L) %/% 2000L), by = chrom]
batches <- split(reg, reg$batch)

cl <- makeCluster(workers)
clusterExport(cl, "avi")
clusterEvalQ(cl, { suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
                   tbx <- open(TabixFile(avi)); NULL })
t0 <- Sys.time()
res <- parLapplyLB(cl, batches, function(b) {
  gr <- GenomicRanges::GRanges(paste0("chr", b$chrom), IRanges::IRanges(b$start, b$end))
  l <- unlist(scanTabix(tbx, param = gr), use.names = FALSE)
  if (!length(l)) return(NULL)
  d <- fread(text = l, sep = "\t", header = FALSE, nThread = 1,
             col.names = c("chrom", "pos", "ref", "alt", "AVI_raw", "AVI_PHRED"),
             colClasses = c("character", "integer", "character", "character", "numeric", "numeric"))
  d[, chrom := sub("^chr", "", chrom)]
  d
})
stopCluster(cl)
dt <- rbindlist(res)
cat(sprintf("rows: %s in %.1f min\n", format(nrow(dt), big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

per_pos <- dt[, .N, by = .(chrom, pos)]
cat("positions returned:", nrow(per_pos), "of", nrow(u), "| with exactly 3 alleles:", per_pos[N == 3, .N], "\n")
cat("duplicated rows:", sum(duplicated(dt, by = c("chrom", "pos", "alt"))), "\n")
missing <- fsetdiff(u[, .(chrom, pos)], per_pos[, .(chrom, pos)])
cat("coding positions without AVI rows:", nrow(missing), "\n")
# Reference-allele agreement with CADD (both GRCh38)
ref_mismatch <- 0L
for (ch in unique(dt$chrom)) {
  cadd_ref <- unique(readRDS(file.path(outdir, sprintf("cadd_coding_%s.rds", ch)))[, .(pos, ref_cadd = ref)])
  ref_mismatch <- ref_mismatch + merge(unique(dt[chrom == ch, .(pos, ref)]), cadd_ref, by = "pos")[ref != ref_cadd, .N]
  saveRDS(dt[chrom == ch], file.path(outdir, sprintf("avi_coding_%s.rds", ch)))
}
cat("positions whose reference allele differs between AVI and CADD:", ref_mismatch, "\n")
fwrite(missing, file.path(outdir, "coding_positions_without_avi.tsv"), sep = "\t")
