# Extracts CADD v1.7 scores for every possible SNV at the MANE Select coding positions built by
# 08_build_mane_cds.R, by tabix queries on the local CADD file (merged CDS regions, parallel).
# Output: one file per chromosome, C:/cadd/data/coding/cadd_coding_<chrom>.rds
#   columns: chrom, pos, ref, alt, RawScore, PHRED   (3 rows per coding position)
# Usage: Rscript 09_extract_cadd_coding.R [workers]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 10L
cadd   <- "C:/cadd/data/whole_genome_SNVs.tsv.gz"
outdir <- "C:/cadd/data/coding"

pos <- readRDS(file.path(outdir, "mane_cds_positions.rds"))
u <- unique(pos[, .(chrom, pos)])
setorder(u, chrom, pos)
# Merge consecutive positions into regions
u[, grp := cumsum(c(1L, diff(pos) != 1L) | c(TRUE, chrom[-1] != chrom[-.N])), ]
reg <- u[, .(start = min(pos), end = max(pos), n = .N), by = .(chrom, grp)]
cat("coding positions:", nrow(u), "| merged regions:", nrow(reg), "\n")

reg[, batch := paste(chrom, (seq_len(.N) - 1L) %/% 2000L), by = chrom]
batches <- split(reg, reg$batch)

cl <- makeCluster(workers)
clusterExport(cl, "cadd")
clusterEvalQ(cl, { suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
                   tbx <- open(TabixFile(cadd)); NULL })
t0 <- Sys.time()
res <- parLapplyLB(cl, batches, function(b) {
  gr <- GenomicRanges::GRanges(b$chrom, IRanges::IRanges(b$start, b$end))
  l <- unlist(scanTabix(tbx, param = gr), use.names = FALSE)
  if (!length(l)) return(NULL)
  fread(text = l, sep = "\t", header = FALSE, nThread = 1,
        col.names = c("chrom", "pos", "ref", "alt", "RawScore", "PHRED"),
        colClasses = c("character", "integer", "character", "character", "numeric", "numeric"))
})
stopCluster(cl)
dt <- rbindlist(res)
cat(sprintf("rows: %s in %.1f min\n", format(nrow(dt), big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Checks: 3 alternative alleles per coding position, no missing positions, no duplicates
per_pos <- dt[, .N, by = .(chrom, pos)]
cat("positions returned:", nrow(per_pos), "of", nrow(u), "| with exactly 3 alleles:", per_pos[N == 3, .N], "\n")
cat("duplicated rows:", sum(duplicated(dt, by = c("chrom", "pos", "alt"))), "\n")
missing <- fsetdiff(u[, .(chrom, pos)], per_pos[, .(chrom, pos)])
cat("coding positions without CADD rows:", nrow(missing), "\n")
if (nrow(missing)) print(head(missing, 20))

for (ch in unique(dt$chrom)) saveRDS(dt[chrom == ch], file.path(outdir, sprintf("cadd_coding_%s.rds", ch)))
fwrite(missing, file.path(outdir, "coding_positions_without_cadd.tsv"), sep = "\t")
cat("Written per-chromosome files to", outdir, "\n")
