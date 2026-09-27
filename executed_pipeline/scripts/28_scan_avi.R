# Genome-wide scan of the AlphaGenome AVI SNV file (amendment A9): every SNV with AVI Phred >= THRESH,
# where the Phred score is rounded to 2 decimals exactly as in the calibration (plan v1). Same
# windowed, resumable tabix scan as 02_rank_top_windows.R; sequence names are converted from
# "chr1" to "1" so that the output matches the CADD naming and chromosome order.
# Usage: Rscript 28_scan_avi.R [THRESH] [WORKERS]
suppressPackageStartupMessages({ library(data.table); library(parallel) })
args    <- commandArgs(trailingOnly = TRUE)
thresh  <- if (length(args) >= 1) as.numeric(args[1]) else 25.18
workers <- if (length(args) >= 2) as.integer(args[2]) else 8L
infile  <- "C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz"
window  <- 2e6
outdir  <- "C:/cadd/avi"
partdir <- file.path(outdir, sprintf("parts_avi_phred%s", thresh)); dir.create(partdir, showWarnings = FALSE)
outfile <- file.path(outdir, sprintf("AVI_PHRED_ge%s_GRCh38.tsv.gz", thresh))
cols <- c("Chrom", "Pos", "Ref", "Alt", "AVI_raw", "AVI_PHRED")
cls  <- c("character", "integer", "character", "character", "numeric", "numeric")

seqs  <- Rsamtools::seqnamesTabix(Rsamtools::TabixFile(infile))
tasks <- CJ(chrom = seqs, start = seq(1, 250e6, by = window), sorted = FALSE)
tasks[, part := file.path(partdir, sprintf("%s_%09d.rds", chrom, as.integer(start)))]
todo  <- tasks[!file.exists(part)]
cat(sprintf("Sequences: %s\nWindows: %d (pending %d), workers: %d, rounded AVI Phred >= %s\n",
            paste(seqs, collapse = ","), nrow(tasks), nrow(todo), workers, thresh))

scan_window <- function(i) {
  t <- todo[i]
  rng <- GenomicRanges::GRanges(t$chrom, IRanges::IRanges(t$start, t$start + window - 1))
  lines <- Rsamtools::scanTabix(tbx, param = rng)[[1]]
  res <- if (length(lines)) {
    d <- data.table::fread(text = lines, sep = "\t", header = FALSE, col.names = cols, colClasses = cls, nThread = 1)
    d <- d[round(AVI_PHRED, 2) >= thresh]
    d[, Chrom := sub("^chr", "", Chrom)]
    d
  } else data.table::data.table()
  saveRDS(list(n = length(lines), d = res), t$part)
  length(lines)
}
if (nrow(todo)) {
  t0 <- Sys.time()
  cl <- makeCluster(workers)
  clusterExport(cl, c("todo", "window", "cols", "cls", "thresh", "infile"))
  clusterEvalQ(cl, { suppressPackageStartupMessages({ library(data.table); library(Rsamtools) }); tbx <- open(TabixFile(infile)); NULL })
  batches <- split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / (workers * 4)))
  for (b in seq_along(batches)) {
    parLapplyLB(cl, batches[[b]], scan_window)
    cat(sprintf("[%s] batch %d/%d, %.1f min elapsed\n", format(Sys.time(), "%H:%M:%S"), b, length(batches),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  stopCluster(cl)
}
parts   <- lapply(tasks$part, readRDS)
scanned <- sum(vapply(parts, `[[`, numeric(1), "n"))
dt      <- rbindlist(lapply(parts, `[[`, "d"))
cat(sprintf("Rows scanned: %s | SNVs with rounded AVI Phred >= %s: %s\n", format(scanned, big.mark = ","), thresh, format(nrow(dt), big.mark = ",")))
chr_levels <- c(as.character(1:22), "X", "Y", "MT")
dt[, chr_ord := match(Chrom, chr_levels, nomatch = length(chr_levels) + 1L)]
setorder(dt, -AVI_PHRED, -AVI_raw, chr_ord, Pos, Ref, Alt)
dt[, Rank := .I]
print(dt[, .N, by = .(bin = cut(AVI_PHRED, c(-Inf, 26.44, 30, 35, 40, 50, Inf), right = FALSE))][order(bin)])
fwrite(dt[, .(Rank, Chrom, Pos, Ref, Alt, AVI_raw, AVI_PHRED)], outfile, sep = "\t")
cat("Written:", outfile, "\n")
