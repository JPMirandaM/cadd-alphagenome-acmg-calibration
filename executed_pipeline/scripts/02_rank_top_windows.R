# Select SNVs with CADD PHRED >= THRESH (v1.7, GRCh38) on Windows using Rsamtools (tabix) + data.table,
# ordered by PHRED (descending), RawScore (descending) and genomic position. With N > 0 only the
# top N are kept.
# The file has ~9 billion rows: it is scanned in parallel over genomic windows, keeping only rows
# with PHRED >= THRESH. Because ALL rows above the threshold are retained, the ranking is exact
# (with N > 0, provided there are >= N candidates; this is checked).
# Each window is saved to parts/, so an interrupted run resumes where it stopped.
#
# Usage: Rscript rank_top_windows.R [file.tsv.gz] [N (0 = all)] [THRESH] [WORKERS]

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

args    <- commandArgs(trailingOnly = TRUE)
infile  <- normalizePath(if (length(args) >= 1) args[1] else "C:/cadd/data/whole_genome_SNVs.tsv.gz")
topn    <- if (length(args) >= 2) as.integer(args[2]) else 0L       # 0 = all variants >= THRESH
thresh  <- if (length(args) >= 3) as.numeric(args[3]) else 25      # PHRED 25 ~ top 0.32%
workers <- if (length(args) >= 4) as.integer(args[4]) else 10L
window  <- 2e6                                                      # ~6 M rows per query

outdir  <- dirname(infile)
partdir <- file.path(outdir, sprintf("parts_phred%s", thresh))
dir.create(partdir, showWarnings = FALSE)
outfile <- file.path(outdir, if (topn > 0) sprintf("top%d_CADDv1.7_GRCh38_PHRED.tsv.gz", topn)
                             else sprintf("PHRED_ge%s_CADDv1.7_GRCh38.tsv.gz", thresh))

cols <- c("Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED")
cls  <- c("character", "integer", "character", "character", "numeric", "numeric")

# Non-overlapping windows [s, s+window-1] up to 250 Mb (covers the longest chromosome);
# empty windows are resolved instantly through the index.
seqs  <- Rsamtools::seqnamesTabix(Rsamtools::TabixFile(infile))
tasks <- CJ(chrom = seqs, start = seq(1, 250e6, by = window), sorted = FALSE)
tasks[, part := file.path(partdir, sprintf("%s_%09d.rds", chrom, as.integer(start)))]
todo  <- tasks[!file.exists(part)]
cat(sprintf("Chromosomes: %s\nWindows: %d (pending %d), workers: %d, PHRED >= %s\n",
            paste(seqs, collapse = ","), nrow(tasks), nrow(todo), workers, thresh))

scan_window <- function(i) {
  t   <- todo[i]
  rng <- GenomicRanges::GRanges(t$chrom, IRanges::IRanges(t$start, t$start + window - 1))
  lines <- Rsamtools::scanTabix(tbx, param = rng)[[1]]
  res <- if (length(lines)) {
    d <- data.table::fread(text = lines, sep = "\t", header = FALSE,
                           col.names = cols, colClasses = cls, nThread = 1)
    d[PHRED >= thresh]
  } else data.table::data.table()
  saveRDS(list(n = length(lines), d = res), t$part)
  length(lines)
}

if (nrow(todo)) {
  t0 <- Sys.time()
  cl <- makeCluster(workers)
  clusterExport(cl, c("todo", "window", "cols", "cls", "thresh", "infile"))
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
    tbx <- open(TabixFile(infile))
    NULL
  })
  # Batches, to report progress
  batches <- split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / (workers * 4)))
  for (b in seq_along(batches)) {
    parLapplyLB(cl, batches[[b]], scan_window)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] batch %d/%d, %.1f min elapsed\n",
                format(Sys.time(), "%H:%M:%S"), b, length(batches), el))
  }
  stopCluster(cl)
}

parts   <- lapply(tasks$part, readRDS)
scanned <- sum(vapply(parts, `[[`, numeric(1), "n"))
dt      <- rbindlist(lapply(parts, `[[`, "d"))
cat(sprintf("Rows scanned: %s | candidates with PHRED >= %s: %s\n",
            format(scanned, big.mark = ","), thresh, format(nrow(dt), big.mark = ",")))

if (nrow(dt) < topn)
  stop(sprintf("Only %d candidates; lower THRESH and rerun (a different parts folder is used).", nrow(dt)))

# Order: PHRED desc, RawScore desc, then genomic position (natural chromosome order, Pos asc).
# Ref/Alt last so the order is deterministic among alleles at the same site.
chr_levels <- c(as.character(1:22), "X", "Y", "MT")
dt[, chr_ord := match(Chrom, chr_levels, nomatch = length(chr_levels) + 1L)]
setorder(dt, -PHRED, -RawScore, chr_ord, Pos, Ref, Alt)

top <- if (topn > 0) dt[seq_len(topn)] else dt
top[, Rank := .I]

if (topn > 0) {
  last       <- top[topn]
  n_tied_out <- dt[(topn + 1L):.N][PHRED == last$PHRED & RawScore == last$RawScore, .N]
  cat(sprintf("Row %d: %s:%d %s>%s PHRED = %s, RawScore = %s\n", topn,
              last$Chrom, last$Pos, last$Ref, last$Alt, last$PHRED, last$RawScore))
  cat("Exact ties (PHRED and RawScore) excluded by position:", n_tied_out, "\n")
}
cat(sprintf("Variants written: %s | PHRED range [%s, %s]\n",
            format(nrow(top), big.mark = ","), min(top$PHRED), max(top$PHRED)))
print(top[, .N, by = .(PHRED_bin = cut(PHRED, c(-Inf, 25, 30, 35, 40, 45, 50, Inf), right = FALSE))
          ][order(PHRED_bin)])

fwrite(top[, .(Rank, Chrom, Pos, Ref, Alt, RawScore, PHRED)], outfile, sep = "\t")
cat("Written:", outfile, "\n")
