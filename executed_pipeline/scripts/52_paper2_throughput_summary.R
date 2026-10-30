# Paper 2, technical throughput test: distribution of the AVI attributions in 2,000 unlabeled gnomAD SNVs,
# used only to fix the "no predicted effect" cut-off and to check the stratum rule before freezing the plan.
# Also checks that the API AVI raw score matches the local Atlas download. Reads no ClinVar label.
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
d <- "C:/cadd/calibration/results/paper2/throughput_test"
x <- fread(file.path(d, "atlas_sample.tsv.gz")); s <- fread(file.path(d, "variants_sample.tsv"))
x <- merge(x, s, by = "variant")
fi <- grep("^fi_", names(x), value = TRUE)
grp <- list(splicing = "fi_MERGED_SPLICING",
            protein = c("fi_ALPHAMISSENSE", "fi_PROTEIN_TERMINATION", "fi_START_LOST", "fi_STOP_LOST"),
            regulatory = grep("^fi_MAX_ABS_", fi, value = TRUE),
            conservation = c("fi_CACTUS_241_WAY", "fi_PHASTCONS_470_WAY"))
stopifnot(length(unlist(grp)) + 2 == length(fi))   # plus IS_INSERTION / IS_DELETION (zero for SNVs)
pos <- function(cols) rowSums(pmax(as.matrix(x[, cols, with = FALSE]), 0))
P <- sapply(grp, pos); x[, pos_total := rowSums(P)]
x[, dominant := colnames(P)[max.col(P, ties.method = "first")]][, share := apply(P, 1, max) / pmax(pos_total, 1e-12)]
x[, sum_check := abs(rowSums(as.matrix(x[, fi, with = FALSE])) - avi_raw)]
cat("Attributions sum to AVI raw: max abs difference", signif(max(x$sum_check), 3), "\n")
cat("\nPositive attribution total, quantiles by region\n")
print(x[, as.list(round(quantile(pos_total, c(.05, .1, .25, .5, .75, .9)), 3)), by = region])
cat("\nDominant group (share >= 0.5) at candidate cut-offs, % by region\n")
for (cut in c(0.05, 0.1, 0.2)) {
  x[, stratum := fifelse(pos_total < cut, "no_effect", fifelse(share >= 0.5, dominant, "mixed"))]
  cat("cut-off", cut, "\n"); print(dcast(x[, .N, by = .(region, stratum)][, pct := round(100 * N / sum(N), 1), by = region], region ~ stratum, value.var = "pct", fill = 0))
}
# API versus local download (AVI raw)
tf <- TabixFile("C:/cadd/avi/alphagenome_variant_impact_score_snvs.tsv.gz"); open(tf); pre <- if (any(grepl("^chr", seqnamesTabix(tf)))) "chr" else ""
k <- x[sample(.N, 200)]; k[, c("chrom", "pos", "ref", "alt") := tstrsplit(sub("^chr", "", variant), "[:>]")]
loc <- rbindlist(lapply(seq_len(nrow(k)), function(i) { r <- scanTabix(tf, param = GRanges(paste0(pre, k$chrom[i]), IRanges(as.integer(k$pos[i]), width = 1)))[[1]]
  if (!length(r)) return(NULL); f <- tstrsplit(r, "\t"); data.table(variant = sprintf("chr%s:%s:%s>%s", sub("^chr", "", f[[1]]), f[[2]], f[[3]], f[[4]]), local_raw = as.numeric(f[[5]])) }))
m <- merge(k[, .(variant, avi_raw)], loc, by = "variant")
cat("\nAPI vs local AVI raw:", nrow(m), "variants matched; max abs difference", signif(max(abs(m$avi_raw - m$local_raw)), 3), "\n")
fwrite(read.delim(paste0(file.path(d, "atlas_sample.tsv.gz"), ".timing.tsv")), file.path(d, "throughput.tsv"), sep = "\t")
