# Paper 2, technical throughput test: alternative stratum rule (mechanistic groups only; conservation as a
# separate "conservation only" stratum) evaluated on the same 2,000 unlabeled gnomAD SNVs. Reads no ClinVar label.
suppressPackageStartupMessages(library(data.table))
d <- "C:/cadd/calibration/results/paper2/throughput_test"
x <- merge(fread(file.path(d, "atlas_sample.tsv.gz")), fread(file.path(d, "variants_sample.tsv")), by = "variant")
fi <- grep("^fi_", names(x), value = TRUE)
grp <- list(splicing = "fi_MERGED_SPLICING", protein = c("fi_ALPHAMISSENSE","fi_PROTEIN_TERMINATION","fi_START_LOST","fi_STOP_LOST"), regulatory = grep("^fi_MAX_ABS_", fi, value = TRUE))
P <- sapply(grp, function(cols) rowSums(pmax(as.matrix(x[, cols, with = FALSE]), 0)))
cons <- pmax(x$fi_CACTUS_241_WAY, 0) + pmax(x$fi_PHASTCONS_470_WAY, 0)
mech <- rowSums(P); dom <- colnames(P)[max.col(P, ties.method = "first")]; sh <- apply(P, 1, max) / pmax(mech, 1e-12)
for (cut in c(0.05, 0.1)) {
  st <- fifelse(mech >= cut, fifelse(sh >= 0.5, dom, "mixed"), fifelse(cons >= cut, "conservation_only", "no_effect"))
  cat("mechanistic cut-off", cut, "\n"); print(dcast(data.table(region = x$region, st)[, .N, by = .(region, st)][, pct := round(100*N/sum(N),1), by = region], region ~ st, value.var = "pct", fill = 0))
}
