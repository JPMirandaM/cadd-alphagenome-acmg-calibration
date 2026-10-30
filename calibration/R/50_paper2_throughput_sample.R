# Paper 2, technical throughput test (allowed before freezing analysis-plan-paper2-v1; reads no ClinVar label).
# Draws 1,000 unlabeled rare gnomAD v4.1 MANE missense SNVs (window-sizing set of the main article) and 1,000
# unlabeled rare gnomAD SNVs of the intronic splice region (A10 set), seed 20261030, for Atlas API queries.
# Usage: Rscript 50_paper2_throughput_sample.R
suppressPackageStartupMessages(library(data.table))
set.seed(20261030)
mis <- as.data.table(readRDS("C:/cadd/data/refsets/gnomad_ref.rds"))
spl <- as.data.table(readRDS("C:/cadd/calibration/results/v17_avi/noncoding/splice_region_intronic/gnomad_set.rds"))
stopifnot(all(is.na(mis$clnsig)))   # window-sizing set: no ClinVar record
s <- rbind(mis[sample(.N, 1000), .(chrom, pos, ref, alt, region = "missense")],
           spl[sample(.N, 1000), .(chrom, pos, ref, alt, region = "splice_region_intronic")])
s[, variant := sprintf("chr%s:%d:%s>%s", chrom, pos, ref, alt)]
out <- "C:/cadd/calibration/results/paper2/throughput_test/variants_sample.tsv"
fwrite(s[, .(variant, region)], out, sep = "\t")
cat("Written", nrow(s), "variants to", out, "\n"); print(s[, .N, by = region])
