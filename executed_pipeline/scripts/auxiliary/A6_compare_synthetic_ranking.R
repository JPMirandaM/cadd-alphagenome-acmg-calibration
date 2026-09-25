# Validation of 02_rank_top_windows.R on the synthetic set built by A1_make_synthetic_test_set.R:
# compares its outputs with a full in-memory sort using the same ordering criteria.
# Usage: Rscript A6_compare_synthetic_ranking.R <test_dir>
suppressPackageStartupMessages(library(data.table))
s <- commandArgs(trailingOnly = TRUE)[1]
top  <- fread(file.path(s, "top20000_CADDv1.7_GRCh38_PHRED.tsv.gz"))
ref  <- fread(file.path(s, "expected.tsv"))                       # top 20,000 from in-memory sort
all  <- fread(file.path(s, "PHRED_ge20_CADDv1.7_GRCh38.tsv.gz"))
src  <- fread(file.path(s, "t.tsv"), skip = 2, header = FALSE)     # full synthetic set
cat("Synthetic SNVs:", nrow(src), "\n")
cat("Top 20,000 identical to in-memory sort:", isTRUE(all.equal(top, ref)), "\n")
cat("All PHRED >= 20: written", nrow(all), "| expected", src[V6 >= 20, .N], "| equal:",
    nrow(all) == src[V6 >= 20, .N], "\n")
cat("First 20,000 of the threshold output identical to top 20,000:", isTRUE(all.equal(all[1:20000], ref)), "\n")
cat("Forced ties (PHRED 99, RawScore 9.5) resolved by position:\n")
print(top[PHRED == 99])
