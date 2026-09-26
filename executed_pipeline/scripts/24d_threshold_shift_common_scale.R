# POST HOC (not pre-specified in amendment A6; added 2026-09-26 to compare predictors on a common
# scale). Threshold shifts between sets are expressed as the percentage of the 2026 gnomAD
# intersection variants (1,888,558) that meet each PP3 (score >= t) or BP4 (score <= t) threshold,
# which removes the differences in score scales between predictors.
# Output: release_series/table_shift_common_scale.tsv
suppressPackageStartupMessages(library(data.table))
rs <- "C:/cadd/calibration/results/v17_avi/release_series"; k <- c("chrom", "pos", "ref", "alt")
sc <- readRDS("C:/cadd/data/refsets/comparator_scores_all_releases.rds")
sc[, `:=`(REVEL = round(REVEL, 3), AM = round(AM, 4), BayesDel = round(BayesDel, 2))]
cols <- c(REVEL = "REVEL", BayesDel = "BayesDel", AlphaMissense = "AM", CADD = "CADD_PHRED", AVI = "AVI_PHRED")
g <- merge(readRDS("C:/cadd/data/refsets/gnomad_ref.rds")[, c(k, "CADD_PHRED", "AVI_PHRED"), with = FALSE], sc, by = k)
g <- g[complete.cases(g[, ..cols])]
th <- fread(file.path(rs, "thresholds_all.tsv"))
th[, pct_gnomAD_PP3 := mapply(function(tl, t) if (is.na(t)) NA_real_ else 100 * mean(g[[cols[[tl]]]] >= t), tool, PP3_discounted)]
th[, pct_gnomAD_BP4 := mapply(function(tl, t) if (is.na(t)) NA_real_ else 100 * mean(g[[cols[[tl]]]] <= t), tool, BP4_discounted)]
w <- function(v) dcast(th[level %in% c("Supporting", "Moderate", "Strong")], tool + level ~ set, value.var = v)
p <- w("pct_gnomAD_PP3"); b <- w("pct_gnomAD_BP4")
p[, `:=`(post_minus_pre = post2022 - pre2022, r2026_minus_r2019 = r2026 - r2019)]
b[, `:=`(post_minus_pre = post2022 - pre2022, r2026_minus_r2019 = r2026 - r2019)]
cat("gnomAD intersection variants:", nrow(g), "\n")
cat("\n== % of gnomAD variants meeting PP3 (higher = less stringent) ==\n"); print(p, digits = 3)
cat("\n== % of gnomAD variants meeting BP4 (higher = less stringent) ==\n"); print(b, digits = 3)
fwrite(rbind(p[, criterion := "PP3"], b[, criterion := "BP4"]), file.path(rs, "table_shift_common_scale.tsv"), sep = "\t")
