# Aim 1 (amendment A8): calibration by gene class (constraint, inheritance, mechanism) for REVEL,
# BayesDel, AlphaMissense, CADD v1.7 and AVI on the 2026 main intersection sets.
# Output: C:/cadd/calibration/results/v17_avi/gene_classes/
suppressPackageStartupMessages(library(data.table))
source("C:/cadd/calibration/R/pejaver_port.R"); source("C:/cadd/calibration/R/pejaver_port_fast.R")
workers <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(workers)) workers <- 10L
out <- "C:/cadd/calibration/results/v17_avi/gene_classes"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
gdir <- "C:/cadd/annot/genes"; k <- c("chrom", "pos", "ref", "alt")
f <- function(p) list.files(gdir, p, full.names = TRUE)[1]

## Gene annotations --------------------------------------------------------------------------------
genes <- unique(readRDS("C:/cadd/data/coding/mane_cds_positions.rds")[, .(gene = gene_name, gene_id = sub("\\..*", "", gene_id))])
con <- fread(f("constraint_metrics"))[mane_select == TRUE, .(gene_id = sub("\\..*", "", gene_id), decile = lof.oe_ci.upper_bin_decile)]
con <- con[!is.na(decile), .(decile = decile[1]), by = gene_id]
genes <- merge(genes, con, by = "gene_id", all.x = TRUE)
genes[, constraint := fcase(decile <= 2, "LOEUF_0-2", decile <= 6, "LOEUF_3-6", decile <= 9, "LOEUF_7-9")]

ok_cls <- c("Definitive", "Strong", "Moderate")
cg <- fread(f("clingen_gene_validity"), skip = "GENE SYMBOL", fill = TRUE)
setnames(cg, make.names(names(cg)))
cg <- cg[CLASSIFICATION %chin% ok_cls & MOI %chin% c("AD", "AR", "XL", "SD", "MT", "UD"), .(moi_cg = paste(sort(unique(MOI)), collapse = "+")), by = .(gene = GENE.SYMBOL)]
gc <- fread(f("gencc_submissions"), select = c("gene_symbol", "moi_title", "classification_title"))
gc <- gc[classification_title %chin% ok_cls]
gc[, moi := fcase(moi_title == "Autosomal dominant", "AD", moi_title == "Autosomal recessive", "AR", default = "other")]
gc <- gc[, .(moi_gc = paste(sort(unique(moi)), collapse = "+")), by = .(gene = gene_symbol)]
inh <- merge(cg, gc, by = "gene", all = TRUE)
inh[, moi := fifelse(!is.na(moi_cg), moi_cg, moi_gc)]
inh[, inheritance := fcase(moi == "AD", "AD_only", moi == "AR", "AR_only", default = NA_character_)]
genes <- merge(genes, inh[, .(gene, inheritance, inh_source = fifelse(!is.na(moi_cg), "ClinGen", "GenCC"))], by = "gene", all.x = TRUE)

g2 <- fread(f("gene2phenotype"))
setnames(g2, make.names(names(g2)))
g2 <- g2[confidence %chin% c("definitive", "strong", "moderate")]
g2[, mech := fcase(molecular.mechanism == "loss of function", "LoF",
                   molecular.mechanism %chin% c("gain of function", "dominant negative"), "GoF_DN", default = "other")]
g2 <- g2[mech != "other", .(mechs = paste(sort(unique(mech)), collapse = "+")), by = .(gene = gene.symbol)]
g2[, mechanism := fcase(mechs == "LoF", "LoF", mechs == "GoF_DN", "GoF_DN", default = NA_character_)]
genes <- merge(genes, g2[, .(gene, mechanism)], by = "gene", all.x = TRUE)
fwrite(genes, file.path(out, "gene_classes.tsv"), sep = "\t")

## Reference sets (2026 main intersection, single-gene SNVs) ---------------------------------------
sc <- readRDS("C:/cadd/data/refsets/comparator_scores.rds")
sc[, `:=`(REVEL = round(REVEL, 3), AM = round(AM, 4), BayesDel = round(BayesDel, 2))]
tools <- list(REVEL = list(col = "REVEL", inc = 0.001), BayesDel = list(col = "BayesDel", inc = 0.01),
              AlphaMissense = list(col = "AM", inc = 0.001), CADD = list(col = "CADD_PHRED", inc = 0.1),
              AVI = list(col = "AVI_PHRED", inc = 0.1))
cols <- vapply(tools, `[[`, "", "col")
clin <- merge(readRDS("C:/cadd/data/refsets/clinvar_ref.rds")[, c(k, "gene", "y", "CADD_PHRED", "AVI_PHRED"), with = FALSE], sc, by = k)
gno  <- merge(readRDS("C:/cadd/data/refsets/gnomad_ref.rds")[, c(k, "gene", "CADD_PHRED", "AVI_PHRED"), with = FALSE], sc, by = k)
clin <- clin[complete.cases(clin[, ..cols]) & !grepl("&", gene, fixed = TRUE)]
gno  <- gno[complete.cases(gno[, ..cols]) & !grepl("&", gene, fixed = TRUE)]
ann <- genes[, .(gene, constraint, inheritance, mechanism)]
clin <- merge(clin, ann, by = "gene", all.x = TRUE); gno <- merge(gno, ann, by = "gene", all.x = TRUE)
cat(sprintf("single-gene intersection: ClinVar %d P/LP, %d B/LB; gnomAD %d; reference genes without a gene_name match in constraint table: %d\n",
            clin[y == 1, .N], clin[y == 0, .N], nrow(gno), uniqueN(clin[is.na(constraint), gene])))

classes <- list(constraint = c("LOEUF_0-2", "LOEUF_3-6", "LOEUF_7-9"), inheritance = c("AD_only", "AR_only"), mechanism = c("LoF", "GoF_DN"))
sizes <- rbindlist(lapply(names(classes), function(ax) rbindlist(lapply(classes[[ax]], function(cl)
  data.table(axis = ax, class = cl, genes = uniqueN(clin[get(ax) == cl, gene]), P_LP = clin[get(ax) == cl & y == 1, .N],
             B_LB = clin[get(ax) == cl & y == 0, .N], gnomAD = gno[get(ax) == cl, .N])))))
sizes[, calibrated := P_LP >= 1000 & B_LB >= 1000]
print(sizes); fwrite(sizes, file.path(out, "class_sizes.tsv"), sep = "\t")

## Calibration per class -----------------------------------------------------------------------------
for (i in which(sizes$calibrated)) for (tn in names(tools)) {
  ax <- sizes$axis[i]; cl <- sizes$class[i]
  o <- file.path(out, "runs", ax, cl, tn)
  if (file.exists(file.path(o, "thresholds.tsv"))) { cat("skip", ax, cl, tn, "\n"); next }
  dir.create(o, recursive = TRUE, showWarnings = FALSE)
  cc <- clin[get(ax) == cl]; gg <- gno[get(ax) == cl]
  t0 <- Sys.time()
  r <- calibrate_fast(cc[[tools[[tn]]$col]], cc$y, gg[[tools[[tn]]$col]], B = 10000, workers = workers,
                      seed = 20260925, increment = tools[[tn]]$inc)
  fwrite(r$thresholds, file.path(o, "thresholds.tsv"), sep = "\t")
  fwrite(r$curves, file.path(o, "posterior_curves.tsv"), sep = "\t")
  fwrite(r$boot_thresholds, file.path(o, "bootstrap_thresholds.tsv"), sep = "\t")
  cat(sprintf("[%s] %s %s %s: %.1f min\n", format(Sys.time(), "%H:%M"), ax, cl, tn, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## Summaries: thresholds with bootstrap intervals; overlap between classes; reclassification --------
lev <- c("Supporting", "Moderate", "Strong", "VeryStrong")
res <- rbindlist(lapply(which(sizes$calibrated), function(i) rbindlist(lapply(names(tools), function(tn) {
  o <- file.path(out, "runs", sizes$axis[i], sizes$class[i], tn)
  th <- fread(file.path(o, "thresholds.tsv"), colClasses = c(level = "character"))
  bt <- fread(file.path(o, "bootstrap_thresholds.tsv"))
  rbindlist(lapply(c("PP3", "BP4"), function(cr) data.table(axis = sizes$axis[i], class = sizes$class[i], tool = tn, crit = cr, level = lev,
    threshold = th[[paste0(cr, "_discounted")]][match(lev, th$level)],
    lo = vapply(lev, function(L) { v <- as.numeric(bt[[paste0(cr, "_", L)]]); if (all(is.na(v))) NA_real_ else quantile(v, 0.025, na.rm = TRUE) }, 0),
    hi = vapply(lev, function(L) { v <- as.numeric(bt[[paste0(cr, "_", L)]]); if (all(is.na(v))) NA_real_ else quantile(v, 0.975, na.rm = TRUE) }, 0))))
}))))
fwrite(res, file.path(out, "thresholds_by_class.tsv"), sep = "\t")
cmp <- res[!is.na(threshold) & level %in% c("Supporting", "Moderate", "Strong"), {
  pr <- if (.N >= 2) combn(.N, 2, simplify = FALSE) else list()
  rbindlist(lapply(pr, function(p) data.table(class_a = class[p[1]], class_b = class[p[2]], thr_a = threshold[p[1]], thr_b = threshold[p[2]],
    intervals_overlap = !(hi[p[1]] < lo[p[2]] | hi[p[2]] < lo[p[1]]))))
}, by = .(axis, tool, crit, level)]
cat("\n== Class comparisons (different = bootstrap 95% intervals do not overlap) ==\n")
print(cmp[order(axis, tool, crit, level)])
fwrite(cmp, file.path(out, "class_comparisons.tsv"), sep = "\t")

glob <- fread("C:/cadd/calibration/results/v17_avi/comparators/summary/thresholds_all_runs.tsv")[,
          .(PP3 = median(PP3_discounted), BP4 = median(BP4_discounted)), by = .(tool, level)]
assign_level <- function(s, pp3, bp4) {
  lv <- rep("Indeterminate", length(s))
  for (L in lev) { t <- pp3[L]; if (!is.na(t)) lv[s >= t] <- paste0("PP3_", L) }
  for (L in lev) { t <- bp4[L]; if (!is.na(t)) lv[s <= t] <- paste0("BP4_", L) }
  lv
}
recl <- rbindlist(lapply(which(sizes$calibrated), function(i) rbindlist(lapply(names(tools), function(tn) {
  ax <- sizes$axis[i]; cl <- sizes$class[i]; cc <- clin[get(ax) == cl]; s <- cc[[tools[[tn]]$col]]
  rc <- res[axis == ax & class == cl & tool == tn]
  a <- assign_level(s, setNames(rc[crit == "PP3", threshold], lev), setNames(rc[crit == "BP4", threshold], lev))
  gt <- glob[tool == tn]; b <- assign_level(s, setNames(gt$PP3, gt$level)[lev], setNames(gt$BP4, gt$level)[lev])
  data.table(axis = ax, class = cl, tool = tn, n = length(s), pct_changed = 100 * mean(a != b),
             pct_changed_PLP = 100 * mean((a != b)[cc$y == 1]), pct_changed_BLB = 100 * mean((a != b)[cc$y == 0]))
}))))
cat("\n== % of class reference variants whose evidence level changes (class-specific vs global thresholds) ==\n")
print(recl, digits = 3)
fwrite(recl, file.path(out, "reclassification.tsv"), sep = "\t")
