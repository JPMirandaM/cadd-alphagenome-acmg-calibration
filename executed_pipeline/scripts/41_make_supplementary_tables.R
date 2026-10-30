# Supplementary Tables S1-S20 of the manuscript (v9, 2026-10-29: renumbered in order of first citation; v7 added S19-S20 for amendment A16; v6 removed the MCPS/BRAVO table), assembled from the analysis outputs.
# Each table is written as TSV to figures/data/SupplementaryTables/ and all of them, with a
# "Contents" sheet giving titles and legends, to figures/Supplementary_Tables_S1-S20.xlsx.
# No new computation: values are copied (and rounded for display) from calibration/results/.
# Usage: Rscript 41_make_supplementary_tables.R
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
R <- "C:/cadd/calibration/results/v17_avi"; F <- "C:/cadd/figures"
out_dir <- file.path(F, "data", "SupplementaryTables"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(...) fread(file.path(R, ...))
num_round <- function(d, digits = 4) { d <- copy(d); for (j in names(d)) if (is.double(d[[j]])) set(d, j = j, value = signif(d[[j]], digits)); d }
tabs <- list(); meta <- list()
add <- function(id, title, legend, d) { tabs[[id]] <<- num_round(as.data.table(d)); meta[[id]] <<- data.table(Table = id, Title = title, Legend = legend) }

## S1-S2 CADD v1.6 reproduction and seed sensitivity ----------------------------------------------
add("S1", "Published versus reproduced CADD v1.6 thresholds",
    "Thresholds of Pejaver et al. (2022, Table 2) and those obtained with our R implementation on the authors' deposited data (10,000 bootstrap replicates, seed 20260925).",
    fread(file.path(F, "TableS1_port_validation.tsv")))
add("S2", "Seed sensitivity of the CADD v1.6 reproduction",
    "Thresholds for five random seeds (10,000 replicates each) with the fraction of bootstrap replicates that did not reach each level.",
    fread(file.path(F, "TableS2_seed_sensitivity.tsv")))

## S3 predictors, precision and increments -------------------------------------------------------
add("S3", "Predictors, score precision and window increments",
    "Score sources and versions, the precision at which scores were used and the window increment of the local posterior method, chosen to match the original study.",
    data.table(Predictor = c("CADD v1.7", "AVI (AlphaGenome Variant Impact)", "REVEL", "BayesDel (no allele frequency)", "AlphaMissense", "CADD v1.6 (Aim 2b)"),
               Source = c("Official GRCh38 whole_genome_SNVs.tsv.gz (MD5 verified)", "AlphaGenome Atlas SNV file (GRCh38)", "REVEL v1.3 official file, GRCh38 positions",
                          "UCSC hg19 bigWig tracks, lifted to GRCh38 with reference-base check", "Zenodo 10813168 hg38 predictions", "Official GRCh38 file, byte ranges covering the reference positions"),
               Score_used = c("PHRED", "Phred", "REVEL score", "BayesDel score", "am_pathogenicity", "PHRED"),
               Precision = c("1 decimal (as distributed)", "2 decimals (rounded)", "3 decimals (as distributed)", "2 decimals (rounded)", "4 decimals (as distributed)", "1 decimal (as distributed)"),
               Window_increment = c(0.1, 0.1, 0.001, 0.01, 0.001, 0.1),
               Multiple_values = c("one value per SNV", "one value per SNV", "maximum over transcripts", "one value per SNV and allele", "maximum over transcripts", "one value per SNV")))

## S4 sensitivity analyses -----------------------------------------------------------------------
p <- rd("summary", "table_sensitivity_PP3.tsv")[, criterion := "PP3"]; b <- rd("summary", "table_sensitivity_BP4.tsv")[, criterion := "BP4"]
add("S4", "Sensitivity analyses of the CADD v1.7 and AVI calibration",
    "Thresholds (seed 20260925) in the main analysis and in the four prespecified sensitivity analyses, namely review status of at least two stars (two_stars), rarity defined as AF below 0.001 (AF_below_0.001), gene-level bootstrap (gene_bootstrap) and exclusion of missense SNVs within 3 nucleotides of an internal exon boundary (no_exon_edge_3nt). Empty cells indicate levels not reached.",
    setnames(setcolorder(rbind(p, b, fill = TRUE), c("criterion", "tool", "level")), c("S1", "S2", "S3", "S4"),
             c("two_stars", "AF_below_0.001", "gene_bootstrap", "no_exon_edge_3nt")))   # descriptive names (v6): "S1-S4" collided with supplementary item numbers

## S5 HGMD exclusion ------------------------------------------------------------------------------
h <- rd("hgmd_excluded", "table_thresholds_hgmd.tsv"); hs <- rd("hgmd_excluded", "set_sizes.tsv")
add("S5", "Thresholds with and without HGMD Public positions",
    paste0("Thresholds of five predictors calibrated on the ClinVar 2019-12-02 and 2026-09-23 intersection sets with all variants (all) and after excluding SNVs at positions of single-nucleotide HGMD Public 2020.4 records (noHGMD), with the published values where available. Set sizes after exclusion were ",
           paste(sprintf("%s %d P/LP, %d B/LB, %d gnomAD", hs$set, hs$P_LP, hs$B_LB, hs$gnomAD), collapse = "; "), "."), h)

## S6 temporal hold-out ---------------------------------------------------------------------------
add("S6", "Temporal hold-out validation",
    "Positive likelihood ratios of each evidence interval of CADD v1.7 and AVI calibrated on the ClinVar release of 27 January 2024 and evaluated on 58,920 variants (22,626 P/LP, 36,294 B/LB) first classified afterwards, with 95% bootstrap intervals (1,000 replicates), the required ratio and whether it was met at the point estimate and at the conservative 95% bound. The second block compares 2024 and 2026 thresholds.",
    rbind(rd("temporal_holdout", "interval_LR_test_set.tsv")[, block := "interval likelihood ratios"],
          rd("temporal_holdout", "thresholds_2024_vs_2026.tsv")[, block := "thresholds 2024 vs 2026"], fill = TRUE))

## S7 comparators -----------------------------------------------------------------------------------
th <- rd("comparators", "summary", "table_thresholds.tsv"); di <- rd("comparators", "summary", "table_discrimination.tsv"); isz <- rd("comparators", "summary", "intersection_sizes.tsv")
add("S7", "Calibration and discrimination of five predictors on the intersection set",
    sprintf("Thresholds (median [range] over five seeds) and AUROC and AUPRC with 95%% bootstrap intervals and paired differences versus AVI (1,000 replicates), on the %d P/LP and %d B/LB ClinVar variants (%d genes) and %d gnomAD variants scored by all five predictors.",
            isz$P_LP, isz$B_LB, isz$genes, isz$gnomAD),
    rbind(th[, block := "thresholds"], di[, block := "discrimination"], fill = TRUE))

## S8 release series and drift ---------------------------------------------------------------------
pp <- rd("release_series", "table_PP3_by_set.tsv")[, criterion := "PP3"]; bb <- rd("release_series", "table_BP4_by_set.tsv")[, criterion := "BP4"]
sc <- rd("release_series", "table_shift_common_scale.tsv"); ss <- rd("release_series", "set_sizes.tsv")
add("S8", "Thresholds across ClinVar releases and before or after November 2022",
    paste0("Thresholds (seed 20260925) of five predictors calibrated on the intersection sets of ClinVar releases 2019-12-02, 2022-11-19, 2024-01-27 and 2026-09-23, and on 2026 variants already present in November 2022 (pre2022) or added later (post2022). The common-scale block expresses each threshold as the percentage of 2026 gnomAD intersection variants meeting it (post hoc). Set sizes (P/LP, B/LB) were ",
           paste(sprintf("%s %d / %d", ss$set, ss$P_LP, ss$B_LB), collapse = "; "), "."),
    rbind(pp[, block := "thresholds"], bb[, block := "thresholds"], sc[, block := "common scale (% of gnomAD variants)"], fill = TRUE))

## S9 gene classes ---------------------------------------------------------------------------------
cs <- rd("gene_classes", "class_sizes.tsv"); tc <- rd("gene_classes", "thresholds_by_class.tsv")
cc <- rd("gene_classes", "class_comparisons.tsv"); rc <- rd("gene_classes", "reclassification.tsv")
add("S9", "Calibration by gene class",
    "Class sizes, thresholds with 2.5th-97.5th percentile bootstrap intervals, pairwise comparisons (classes considered different when intervals do not overlap) and the percentage of reference variants whose evidence level changes with class-specific instead of global thresholds.",
    rbind(cs[, block := "class sizes"], tc[, block := "thresholds"], cc[, block := "pairwise comparisons"], rc[, block := "reclassification"], fill = TRUE))

## S10 VUS -----------------------------------------------------------------------------------------
add("S10", "Computational evidence for ClinVar missense variants of uncertain significance",
    "Number and percentage of 1,208,461 ClinVar missense VUS (at least one star, MANE Select, disease genes) at each evidence level for each predictor, with ACMG/AMP points, using main-analysis thresholds (CADD, AVI) or intersection thresholds (all five), for all VUS and for VUS with AF below 0.01; the second block cross-tabulates CADD and AVI levels.",
    rbind(rd("vus", "table_vus_evidence.tsv")[, block := "evidence levels"], rd("vus", "table_vus_cadd_vs_avi.tsv")[, block := "CADD x AVI"], fill = TRUE))

## S11 splice region -------------------------------------------------------------------------------
add("S11", "Calibration of CADD v1.7 and AVI in the intronic splice region",
    "Thresholds (one seed, 10,000 bootstrap replicates, missense prior 0.0441) for SNVs 3 to 8 nucleotides into introns of MANE Select transcripts, from 2,392 P/LP and 68,917 B/LB ClinVar variants and 175,221 rare gnomAD variants at the same positions and genes, compared with the missense thresholds.",
    rbind(rd("noncoding", "splice_region_intronic", "thresholds.tsv")[, block := "splice-region thresholds"],
          rd("noncoding", "splice_region_intronic", "comparison_with_missense.tsv")[, block := "comparison with missense"], fill = TRUE))

## S12 noncoding classes and MPRA ------------------------------------------------------------------
ev <- rd("noncoding", "evaluation_missense_thresholds_by_class.tsv")
add("S12", "Missense-calibrated thresholds in noncoding and synonymous variant classes, and MPRA",
    "Reference variants per region class; AUROC of CADD and AVI per class; likelihood ratios (with 1,000-replicate bootstrap intervals) of each missense-calibrated evidence interval in each class and whether the conservative bound meets the required ratio; AUROC for functional versus non-functional SNVs in saturation-mutagenesis MPRA data (Kircher et al. 2019).",
    rbind(rd("noncoding", "class_sizes.tsv")[, block := "class sizes"], ev[, block := "interval likelihood ratios"],
          rd("noncoding", "mpra_auroc.tsv")[, block := "MPRA AUROC"], fill = TRUE))

## Manuscript v6 (2026-10-05): the population table (MCPS and BRAVO) is no longer part of this article
## (population analyses reserved for a separate study); later tables are renumbered S13-S18.

## S13 MAVE score sets -----------------------------------------------------------------------------
qc <- rd("mave", "scoreset_qc.tsv"); inv <- unique(rd("mave", "mavedb_human_scoreset_inventory.tsv")[, .(urn, title, licence, pmid)], by = "urn")
add("S13", "MaveDB score sets used for validation",
    "The 89 MaveDB score sets with normal and abnormal functional classes, with gene, calibration used, mapping level, numbers of scored and classified variants, protein reference-residue match rate, licence and primary publication (PubMed ID); the second block gives classified MANE missense SNVs per gene after excluding discordant classifications. Please cite MaveDB and the primary publication of each score set.",
    rbind(merge(qc, inv, by = "urn", all.x = TRUE)[, block := "score sets"], rd("mave", "labels_per_gene.tsv")[, block := "labelled SNVs per gene"], fill = TRUE))

## S14 CADD v1.6 vs v1.7 ---------------------------------------------------------------------------
cd <- file.path("cadd16_vs_17")
add("S14", "CADD v1.6 versus v1.7 on identical reference data",
    "Thresholds (median and range over five seeds) of CADD v1.6 and v1.7 calibrated on the same 2026 reference sets, AUROC with 95% bootstrap intervals and the paired difference, the percentage of reference variants whose evidence level changes between versions, and the full cross-tabulation of evidence levels (each version at its own thresholds).",
    rbind(rd(cd, "thresholds_median.tsv")[, block := "thresholds"], rd(cd, "discrimination.tsv")[, block := "discrimination"],
          rd(cd, "evidence_changes.tsv")[, block := "evidence level changes"], rd(cd, "evidence_crosstab.tsv")[, block := "cross-tabulation (y = 1 P/LP, 0 B/LB)"], fill = TRUE))

## S15 analysis plan and amendments -----------------------------------------------------------------
add("S15", "Prespecified analysis plan and dated amendments",
    "The analysis plan was fixed on 25 September 2026 before any CADD v1.7 or AVI calibration result and tagged in the version-controlled repository (analysis-plan-v1). Each later analysis was specified in a dated amendment committed before it was run. The Timing column states whether the amendment was motivated by results already seen, and the last column whether the analysis is reported in this article. Post hoc analyses that were not specified in any amendment are listed at the end.",
    data.table(
      Item = c(sprintf("Amendment %d", 1:15), "Post hoc", "Post hoc"),
      Date = c(rep("2026-09-25", 6), rep("2026-09-26", 5), rep("2026-09-27", 3), "2026-10-04", "2026-09-26", "2026-09-27"),
      Analysis = c("Paired bootstrap comparison of AUROC and AUPRC between AVI and CADD v1.7",
                   "Temporal hold-out validation (calibration on ClinVar 2024-01-27, test on later variants)",
                   "Calibration of REVEL, BayesDel and AlphaMissense on the five-predictor intersection",
                   "Population frequency and ancestry of PP3-level variants in gnomAD v4.1",
                   "Computational evidence for ClinVar missense variants of uncertain significance",
                   "Threshold drift across ClinVar releases and before/after November 2022 (AVI as negative control)",
                   "Exclusion of HGMD Public positions from the 2019 and 2026 reference sets",
                   "Calibration by gene class (constraint, inheritance, mechanism)",
                   "Genome-wide catalogue of SNVs at PP3 levels",
                   "Noncoding variants (ClinVar region classes, calibration of classes with enough variants, MPRA)",
                   "SpliceVarDB analysis deferred to a later study (database not reachable)",
                   "Population frequencies in the Mexico City Prospective Study and TOPMed BRAVO",
                   "Validation with MAVE functional classes from MaveDB",
                   "CADD v1.6 versus v1.7 on identical reference data",
                   "Calibrated evidence for VUS of January 2024 resolved by September 2026 (requested by co-author review)",
                   "Drift thresholds expressed on a common scale (percentage of gnomAD variants meeting each threshold)",
                   "Source of MCPS genotypes (array, exome or genome) for variants with frequency conflicts"),
      Timing = c("after the main results", rep("before running", 4), "before running; motivated by the comparator results",
                 "before running; motivated by the drift results", rep("before running", 3), "scope change before running",
                 "before any MCPS or BRAVO result", "before any MAVE score was downloaded", "before any CADD v1.6 score was read", "before any result; motivated by co-author review",
                 "post hoc", "post hoc"),
      Reported_in_this_article = c("yes", "yes", "yes", "no (population analyses reserved for a separate study)", "yes", "yes", "yes", "yes",
                                   "no (companion resource article)", "yes", "no (deferred)", "no (population analyses reserved for a separate study)",
                                   "yes", "yes", "yes", "yes", "no (population analyses reserved for a separate study)")))

## S16-S18 (added for manuscript v5, 2026-10-04; renumbered in v6) ------------------------------------------------------------
add("S16", "Calibrated evidence intervals on MAVE functional classes",
    "For each predictor and ClinVar-calibrated evidence interval, the likelihood ratio P(interval | functionally abnormal) / P(interval | functionally normal) with 95% bootstrap interval (2,000 replicates) and abnormal/normal variants in the interval, and the pooled AUROC, on 93,154 MANE missense SNVs in 39 genes (18,936 abnormal, 74,218 normal) from 89 MaveDB score sets. Functional abnormality is not pathogenicity; the ratios test direction and ordering, not the ClinGen targets.",
    fread(file.path(F, "Table2_mave.tsv")))
vr <- file.path("vus_resolved")
add("S17", "Calibrated evidence for ClinVar VUS of January 2024 that were resolved by September 2026",
    "Amendment A15. VUS (>= 1 star, MANE Select missense, disease gene) in the ClinVar release of 27 January 2024, followed to the release of 23 September 2026. Thresholds: r2024 = calibrated on the January 2024 release only (five-predictor intersection, primary); main2026 and inter2026 = 2026 thresholds (secondary). Blocks: set sizes; percentage of variants with PP3 or BP4 evidence and mean ACMG/AMP points by 2026 outcome; concordance of the evidence with the 2026 classification; likelihood ratio of each interval, P(interval | resolved P/LP) / P(interval | resolved B/LB), with 95% bootstrap interval (2,000 replicates), the ratio required and whether the conservative bound meets it. Subset AF<0.01 (gnomAD v4.1 group maximum) is primary.",
    rbind(rd(vr, "set_sizes.tsv")[, block := "set sizes"], rd(vr, "table_direction.tsv")[, block := "direction by outcome"],
          rd(vr, "table_concordance.tsv")[, block := "concordance (r2024, AF<0.01)"], rd(vr, "table_interval_LR.tsv")[, block := "interval likelihood ratios"], fill = TRUE))
add("S18", "Illustrative resolved VUS",
    "Variants selected by the rule fixed in amendment A15 before inspection: VUS in January 2024 whose September 2026 record was reviewed by an expert panel or practice guideline, up to three per direction with the largest absolute AVI points (r2024 thresholds), ties broken by genomic position. Scores and evidence levels of all five predictors with the r2024 thresholds.",
    fread(file.path(R, vr, "table_examples.tsv"), colClasses = list(character = "chrom")))

## S19-S20 (added for manuscript v7, 2026-10-29; amendment A16) ---------------------------------------
vc <- file.path(R, "vus_candidates")
add("S19", "ClinVar missense VUS reviewed by an expert panel with concordant calibrated computational evidence",
    "Amendment A16. VUS in the ClinVar release of 23 September 2026 (MANE Select missense, disease gene) whose record was reviewed by an expert panel or practice guideline and for which REVEL, BayesDel and AlphaMissense (published ClinGen thresholds) and AVI (Table 1) all gave evidence in the same direction. Comp points, the smallest absolute ACMG/AMP points among the three missense predictors. Scenario, the other evidence with which these points would complete a likely pathogenic (>= +6 points) or likely benign (<= -1 point) classification under the ACMG/AMP point system; it is not a reclassification, because ClinVar does not record the points of the remaining evidence and expert panels may already have applied their own computational criteria. Ev CADD16, CADD v1.7 scores read with the published CADD v1.6 thresholds; Flip AVI, evidence in opposite directions under those thresholds and under calibrated AVI. AF max, gnomAD v4.1 group maximum allele frequency. Ordered by absolute points and genomic position. The 483,690 VUS with concordant evidence at any review status are in the repository (calibration/results/v17_avi/vus_candidates/vus_concordant_evidence.tsv.gz).",
    fread(file.path(vc, "vus_concordant_expert_panel.tsv"), colClasses = list(character = "chrom")))
add("S20", "Concordant calibrated evidence and transferred thresholds for all ClinVar missense VUS",
    "Amendment A16. The 1,208,430 rare (gnomAD v4.1 AF max < 0.01) MANE Select missense VUS with at least one star in disease genes, ClinVar release of 23 September 2026. Blocks: thresholds used (CADD and AVI, Table 1; REVEL, BayesDel and AlphaMissense, published ClinGen values; CADD16, published CADD v1.6 values applied to CADD v1.7 scores); number of VUS by concordant ACMG/AMP points, for all VUS and for those reviewed by an expert panel or practice guideline (3-4 stars), and number whose evidence level or direction differs between transferred CADD v1.6 thresholds and calibrated CADD v1.7 or AVI thresholds; cross-tabulations of evidence levels under transferred CADD v1.6 thresholds against calibrated CADD v1.7 and against calibrated AVI.",
    rbind(rd("vus_candidates", "thresholds_used.tsv")[, block := "thresholds used"], rd("vus_candidates", "summary.tsv")[, block := "counts"],
          rd("vus_candidates", "crosstab_cadd16_vs_cadd17.tsv")[, block := "CADD v1.6 transferred versus CADD v1.7 calibrated"],
          rd("vus_candidates", "crosstab_cadd16_vs_avi.tsv")[, block := "CADD v1.6 transferred versus AVI calibrated"], fill = TRUE))

## Numbering in order of first citation (manuscript v9, 2026-10-29) -------------------------------
# The sections above keep the ids under which each table was added. The manuscript numbers supplementary
# tables in the order in which they are first cited in the text, so the published ids follow this map
# (old id = new id). Table 3 and Table 4 of the main text were swapped for the same reason.
renum <- c(S15 = "S1", S3 = "S2", S1 = "S3", S2 = "S4", S4 = "S5", S6 = "S6", S7 = "S7", S13 = "S8", S16 = "S9",
           S17 = "S10", S18 = "S11", S10 = "S12", S20 = "S13", S19 = "S14", S11 = "S15", S12 = "S16", S14 = "S17",
           S5 = "S18", S8 = "S19", S9 = "S20")
stopifnot(setequal(names(renum), names(tabs)), !anyDuplicated(renum))
tabs <- setNames(tabs[names(renum)], renum); meta <- setNames(meta[names(renum)], renum)
for (id in names(meta)) meta[[id]]$Table <- id
# cross-references to supplementary tables inside legends and cells follow the new numbers
relabel <- function(x) { m <- gregexpr("Table S[0-9]+", x); regmatches(x, m) <- lapply(regmatches(x, m), function(v) paste0("Table ", renum[sub("^Table ", "", v)])); x }
for (id in names(meta)) meta[[id]]$Legend <- relabel(meta[[id]]$Legend)
for (id in names(tabs)) for (j in names(tabs[[id]])[vapply(tabs[[id]], is.character, TRUE)]) set(tabs[[id]], j = j, value = relabel(tabs[[id]][[j]]))

## Write -------------------------------------------------------------------------------------------
wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold")
contents <- rbindlist(meta)
addWorksheet(wb, "Contents"); writeData(wb, "Contents", contents, headerStyle = hs)
setColWidths(wb, "Contents", cols = 1:3, widths = c(8, 60, 120)); addStyle(wb, "Contents", createStyle(wrapText = TRUE, valign = "top"), rows = 2:(nrow(contents) + 1), cols = 1:3, gridExpand = TRUE)
# Layout (v6, 2026-10-05, authors found the sheets hard to read): tables that combine several blocks
# with different columns are no longer written as one wide sheet full of empty cells. Each block is
# written as its own sub-table (bold subheading, its own header row, only the columns it uses),
# separated by a blank row; the title and legend sit above, column names are made readable, and
# column widths follow the content. The TSV copies keep the long machine-readable format.
title_st <- createStyle(textDecoration = "bold", fontSize = 12)
leg_st <- createStyle(wrapText = TRUE, valign = "top")
sub_st <- createStyle(textDecoration = "bold", fontSize = 11)
head_st <- createStyle(textDecoration = "bold", fgFill = "#E7E6E6", border = "Bottom", wrapText = TRUE, valign = "center")
nice <- function(x) { x <- sub("^lo$|_lo$", " 95% CI lower", x); x <- sub("^hi$|_hi$", " 95% CI upper", x); x <- gsub("(^|_| )P_?LP($|_| )", "\\1P/LP\\2", x); x <- gsub("(^|_| )B_?LB($|_| )", "\\1B/LB\\2", x); trimws(gsub("_", " ", x)) }
for (id in names(tabs)) {
  d <- tabs[[id]]
  fwrite(d, file.path(out_dir, sprintf("Table%s.tsv", id)), sep = "\t")
  blocks <- if ("block" %in% names(d)) split(d, factor(d$block, levels = unique(d$block))) else list(main = d)
  addWorksheet(wb, id, gridLines = FALSE)
  writeData(wb, id, paste0("Table ", id, ". ", meta[[id]]$Title), startRow = 1); addStyle(wb, id, title_st, rows = 1, cols = 1)
  writeData(wb, id, meta[[id]]$Legend, startRow = 2)
  ncol_max <- max(vapply(blocks, function(b) sum(vapply(b, function(v) any(!is.na(v) & v != ""), TRUE)), 0L))
  mergeCells(wb, id, cols = 1:max(6, ncol_max), rows = 2); addStyle(wb, id, leg_st, rows = 2, cols = 1)
  setRowHeights(wb, id, rows = 2, heights = 15 * ceiling(nchar(meta[[id]]$Legend) / 140))
  r <- 4; widths <- numeric(0)
  for (nm in names(blocks)) {
    b <- copy(blocks[[nm]]); if ("block" %in% names(b)) b[, block := NULL]
    keep <- names(b)[vapply(b, function(v) any(!is.na(v) & as.character(v) != ""), TRUE)]
    b <- b[, ..keep]
    idc <- intersect(c("thresholds", "set", "analysis", "subset", "class", "criterion", "tool", "predictor", "level", "evidence", "interval", "outcome"), names(b))
    setcolorder(b, c(idc, setdiff(names(b), idc)))   # identifying columns first
    setnames(b, vapply(nice(names(b)), function(s) paste0(toupper(substr(s, 1, 1)), substring(s, 2)), ""))
    if (length(blocks) > 1) { writeData(wb, id, paste0(toupper(substr(nm, 1, 1)), substring(nm, 2)), startRow = r); addStyle(wb, id, sub_st, rows = r, cols = 1); r <- r + 1 }
    writeData(wb, id, b, startRow = r, headerStyle = head_st, na.string = "")
    w <- pmin(40, pmax(8, vapply(seq_along(b), function(j) max(nchar(c(names(b)[j], as.character(b[[j]]))), na.rm = TRUE) + 2, 0)))
    widths <- pmax(c(widths, rep(0, max(0, length(w) - length(widths)))), c(w, rep(0, max(0, length(widths) - length(w)))))
    r <- r + nrow(b) + 2
  }
  setColWidths(wb, id, cols = seq_along(widths), widths = widths)
  if (length(blocks) == 1) freezePane(wb, id, firstActiveRow = 5)
}
fwrite(contents, file.path(out_dir, "contents.tsv"), sep = "\t")
saveWorkbook(wb, file.path(F, "Supplementary_Tables_S1-S20.xlsx"), overwrite = TRUE)
cat("Written", length(tabs), "tables:", file.path(F, "Supplementary_Tables_S1-S20.xlsx"), "\n")
print(contents[, .(Table, Title, rows = vapply(tabs, nrow, 0L))])
