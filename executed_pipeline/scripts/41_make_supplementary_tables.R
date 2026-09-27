# Supplementary Tables S1-S16 of the manuscript, assembled from the analysis outputs.
# Each table is written as TSV to figures/data/SupplementaryTables/ and all of them, with a
# "Contents" sheet giving titles and legends, to figures/Supplementary_Tables_S1-S16.xlsx.
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
    "Thresholds (seed 20260925) in the main analysis and in the pre-specified sensitivity analyses S1 (review status of at least two stars), S2 (rarity defined as AF below 0.001), S3 (gene-level bootstrap) and S4 (exclusion of missense SNVs within 3 nucleotides of an internal exon boundary). Empty cells indicate levels not reached.",
    setcolorder(rbind(p, b, fill = TRUE), c("criterion", "tool", "level")))

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

## S13 population resources ------------------------------------------------------------------------
add("S13", "PP3-level missense variants in the Mexico City Prospective Study and TOPMed BRAVO",
    "Summary counts per predictor and PP3 level in MCPS (overall and Indigenous Mexican ancestry, IMX) and BRAVO, and the 74 PP3-level SNVs with MCPS AF above 1% and at or below 1% in every gnomAD v4.1 group, with MCPS variant source (array, exome or genome; added post hoc). MCPS data from the Regeneron Genetics Center and collaborators (Ziyatdinov et al. 2023).",
    rbind(rd("mcps", "pp3_levels_mcps_summary.tsv")[, block := "MCPS summary"], rd("bravo", "pp3_levels_bravo_summary.tsv")[, block := "BRAVO summary"],
          fread(file.path(R, "mcps", "mcps_only_frequency_conflicts.tsv"), colClasses = list(character = "chrom"))[, block := "MCPS-only frequency conflicts"], fill = TRUE))

## S14 MAVE score sets -----------------------------------------------------------------------------
qc <- rd("mave", "scoreset_qc.tsv"); inv <- unique(rd("mave", "mavedb_human_scoreset_inventory.tsv")[, .(urn, title, licence, pmid)], by = "urn")
add("S14", "MaveDB score sets used for validation",
    "The 89 MaveDB score sets with normal and abnormal functional classes, with gene, calibration used, mapping level, numbers of scored and classified variants, protein reference-residue match rate, licence and primary publication (PubMed ID); the second block gives classified MANE missense SNVs per gene after excluding discordant classifications. Please cite MaveDB and the primary publication of each score set.",
    rbind(merge(qc, inv, by = "urn", all.x = TRUE)[, block := "score sets"], rd("mave", "labels_per_gene.tsv")[, block := "labelled SNVs per gene"], fill = TRUE))

## S15 CADD v1.6 vs v1.7 ---------------------------------------------------------------------------
cd <- file.path("cadd16_vs_17")
add("S15", "CADD v1.6 versus v1.7 on identical reference data",
    "Thresholds (median and range over five seeds) of CADD v1.6 and v1.7 calibrated on the same 2026 reference sets, AUROC with 95% bootstrap intervals and the paired difference, the percentage of reference variants whose evidence level changes between versions, and the full cross-tabulation of evidence levels (each version at its own thresholds).",
    rbind(rd(cd, "thresholds_median.tsv")[, block := "thresholds"], rd(cd, "discrimination.tsv")[, block := "discrimination"],
          rd(cd, "evidence_changes.tsv")[, block := "evidence level changes"], rd(cd, "evidence_crosstab.tsv")[, block := "cross-tabulation (y = 1 P/LP, 0 B/LB)"], fill = TRUE))

## S16 analysis plan and amendments -----------------------------------------------------------------
add("S16", "Prespecified analysis plan and dated amendments",
    "The analysis plan was fixed on 25 September 2026 before any CADD v1.7 or AVI calibration result and tagged in the version-controlled repository (analysis-plan-v1). Each later analysis was specified in a dated amendment committed before it was run. The last column states whether the amendment was motivated by results already seen. Post hoc analyses that were not specified in any amendment are listed at the end.",
    data.table(
      Item = c(sprintf("Amendment %d", 1:14), "Post hoc", "Post hoc"),
      Date = c(rep("2026-09-25", 6), rep("2026-09-26", 5), rep("2026-09-27", 3), "2026-09-26", "2026-09-27"),
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
                   "Drift thresholds expressed on a common scale (percentage of gnomAD variants meeting each threshold)",
                   "Source of MCPS genotypes (array, exome or genome) for variants with frequency conflicts"),
      Timing = c("after the main results", rep("before running", 4), "before running; motivated by the comparator results",
                 "before running; motivated by the drift results", rep("before running", 3), "scope change before running",
                 "before any MCPS or BRAVO result", "before any MAVE score was downloaded", "before any CADD v1.6 score was read",
                 "post hoc", "post hoc")))

## Write -------------------------------------------------------------------------------------------
wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold")
contents <- rbindlist(meta)
addWorksheet(wb, "Contents"); writeData(wb, "Contents", contents, headerStyle = hs)
setColWidths(wb, "Contents", cols = 1:3, widths = c(8, 60, 120)); addStyle(wb, "Contents", createStyle(wrapText = TRUE, valign = "top"), rows = 2:(nrow(contents) + 1), cols = 1:3, gridExpand = TRUE)
for (id in names(tabs)) {
  d <- tabs[[id]]; if ("block" %in% names(d)) setcolorder(d, c("block", setdiff(names(d), "block")))
  addWorksheet(wb, id); writeData(wb, id, paste0("Table ", id, ". ", meta[[id]]$Title), startRow = 1)
  addStyle(wb, id, hs, rows = 1, cols = 1); writeData(wb, id, d, startRow = 3, headerStyle = hs)
  freezePane(wb, id, firstActiveRow = 4); setColWidths(wb, id, cols = seq_len(ncol(d)), widths = "auto")
  fwrite(d, file.path(out_dir, sprintf("Table%s.tsv", id)), sep = "\t")
}
fwrite(contents, file.path(out_dir, "contents.tsv"), sep = "\t")
saveWorkbook(wb, file.path(F, "Supplementary_Tables_S1-S16.xlsx"), overwrite = TRUE)
cat("Written", length(tabs), "tables:", file.path(F, "Supplementary_Tables_S1-S16.xlsx"), "\n")
print(contents[, .(Table, Title, rows = vapply(tabs, nrow, 0L))])
