# Results (draft, generated from executed analyses)

> Numbers are copied from `calibration/results/v17_avi/summary/` (script
> `executed_pipeline/scripts/15_summarize_calibration.R`, log `logs/15_summarize_calibration.log`).
> All analyses follow the pre-specified plan (`STUDY_PLAN.md`, tag `analysis-plan-v1`).
> Interpretation is provisional until the remaining aims are completed.

## R1. Calibration of CADD v1.7 and AVI (main analysis)

Reference set: 59,979 P/LP and 77,128 B/LB rare missense ClinVar variants (4,826 genes); gnomAD
set 1,955,995 variants. Thresholds: one-sided 95% bootstrap bound, B = 10,000, median [range]
over five seeds. PP3: score ≥ threshold; BP4: score ≤ threshold.

| Predictor | PP3 Supporting | PP3 Moderate | PP3 Strong | BP4 Supporting | BP4 Moderate | BP4 Strong |
|---|---|---|---|---|---|---|
| CADD v1.7 PHRED | 25.6 | 28.0 | not reached | 22.6 [22.6–22.7] | 19.59 | 8.341 |
| AVI Phred | 25.18 | 26.44 | not reached (reached in 58.7% of replicates) | 22.95 | 20.52 [20.52–20.53] | 9.15 [9.15–9.17] |
| *CADD v1.6 (Pejaver et al. 2022, for reference)* | *25.3* | *28.1* | *not reached* | *22.7* | *17.3* | *0.15* |

Very Strong was not reached at either end by either predictor. Thresholds were essentially
invariant across seeds (maximum range 0.1 PHRED units), in contrast with the CADD v1.6
reproduction (BP4 Strong 0.15–0.31), reflecting a reference set ~10 times larger.
The AVI pathogenic posterior plateaued at ~0.60, just below the Strong posterior (0.607); CADD
plateaued at ~0.37 (Figure 1).

## R2. Sensitivity analyses (seed 20260925)

| Predictor | Level | Main | S1 ≥2 stars | S2 AF<0.001 | S3 gene bootstrap | S4 no exon edge |
|---|---|---|---|---|---|---|
| CADD | PP3 Supporting | 25.6 | 25.9 | 25.6 | 25.9 | 25.5 |
| CADD | PP3 Moderate | 28.0 | 28.2 | 28.0 | 28.2 | 27.8 |
| CADD | BP4 Supporting | 22.7 | 22.6 | 22.6 | 22.4 | 22.5 |
| CADD | BP4 Moderate | 19.59 | 20.60 | 19.56 | 19.49 | 19.54 |
| CADD | BP4 Strong | 8.341 | 1.74 | 8.342 | 7.773 | 8.302 |
| AVI | PP3 Supporting | 25.18 | 25.04 | 25.20 | 25.32 | 25.13 |
| AVI | PP3 Moderate | 26.44 | 25.99 | 26.57 | 26.60 | 26.37 |
| AVI | PP3 Strong | not reached | **30.23** | not reached | not reached | not reached |
| AVI | BP4 Supporting | 22.95 | 22.91 | 22.97 | 22.89 | 22.91 |
| AVI | BP4 Moderate | 20.52 | 20.81 | 20.45 | 20.16 | 20.36 |
| AVI | BP4 Strong | 9.15 | 9.17 | 9.08 | 8.98 | 9.01 |

Supporting and Moderate thresholds varied by ≤ 1.1 PHRED units across all sensitivity analyses
(largest change: CADD BP4 Moderate 19.59 → 20.60 with ≥ 2-star variants).
Restricting to AF < 0.001 (S2) did not change AVI's BP4 thresholds (≤ 0.07 units), giving no
evidence that the overlap with AVI's frequency-based training labels favoured AVI. Gene-level
resampling (S3) changed thresholds by ≤ 0.6 units. With ≥ 2-star variants (S1), AVI reached PP3
Strong (≥ 30.23) and CADD's BP4 Strong moved to 1.74 (fewer B/LB variants at low scores).

## R3. Discrimination (main reference set; 95% bootstrap intervals, 1,000 replicates)

| Predictor | AUROC | AUPRC (prevalence 0.4375) |
|---|---|---|
| CADD v1.7 | 0.902 [0.900–0.903] | 0.849 [0.846–0.853] |
| AVI | 0.940 [0.939–0.942] | 0.917 [0.914–0.919] |

Paired comparison (added after analysis plan v1; 1,000 paired bootstrap replicates, seed
20260925; script `17_paired_discrimination.R`): AVI − CADD difference in AUROC +0.039
[0.037–0.040] and in AUPRC +0.068 [0.065–0.071]; the difference was ≤ 0 in none of the 1,000
replicates for either metric (bootstrap p < 0.001).

## R4. Evidence levels assigned to the reference variants (main median thresholds; in-sample)

| Evidence | CADD P/LP | CADD B/LB | AVI P/LP | AVI B/LB |
|---|---|---|---|---|
| PP3 Moderate | 46.2% | 4.8% | 67.9% | 4.0% |
| PP3 Supporting | 26.1% | 6.6% | 12.1% | 3.1% |
| Indeterminate | 21.2% | 19.6% | 12.7% | 13.1% |
| BP4 Supporting | 4.1% | 16.1% | 4.3% | 18.7% |
| BP4 Moderate | 2.1% | 29.1% | 2.7% | 38.8% |
| BP4 Strong | 0.3% | 23.8% | 0.3% | 22.4% |

These proportions are computed on the calibration set itself; out-of-sample performance is
assessed in R5.

## R6. Head-to-head comparison with already-recommended predictors (Aim 6; amendment A3)

Intersection of reference variants scored by all five predictors: 59,034 P/LP and 75,306 B/LB
(4,725 genes); gnomAD set 1,888,558. Five seeds × B = 10,000 per predictor (script
`21_calibrate_comparators.R`). Thresholds (median [range]; PP3 score ≥, BP4 score ≤):

| Predictor | PP3 Sup | PP3 Mod | PP3 Strong | BP4 Sup | BP4 Mod | BP4 Strong |
|---|---|---|---|---|---|---|
| REVEL | 0.642 | 0.752 [0.752–0.753] | 0.898 | 0.396 [0.395–0.396] | 0.272 | 0.087 |
| BayesDel noAF | 0.14 | 0.27 | 0.48 | −0.13 | −0.25 | −0.54 |
| AlphaMissense | 0.6698 [0.6698–0.67] | 0.8605 [0.8603–0.8607] | 0.9898 | 0.2367 [0.2367–0.2368] | 0.14 | 0.0763 |
| AVI | 25.19 | 26.44 | not reached | 22.98 | 20.61 | 9.28 [9.26–9.28] |
| CADD v1.7 | 25.6 | 28.0 | not reached | 22.7 [22.6–22.7] | 19.59 | 8.383 |

Very Strong was not reached by any predictor. CADD and AVI thresholds on the intersection are
within 0.3 units of those on the full main set.
Comparison with published thresholds (Pejaver et al. 2022, Table 2; values provided by the
author from the article on 2026-09-25):

| Predictor | Level | Published (ClinVar 2019) | This study (ClinVar 2026, intersection) |
|---|---|---|---|
| BayesDel | PP3 Supporting | ≥ 0.13 | ≥ 0.14 |
| BayesDel | PP3 Moderate | ≥ 0.27 | ≥ 0.27 |
| BayesDel | PP3 Strong | ≥ 0.50 | ≥ 0.48 |
| BayesDel | BP4 Supporting | ≤ −0.18 | ≤ −0.13 |
| BayesDel | BP4 Moderate | ≤ −0.36 | ≤ −0.25 |
| BayesDel | BP4 Strong | not reached | ≤ −0.54 |
| REVEL | PP3 Supporting | ≥ 0.644 | ≥ 0.642 |
| REVEL | PP3 Moderate | ≥ 0.773 | ≥ 0.752 |
| REVEL | PP3 Strong | ≥ 0.932 | ≥ 0.898 |
| REVEL | BP4 Supporting | ≤ 0.290 | ≤ 0.396 |
| REVEL | BP4 Moderate | ≤ 0.183 | ≤ 0.272 |
| REVEL | BP4 Strong | ≤ 0.016 | ≤ 0.087 |
| REVEL | BP4 Very Strong | ≤ 0.003 | not reached |

PP3 Supporting thresholds were reproduced almost exactly, whereas higher PP3 levels and all BP4
levels were less stringent with current data. Candidate explanations (to be tested, proposed
amendment A6): classifications after the 2022 ClinGen recommendations using these predictors as
PP3/BP4 evidence (circularity); inclusion of predictor training variants and HGMD variants, which
the original study excluded; changes in the composition of ClinVar.
AlphaMissense (Bergquist et al. 2025; values provided by the author from the article on
2026-09-25) versus this study: PP3 Very Strong ≥ 0.990 vs not reached; PP3 Strong ≥ 0.972 vs
≥ 0.9898; PP3 Moderate ≥ 0.906 vs ≥ 0.8605; PP3 Supporting ≥ 0.792 vs ≥ 0.6698; BP4 Supporting
≤ 0.169 vs ≤ 0.2367; BP4 Moderate ≤ 0.099 vs ≤ 0.14; BP4 Strong ≤ 0.070 vs ≤ 0.0763. Lower-level
thresholds were again less stringent with current data, while PP3 Strong was more stringent and
Very Strong was not reached.

Discrimination on the intersection (95% bootstrap CI) and paired differences versus AVI (1,000
paired replicates):

| Predictor | AUROC | AUPRC | ΔAUROC vs AVI [95% CI] | ΔAUPRC vs AVI [95% CI] |
|---|---|---|---|---|
| REVEL | 0.961 [0.960–0.962] | 0.952 [0.951–0.954] | +0.020 [0.019–0.021] | +0.035 [0.033–0.037] |
| BayesDel | 0.957 [0.956–0.958] | 0.946 [0.945–0.948] | +0.017 [0.015–0.018] | +0.029 [0.027–0.031] |
| AlphaMissense | 0.949 [0.948–0.950] | 0.935 [0.933–0.937] | +0.008 [0.007–0.009] | +0.018 [0.016–0.019] |
| AVI | 0.941 [0.939–0.942] | 0.918 [0.915–0.920] | — | — |
| CADD v1.7 | 0.901 [0.900–0.903] | 0.849 [0.846–0.852] | −0.039 [−0.040 to −0.038] | −0.068 [−0.071 to −0.066] |

All paired differences had the same sign in all 1,000 replicates.

Evidence levels on the intersection (in-sample; % of P/LP / % of B/LB):

| Predictor | PP3 Strong | PP3 Moderate | PP3 Supporting | Indeterminate | BP4 Sup | BP4 Mod | BP4 Strong |
|---|---|---|---|---|---|---|---|
| REVEL | 45.6 / 0.5 | 27.1 / 2.0 | 10.1 / 2.6 | 11.3 / 12.0 | 3.1 / 12.5 | 2.5 / 38.3 | 0.4 / 32.0 |
| BayesDel | 41.1 / 0.5 | 29.3 / 2.2 | 11.4 / 3.1 | 12.8 / 14.1 | 2.8 / 11.6 | 2.4 / 38.3 | 0.3 / 30.2 |
| AlphaMissense | 33.4 / 0.5 | 34.0 / 2.8 | 12.2 / 3.1 | 14.7 / 14.7 | 3.1 / 13.3 | 2.3 / 35.9 | 0.4 / 29.7 |
| AVI | — | 68.0 / 4.0 | 12.0 / 3.0 | 12.7 / 12.9 | 4.3 / 18.3 | 2.7 / 39.1 | 0.3 / 22.7 |
| CADD v1.7 | — | 46.1 / 4.8 | 26.2 / 6.6 | 20.9 / 19.0 | 4.3 / 16.7 | 2.1 / 29.0 | 0.3 / 23.9 |

Interpretation (provisional): for missense variants, the missense-specialised predictors
(REVEL, BayesDel, AlphaMissense) discriminate better than the genome-wide AVI and reach PP3
Strong, which AVI and CADD do not; AVI clearly outperforms CADD v1.7. AVI's advantage over the
missense-specialised tools, if any, would be expected outside missense variants (deferred Aim 9).

## R7. Population frequency of variants with computational evidence (Aim 8; amendment A4)

Universe: 8,214,534 MANE missense SNVs observed in gnomAD v4.1. PP3 levels with the main-analysis
median thresholds (CADD ≥ 25.6 / ≥ 28.0; AVI ≥ 25.18 / ≥ 26.44). Script `22_ancestry.R`.

| Predictor | PP3 level | Observed in gnomAD | Group AF > 1% | Group AF > 5% |
|---|---|---|---|---|
| CADD v1.7 | Supporting or higher | 1,756,432 | 7,399 (0.42%) | 1,529 (0.087%) |
| CADD v1.7 | Moderate | 816,315 | 2,914 (0.36%) | 553 (0.068%) |
| AVI | Supporting or higher | 1,430,434 | 3,950 (0.28%) | 857 (0.060%) |
| AVI | Moderate | 930,564 | 2,064 (0.22%) | 447 (0.048%) |

Among PP3-level variants common (> 1%) in at least one group, the group with the highest
frequency was rarely NFE (3.5–3.7% of such variants for both predictors and levels); the most
frequent grpmax groups were AFR (17–22%), AMI (15–17%), FIN (13–14%), MID (11–13%), EAS (11–12%)
and ASJ (11–12%). Variants common (> 1%) in only one group ("group-specific common"), CADD
Supporting+ / AVI Supporting+: AFR 1,286 / 567; AMI 710 / 445; MID 635 / 311; EAS 639 / 306;
FIN 593 / 352; ASJ 557 / 286; SAS 303 / 164; AMR 211 / 117; NFE 72 / 41. AVI assigned PP3 to
fewer common variants than CADD, consistent with its frequency-based training labels.

Per-group proportions of PP3-level variants exceeding 1% (pre-specified outcome 3; full table
`calibration/results/v17_avi/ancestry/table_per_group.tsv`) ranged from 0.08–0.13% (NFE) to
2.3–3.2% (ASJ, MID) and ~40–44% (AMI). **Caution (post hoc note):** this proportion is strongly
affected by each group's sample size — in small cohorts (AMI, ASJ, MID) only variants of
appreciable frequency are observed at all, and AF > 1% is reached with few allele counts — so it
should not be compared across groups as a measure of risk; the grpmax distribution and the
group-specific counts are less affected. Allele numbers are not available in the VEP cache (see A4).

Interpretation (provisional): variants that would receive PP3 evidence but are common in some
population are concentrated in non-European and founder populations; frequency-based criteria
(BA1/BS1) protect against these only when the population is represented in reference databases,
which motivates the future local-cohort analysis (Aim 5).

## R8. Computational evidence for variants of uncertain significance (Aim 7; amendment A5)

ClinVar 2026-09-23: 2,260,073 SNV VUS (≥ 1 star), of which 1,208,461 are MANE missense in disease
genes (1,208,430 with AF < 0.01). Script `23_vus_yield.R`. Percent of VUS per evidence level
(ACMG/AMP points in brackets):

| Evidence | CADD v1.7 (main) | AVI (main) | REVEL | BayesDel | AlphaMissense |
|---|---|---|---|---|---|
| PP3 Strong [+4] | — | — | 4.2 | 3.0 | 4.3 |
| PP3 Moderate [+2] | 12.9 | 16.3 | 7.9 | 8.6 | 10.5 |
| PP3 Supporting [+1] | 14.4 | 7.3 | 6.5 | 7.8 | 7.6 |
| Indeterminate [0] | 26.8 | 19.8 | 18.6 | 22.4 | 23.1 |
| BP4 Supporting [−1] | 16.4 | 19.4 | 14.2 | 13.8 | 15.4 |
| BP4 Moderate [−2] | 21.5 | 28.9 | 31.2 | 31.9 | 28.4 |
| BP4 Strong [−4] | 8.0 | 8.3 | 17.4 | 12.4 | 10.6 |

(REVEL, BayesDel and AlphaMissense with the Aim 6 intersection thresholds; CADD and AVI with
main-analysis thresholds; their intersection-threshold results differ by ≤ 0.4 points.)

CADD v1.7 vs AVI: the two predictors gave evidence in opposite directions for 46,377 VUS (3.8%):
36,148 with PP3 by CADD and BP4 by AVI, and 10,229 the reverse (cross-table
`calibration/results/v17_avi/vus/table_vus_cadd_vs_avi.tsv`).

Interpretation (provisional): most missense VUS would receive benign-side computational
evidence, and roughly one in five to one in four would receive pathogenic-side evidence. These
proportions depend on the thresholds; the less stringent BP4 thresholds estimated with current
data (R6) would increase BP4 evidence relative to the published thresholds, which must be
considered before any clinical interpretation. Computational evidence alone does not reclassify a
variant.

## R9. Threshold drift across ClinVar releases and circularity (amendment A6)

Reference sets (variants scored by all five predictors; plan-v1 criteria): ClinVar 2019-12-02
17,330 P/LP / 11,379 B/LB; 2022-11-19 32,868 / 32,044; 2024-01-27 40,282 / 49,090; 2026-09-23
59,034 / 75,306; within 2026, pre-2022 (present in the 2022-11-19 set) 29,027 / 23,923 and
post-2022 (absent from it) 30,007 / 51,383. One seed, B = 10,000 per set and predictor (scripts
`24a`–`24c`; tables `calibration/results/v17_avi/release_series/`).

**1. The discrepancy with published thresholds is not temporal.** Calibrated on the same ClinVar
release as the original study (2019-12-02), our reference-set construction already gave REVEL BP4
thresholds of ≤ 0.400 (Supporting), ≤ 0.270 (Moderate), ≤ 0.045 (Strong) and ≤ 0.001 (Very
Strong), versus published ≤ 0.290, ≤ 0.183, ≤ 0.016 and ≤ 0.003, and BayesDel BP4 ≤ −0.12 and
≤ −0.22 versus ≤ −0.18 and ≤ −0.36; PP3 thresholds on the 2019 release were close to published
(REVEL 0.682 / 0.780 / 0.928 vs 0.644 / 0.773 / 0.932). The difference is therefore attributable
mainly to reference-set construction — in particular the original exclusion of HGMD variants and of
predictor-training variants, which we could not apply, and other filters — rather than to changes
over time. Absolute thresholds for REVEL, BayesDel and AlphaMissense from this study should not be
used interchangeably with the published ClinGen thresholds.

**2. Drift and circularity (post hoc common scale).** To compare predictors with different score
scales, each threshold was expressed as the percentage of the 2026 gnomAD intersection variants
(1,888,558) meeting it (higher = less stringent); this metric was added after A6 was registered
(script `24d`, post hoc). Change in PP3 (percentage points of gnomAD variants):

| Predictor (ClinGen recommendation) | Supporting, 2026 − 2019 | Supporting, post − pre | Moderate, 2026 − 2019 | Moderate, post − pre | Strong, 2026 − 2019 | Strong, post − pre |
|---|---|---|---|---|---|---|
| REVEL (2022) | +1.95 | +1.44 | +1.10 | +0.90 | +1.02 | +0.79 |
| BayesDel (2022) | +2.35 | +1.21 | +1.57 | +0.40 | +0.52 | +0.28 |
| CADD (2022) | +1.44 | +3.59 | +1.78 | +3.83 | — | — |
| AlphaMissense (2025) | −0.04 | +1.15 | +0.78 | −0.24 | −0.34 | −0.20 |
| AVI (2026; negative control) | +0.72 | +1.24 | −0.12 | −0.40 | — | — |

At the Supporting level all predictors, including the negative control, relaxed by a similar
amount between pre- and post-2022 variants (+1.2 to +1.4 points; CADD +3.6), which points to a
change in the composition of ClinVar. At Moderate and Strong, and over the 2019–2026 series, PP3
relaxation was present for the predictors recommended in 2022 (REVEL, BayesDel, CADD; +0.3 to +3.8
points) and absent or smaller for AlphaMissense and AVI (−0.4 to +0.8), a pattern consistent with a
modest circularity on the pathogenic side. BP4 thresholds did not show this pattern: post-2022
variants gave more stringent BP4 thresholds for most predictors, including AVI.
Limitations: one seed per set without bootstrap intervals for the differences; the pre/post split
confounds time of classification with other differences between the two groups of variants.

## R10. HGMD exclusion (amendment A7)

HGMD Public release 2020.4 in the Ensembl VEP 116 cache: 225,016 records, 187,581 unique
accessions at 170,107 unique positions (CM missense/nonsense 117,486 accessions; CD 30,860; CS
18,547; CI 12,769; others < 4,000 each); 209,566 single-nucleotide records (159,235 positions),
used for the position-level exclusion (scripts `25a`, `25b`). For reference (author), HGMD
Professional 2026.4 contains 584,565 variants and the current public release 354,396.

Reference sets after exclusion (intersection of the five predictors): ClinVar 2019-12-02 17,330 →
5,395 P/LP and 11,379 → 9,681 B/LB (the original study had 5,704 P/LP and 6,059 B/LB); 2026-09-23
59,034 → 26,254 P/LP and 75,306 → 72,077 B/LB; gnomAD sets reduced by < 2%.

With the 2019 release and HGMD positions excluded, the published thresholds were reproduced
closely: REVEL BP4 ≤ 0.288 / 0.188 / 0.014 / 0.008 (published 0.290 / 0.183 / 0.016 / 0.003), PP3
≥ 0.641 / 0.759 / 0.921 (published 0.644 / 0.773 / 0.932); BayesDel BP4 ≤ −0.18 / −0.36
(published −0.18 / −0.36), PP3 ≥ 0.12 / 0.23 / 0.50 (published 0.13 / 0.27 / 0.50);
AlphaMissense moved towards but did not reach the values of Bergquist et al. (BP4 ≤ 0.187 / 0.115
/ 0.018 vs 0.169 / 0.099 / 0.070; PP3 ≥ 0.743 / 0.889 / 0.985 vs 0.792 / 0.906 / 0.972). The
discrepancy reported in R6/R9 is therefore explained mainly by the inclusion of HGMD-listed
variants in our reference sets.

With the 2026 release and HGMD excluded: REVEL PP3 ≥ 0.614 / 0.734 / 0.897 and BP4 ≤ 0.339 /
0.206 / 0.036; BayesDel PP3 ≥ 0.11 / 0.25 / 0.49 and BP4 ≤ −0.17 / −0.34 / −0.58; AlphaMissense
PP3 ≥ 0.694 / 0.879 / 0.990 and BP4 ≤ 0.242 / 0.138 / 0.069. PP3 thresholds of REVEL were more
lenient than published even without HGMD, in line with the PP3 drift described in R9.

CADD v1.7 and AVI were robust to the exclusion at the Supporting and Moderate levels (2026: CADD
PP3 25.6 / 28.0 unchanged, BP4 22.5 / 19.43 vs 22.7 / 19.59; AVI PP3 25.30 / 26.53 vs 25.19 /
26.44, BP4 23.02 / 20.13 vs 22.98 / 20.61); their BP4 Strong thresholds, which lie in sparse score
regions, changed more (CADD 8.38 → 2.76; AVI 9.28 → 8.52).

## R11. Calibration by gene class (Aim 1; amendment A8)

Single-gene 2026 intersection sets (59,014 P/LP, 75,303 B/LB); 296 reference genes lacked a
LOEUF decile. All seven classes met the minimum size (script `26_gene_class_stratification.R`;
Figure 2; tables `calibration/results/v17_avi/gene_classes/`):

| Axis | Class | Genes | P/LP | B/LB | gnomAD |
|---|---|---|---|---|---|
| Constraint | LOEUF deciles 0–2 | 1,829 | 28,457 | 42,032 | 899,214 |
| Constraint | LOEUF deciles 3–6 | 1,973 | 19,158 | 22,905 | 750,992 |
| Constraint | LOEUF deciles 7–9 | 624 | 4,838 | 4,592 | 170,337 |
| Inheritance | AD only | 1,202 | 23,713 | 32,326 | 432,895 |
| Inheritance | AR only | 2,227 | 17,687 | 23,497 | 908,247 |
| Mechanism | Loss of function | 1,800 | 33,633 | 38,451 | 707,521 |
| Mechanism | Gain of function / dominant negative | 158 | 4,992 | 4,574 | 50,378 |

PP3 Supporting thresholds by class (score ≥):

| Predictor | LOEUF 0–2 | LOEUF 3–6 | LOEUF 7–9 | AD only | AR only | LoF | GoF/DN |
|---|---|---|---|---|---|---|---|
| REVEL | 0.663 | 0.631 | 0.581 | 0.666 | 0.584 | 0.667 | 0.745 |
| BayesDel | 0.16 | 0.11 | 0.10 | 0.16 | 0.08 | 0.16 | 0.23 |
| AlphaMissense | 0.788 | 0.515 | 0.493 | 0.829 | 0.496 | 0.695 | 0.775 |
| AVI | 25.91 | 24.60 | 24.20 | 26.07 | 24.44 | 25.30 | 26.21 |
| CADD v1.7 | 26.2 | 25.1 | 24.4 | 26.5 | 25.1 | 25.8 | 26.7 |

BP4 Moderate thresholds by class (score ≤): REVEL 0.272 / 0.251 / 0.203 / 0.279 / 0.244 / 0.273 /
0.261; BayesDel −0.23 / −0.29 / −0.36 / −0.21 / −0.30 / −0.22 / −0.29; AlphaMissense 0.198 / 0.120
/ 0.109 / 0.202 / 0.116 / 0.140 / 0.135; AVI 22.59 / 19.22 / 14.23 / 21.93 / 19.17 / 20.72 / 19.56;
CADD 21.2 / 18.34 / 14.45 / 21.1 / 18.27 / 20.1 / 17.06 (same class order).

The direction was consistent across all five predictors: genes that are more constrained,
associated with dominant inheritance, or with gain-of-function/dominant-negative mechanisms
required higher scores for PP3 and allowed BP4 at higher scores; tolerant and recessive genes
required lower scores for PP3 and much lower scores for BP4. Bootstrap 95% intervals did not
overlap in 8–13 of 15–18 constraint comparisons, 4–5 of 5–6 inheritance comparisons and 1–3 of 3–6
mechanism comparisons per predictor. Using class-specific instead of global thresholds changed the
evidence level of 3–8% (LoF genes) to 26–44% (tolerant genes) of the reference variants, and of
17–54% in GoF/DN genes, depending on the predictor (`reclassification.tsv`).

Interpretation (provisional): a single genome-wide threshold does not fit all gene classes, for
genome-wide (CADD, AVI) and missense-specialised predictors alike; the differences are large
enough to change the evidence assigned to a substantial fraction of variants in some classes.
Limitation: the same prior (0.0441) was used for every class, although the prior probability of
pathogenicity likely differs between classes; class-specific priors would shift the thresholds.

## R12. Genome-wide catalogue (Aim 4; amendment A9)

Full scans: both files contain the same 8,812,917,339 SNVs; 24,926,994 SNVs have CADD v1.7 PHRED
≥ 25.6 and 26,743,913 have AVI Phred ≥ 25.18 (scripts `02`, `28`). Their union, 32,539,174 SNVs,
forms the catalogue (script `29`; 24 per-chromosome files, 0.45 GB, kept outside the repository;
296.6 min): 24,111,910 (74.1%) MANE Select missense (calibrated context), 4,446,346 other MANE
coding SNVs, 3,928,565 in gene bodies outside MANE CDS, and 52,353 intergenic (0.16%). Both
predictors assigned PP3 Moderate to 11,596,538 SNVs; 15.0 M SNVs reached PP3 with only one of them.
The per-gene summary covers 18,842 genes (median 919 catalogued missense SNVs per gene).
**Correction (2026-09-27):** A9 defines "common" by the gnomAD group AF, but script 29 counted
common variants with the overall AF (max of exome and genome AF); the counts were recomputed
(script `29b_gene_summary_group_af.R`; both columns kept in `gene_summary.tsv`). PP3-level missense
SNVs common in some gnomAD group: 8,879 with group AF > 1% and 1,903 with group AF > 5% (2,346 with
overall AF > 1%, the figure implied by the first version). Genes with most such variants (group AF
> 1%): TTN (31), PABPC1 (21), KCNJ12 (18), HSPG2 (12), SLC25A5 (12), DNAH7, OBSCN, OTOF, TNXB
(11 each); with group AF > 5%: PABPC1 (13), KCNJ12 (10), SLC25A5 (9). By overall AF (first
version): PABPC1 (13), KCNJ12 (11), SLC25A5 (11), TTN (10), GXYLT1 (8), RASA4 (6). PABPC1, KCNJ12
and SLC25A5 have processed pseudogenes or close paralogs, suggesting that such "common" variants may
partly reflect paralogous sequence variation or mapping artefacts in gnomAD rather than common
pathogenic variation — a caution for users of the catalogue; large genes (TTN, HSPG2, OBSCN, NEB,
dyneins) contribute by size.

## R13. Noncoding and synonymous variants (Aim 9; amendments A10, A11)

ClinVar 2026-09-23 SNVs with ≥ 1 star outside MANE Select CDS: 517,556; after region assignment,
plan-v1 rules (rare, disease genes) and scoring: 1,123,855 noncoding/synonymous reference SNVs
(script `30_noncoding_clinvar_mpra.R`; `calibration/results/v17_avi/noncoding/`). Execution notes:
run 3 was lost when the session ended after > 4 h of sequential lookups; run 4 failed at the MPRA
join (SNVs tested in more than one element); both fixed (cached parallel lookups; deduplicated join).

Reference variants per class (P/LP / B/LB): canonical splice 40,194 / 135; intronic splice region
(±3–8) 2,392 / 68,917; 5′UTR 126 / 3,915; 3′UTR 54 / 10,337; near-intronic (9–50) 635 / 241,278;
deep intronic 810 / 25,431; upstream ≤ 5 kb 169 / 2,238; downstream ≤ 5 kb 104 / 657; intergenic
157 / 2,414; synonymous exon-edge 435 / 10,981; other synonymous 295 / 609,779. Only the intronic
splice region met the A10 criterion for calibration (≥ 1,000 of each).

**Calibration of the intronic splice region (A10; script `30b`; one seed, B = 10,000).** gnomAD set:
175,221 rare SNVs observed in gnomAD at splice-region intronic positions (709,506 candidate
positions in the 3,703 reference genes), absent from ClinVar.

| Predictor | PP3 Sup | PP3 Mod | PP3 Strong | BP4 Sup | BP4 Mod | BP4 Strong |
|---|---|---|---|---|---|---|
| CADD v1.7 (splice region) | ≥ 15.45 | ≥ 17.36 | ≥ 20.10 | ≤ 11.52 | ≤ 6.94 | not reached (31.6% of replicates fail) |
| CADD v1.7 (missense, R1) | ≥ 25.6 | ≥ 28.0 | not reached | ≤ 22.6 | ≤ 19.59 | ≤ 8.34 |
| AVI (splice region) | ≥ 21.30 | ≥ 22.44 | ≥ 24.09 | ≤ 17.50 | ≤ 14.01 | ≤ 0.11 (sparse region) |
| AVI (missense, R1) | ≥ 25.18 | ≥ 26.44 | not reached | ≤ 22.95 | ≤ 20.52 | ≤ 9.15 |

Both predictors reached PP3 strong in the splice region, with thresholds ~5–10 units lower than for
missense variants. Caveat (declared in A10): the missense prior (0.0441) was used; the true prior in
this region likely differs.

**Missense thresholds applied to other classes (interval LRs, 1,000 bootstrap replicates;
`evaluation_missense_thresholds_by_class.tsv`).** PP3 intervals met the required LR in almost every
class where they contained variants (e.g. near-intronic AVI PP3 moderate LR 717; deep intronic CADD
421), but captured few P/LP variants. BP4 intervals, by contrast, failed in almost every noncoding
class, and the BP4 supporting interval had LR > 1 — i.e. was enriched for P/LP variants — in most
classes (near-intronic AVI 15.8, CADD 12.0; splice region CADD 59.9; deep intronic AVI 8.8, CADD 7.5;
other synonymous CADD 204). Pathogenic noncoding variants score lower than pathogenic missense
variants, so missense BP4 thresholds would give false benign evidence outside missense variation.
Canonical splice variants (40,194 P/LP, 135 B/LB) could not be evaluated meaningfully.

**Discrimination (AUROC, P/LP vs B/LB) by class, AVI vs CADD:** splice region 0.991 vs 0.962;
synonymous exon-edge 0.984 vs 0.978; near-intronic 0.956 vs 0.779; deep intronic 0.908 vs 0.834;
3′UTR 0.880 vs 0.837; upstream 0.874 vs 0.805; other synonymous 0.835 vs 0.781; 5′UTR 0.696 vs 0.659;
canonical splice 0.638 vs 0.596; intergenic 0.958 vs 0.944; downstream 0.942 vs 0.961.

**MPRA (Kircher et al. 2019; 39,170 SNVs in 20 elements/conditions, 6,684 functional).** AUROC AVI
0.679 vs CADD 0.565; difference 0.114 [95% CI 0.108–0.120]; AVI higher in 23 of 28 element conditions
with enough variants (e.g. PKLR 0.81–0.83 vs 0.64–0.66, MSMB 0.76 vs 0.54, HBB 0.79 vs 0.64);
both near 0.5 or below for SORT1 and BCL11A.

Interpretation (provisional): outside missense variants, AVI is consistently more informative than
CADD; missense-calibrated thresholds must not be used for BP4 in noncoding classes, while class-
specific calibration (feasible for the splice region) yields much lower PP3 thresholds, including
PP3 strong.

## R17. CADD v1.6 versus v1.7 on identical reference data (Aim 2b; amendment A14)

CADD v1.6 GRCh38 scores were extracted for all 2,093,102 reference SNVs from 6,419 byte ranges of the
official file (9.38 GB of 86.6 GB; 938,158,116 records read; scripts `35a`–`35c`); every ClinVar and
gnomAD reference SNV had a v1.6 score, so the comparison set equals the main set (59,979 P/LP,
77,128 B/LB, 1,955,995 gnomAD). Spearman correlation v1.6–v1.7 on the ClinVar set 0.937. Both versions
calibrated with five seeds, B = 10,000 (script `35d`; `calibration/results/v17_avi/cadd16_vs_17/`).

| Level | CADD v1.6, 2026 data | CADD v1.7, 2026 data (this comparison) | CADD v1.6, published (2019 data) |
|---|---|---|---|
| PP3 Supporting | ≥ 25.2 | ≥ 25.6 | ≥ 25.3 |
| PP3 Moderate | ≥ 27.5 | ≥ 28.0 | ≥ 28.1 |
| BP4 Supporting | ≤ 22.9 | ≤ 22.6 [22.6–22.7] | ≤ 22.7 |
| BP4 Moderate | ≤ 21.1 | ≤ 19.59 | ≤ 17.3 |
| BP4 Strong | ≤ 7.35 [7.34–7.35] | ≤ 8.34 | ≤ 0.15 |

PP3 strong was not reached by either version. The v1.7 thresholds reproduce the main analysis.

- Discrimination: AUROC v1.6 0.904 [0.902–0.906], v1.7 0.902 [0.900–0.903]; paired difference
  v1.7 − v1.6 −0.002 [−0.003 to −0.001]. On rare missense ClinVar variants, CADD v1.7 is not more
  discriminative than v1.6.
- Evidence levels (each version at its own median thresholds): 38.5% of P/LP (23,106) and 33.0% of
  B/LB (25,444) reference variants change level between versions; 90% of P/LP and 88% of B/LB changes
  are to an adjacent level; the direction of evidence is reversed (PP3 ↔ BP4) for 333 P/LP and 535
  B/LB variants.
- Separating version from data: with the same (2026) data, the version change moved thresholds by
  ≤ 0.5 units except BP4 moderate (21.1 → 19.59) and BP4 strong (7.35 → 8.34); with the same version
  (v1.6), the change of reference data from the 2019 study to 2026 moved BP4 moderate from 17.3 to 21.1
  and BP4 strong from 0.15 to 7.35, consistent with the reference-set effects of R9–R10.

Interpretation (provisional): the v1.7 update did not improve missense discrimination, and although
thresholds shifted little, a third of reference variants would receive a different evidence level
under v1.7 than under v1.6, so thresholds must be re-derived for each model version. The larger BP4
differences with the published values are due to the reference data rather than to the version.

## R14. Mexico City Prospective Study frequencies (Aim 8 extension; amendment A12(a))

Data: RGC-MCPS allele-frequency VCFs (release 20230130; ~140,000 adults; frequencies overall and
by local ancestry: AFR, EUR, IMX = Indigenous Mexican). Script `32_mcps_frequencies.R`. Universe:
24,084,172 PP3-level MANE Select missense SNVs of the Aim 4 catalogue (CADD or AVI); 678,887 were
observed in MCPS (AC > 0; source: exome 659,603, array 14,120, genome 5,164). Levels are exclusive
(Supporting = Supporting only). "gnomAD group max" = maximum over gnomAD v4.1 genetic ancestry groups.

| Predictor | PP3 level | SNVs | Observed in MCPS | MCPS AF > 1% | MCPS AF > 5% | MCPS AF > 1% and gnomAD group max ≤ 1% | IMX AF > 1% | IMX AF > 1% and gnomAD group max ≤ 1% |
|---|---|---|---|---|---|---|---|---|
| CADD v1.7 | Supporting | 9,968,836 | 287,521 | 676 | 309 | 36 | 505 | 73 |
| CADD v1.7 | Moderate | 8,566,291 | 280,389 | 349 | 138 | 30 | 274 | 69 |
| AVI | Supporting | 5,685,222 | 144,585 | 247 | 119 | 15 | 189 | 36 |
| AVI | Moderate | 13,135,064 | 267,222 | 238 | 97 | 18 | 178 | 38 |

- Variants observed in MCPS but absent from gnomAD v4.1: 59,048–126,645 per predictor/level.
- MCPS-only frequency conflicts (PP3 by either predictor, MCPS AF > 1%, every gnomAD group ≤ 1%):
  74 SNVs (`calibration/results/v17_avi/mcps/mcps_only_frequency_conflicts.tsv`). Most have their
  highest gnomAD frequency in the admixed American group (0.6–1%) and roughly 1.5-fold higher
  frequency in the Indigenous Mexican ancestry component of MCPS (2–3%); several are classified
  benign/likely benign in ClinVar. Using the Indigenous Mexican ancestry frequency, 36–73 PP3-level
  SNVs per predictor/level exceed 1% while being ≤ 1% in every gnomAD group.
- Quality check, ClinVar P/LP (≥ 1 star) catalogue SNVs: 6,009 of 52,020 observed in MCPS; only 2
  with AF > 1% (1 > 5%). The latter, ALDH7A1 p.Asn195Ser (ClinVar Pathogenic), has MCPS AF 0.28 with
  similar frequency in all three ancestry components but gnomAD AF 2 × 10⁻⁵, and was called only
  from the genotyping array: most likely an array genotyping artefact.
- Post hoc (not pre-specified): source of the 74 conflicting variants — exome 54, array 19,
  genome 1; array-only calls should be interpreted with caution.
- Concordance with gnomAD v4.1 admixed American AF (252,190 SNVs with AF > 0 in both): Spearman
  ρ = 0.52 (overall MCPS AF); ρ = 0.10 for the Indigenous Mexican component.
- Execution note: a first run compared MCPS with the gnomAD *overall* AF instead of the pre-specified
  group maximum (354 apparent conflicts); corrected before interpretation (log
  `32_mcps_frequencies_run1_overallAF_bug.log`).

Interpretation (provisional): a small number of variants that CADD or AVI place at PP3 levels are
common in the Mexican population, mostly Indigenous-ancestry alleles that gnomAD captures only
partially; in Mexican patients, frequency criteria based on gnomAD alone would miss them.
Licence: RGC-MCPS summary data used under its licence agreement (citation of Ziyatdinov et al. 2023
and acknowledgement of RGC, University of Oxford, UNAM and INMEGEN required); only aggregate counts
and variant identifiers are reported.

## R15. TOPMed BRAVO frequencies (amendment A12(b))

Data: BRAVO Freeze 10 public VCFs (overall AC/AF only; GRCh38); FILTER = PASS sites. Script
`34_bravo_frequencies.R`. Of the 24,084,172 PP3-level MANE missense SNVs, 1,421,923 are at a BRAVO
site (1,405,918 PASS).

| Predictor | PP3 level | Observed (PASS) | AF > 1% | AF > 5% | AF > 1% and gnomAD group max ≤ 1% | With homozygotes |
|---|---|---|---|---|---|---|
| CADD v1.7 | Supporting | 607,740 | 967 | 297 | 0 | 11,285 |
| CADD v1.7 | Moderate | 538,224 | 481 | 124 | 0 | 6,734 |
| AVI | Supporting | 318,672 | 353 | 111 | 0 | 5,906 |
| AVI | Moderate | 574,209 | 329 | 85 | 0 | 6,686 |

ClinVar P/LP (≥ 1 star): 8,856 observed, 4 with AF > 1%, none > 5%. Spearman correlation of BRAVO
AF with gnomAD overall AF: 0.67 (1,125,239 SNVs). No PP3-level SNV was common in BRAVO while ≤ 1% in
every gnomAD group — expected, because gnomAD v4.1 genomes include many TOPMed samples (the two
resources are not independent; A12). BRAVO therefore adds no population-specific frequency
information beyond gnomAD groups for these variants, in contrast with MCPS (R14).

## R16. Validation with multiplexed assays of variant effect (Aim 3; amendment A13)

MaveDB (accessed 2026-09-27): 89 human score sets with investigator-provided (or other) normal /
abnormal functional classes passed all checks (38 mapped at nucleotide level through MaveDB's GRCh38
genomic mappings; 51 protein-level sets mapped through the RefSeq protein change, with 99.2–100% of
their reference residues matching MANE Select). After excluding 5,259 SNVs classified discordantly
by different score sets: 93,154 MANE Select missense SNVs in 39 genes (18,936 functionally abnormal,
74,218 normal), all scored by the five predictors (scripts `33a`–`33c`; tables
`calibration/results/v17_avi/mave/`).

Likelihood ratio LR = P(interval | abnormal) / P(interval | normal) per evidence interval of the
calibrated thresholds (95% bootstrap CI, B = 2,000):

| Interval | CADD v1.7 | AVI | REVEL | BayesDel | AlphaMissense |
|---|---|---|---|---|---|
| PP3 Strong | — | — | 4.35 [4.21–4.49] | 6.09 [5.86–6.34] | 4.29 [4.15–4.43] |
| PP3 Moderate | 2.50 [2.44–2.56] | 2.86 [2.81–2.92] | 1.98 [1.92–2.05] | 1.93 [1.87–2.00] | 2.19 [2.13–2.25] |
| PP3 Supporting | 1.64 [1.59–1.69] | 1.46 [1.39–1.52] | 1.19 [1.14–1.24] | 1.53 [1.47–1.59] | 1.18 [1.13–1.24] |
| Indeterminate | 0.89 [0.87–0.92] | 0.72 [0.69–0.75] | 0.68 [0.66–0.70] | 0.62 [0.60–0.63] | 0.64 [0.61–0.66] |
| BP4 Supporting | 0.41 [0.39–0.44] | 0.39 [0.37–0.42] | 0.45 [0.42–0.48] | 0.41 [0.38–0.43] | 0.34 [0.32–0.36] |
| BP4 Moderate | 0.23 [0.22–0.25] | 0.22 [0.21–0.23] | 0.30 [0.28–0.32] | 0.25 [0.23–0.26] | 0.20 [0.18–0.21] |
| BP4 Strong | 0.13 [0.11–0.15] | 0.12 [0.10–0.14] | 0.10 [0.08–0.11] | 0.08 [0.07–0.10] | 0.15 [0.13–0.17] |

- The pre-specified expectation was met for every predictor: LR > 1 for every PP3 interval and < 1
  for every BP4 interval, with LRs ordered by evidence strength (monotonic) in all five predictors.
  Proportion of abnormal variants: 39–42% at PP3 Moderate (CADD/AVI) and 52–61% at PP3 Strong
  (missense-specialised tools) versus 2–4% at BP4 Strong (baseline 20.3%).
- Pooled AUROC (abnormal vs normal): AlphaMissense 0.782 [0.778–0.785], BayesDel 0.766
  [0.763–0.770], AVI 0.764 [0.760–0.767], REVEL 0.757 [0.754–0.761], CADD v1.7 0.725
  [0.721–0.729]. AVI − CADD paired difference +0.038 [0.036–0.041]; AVI's AUROC exceeded CADD's in
  32 of 36 genes with ≥ 30 abnormal and ≥ 30 normal SNVs. On this ClinVar-independent benchmark AVI
  performed on par with REVEL and BayesDel, unlike on ClinVar (R6), where REVEL and BayesDel led by
  0.017–0.020.
- Sensitivity (A13): excluding the 3,450 SNVs of the ClinVar calibration set left the results
  essentially unchanged (AUROC CADD 0.719, AVI 0.754; same LR pattern).
- Per-gene AUROC varied widely (e.g. TP53 0.91–0.98, BRCA1 0.84–0.91; FGFR4, TPK1 and SCN5A ≤ 0.59
  for all predictors), reflecting assay-specific phenotypes (for example, gain-of-function or
  assay-limited readouts).
- Caveat (pre-specified): functional abnormality is not pathogenicity, so these LRs are not
  comparable with the ClinGen LR targets; they confirm the direction and ordering of the evidence
  intervals on data independent of ClinVar.

Interpretation (provisional): the ClinVar-calibrated intervals of all five predictors behave as
expected on independent functional data; the ranking of predictors differs between ClinVar and
MAVE benchmarks (AlphaMissense first and AVI comparable with REVEL/BayesDel on MAVE), consistent
with ClinVar-based comparisons favouring predictors whose evidence has contributed to ClinVar
classifications or whose training overlaps it.

## R5. Temporal hold-out validation (Aim 2c; amendment A2)

Consistency check: the 2026 reference sets rebuilt from the all-MANE gnomAD table were identical
to those of the main analysis (137,107 ClinVar and 1,955,995 gnomAD variants; keys, labels,
review status, AF and scores). Calibration set: ClinVar 2024-01-27, 40,870 P/LP and 50,278 B/LB
(4,171 genes); gnomAD set 2,009,577. Test set: 22,626 P/LP and 36,294 B/LB variants of the 2026
main set absent from the 2024 set; 11 variants whose label changed between releases were excluded
(6 P→B, 5 B→P). Script `18_temporal_holdout.R`; Figure S3.

Interval-based positive likelihood ratios on the test set, thresholds calibrated on 2024 data
(95% bootstrap intervals, 1,000 replicates):

| Interval | Required LR | CADD v1.7 LR [95% CI] | AVI LR [95% CI] |
|---|---|---|---|
| PP3 Moderate | ≥ 5.79 | 11.2 [10.5–11.8] | 18.1 [17.3–19.1] |
| PP3 Supporting | ≥ 2.41 | 5.44 [5.19–5.73] | 4.49 [4.18–4.81] |
| Indeterminate | — | 1.43 [1.39–1.48] | 1.16 [1.11–1.21] |
| BP4 Supporting | ≤ 0.416 | 0.304 [0.286–0.323] | 0.264 [0.246–0.281] |
| BP4 Moderate | ≤ 0.173 | 0.060 [0.055–0.065] | 0.073 [0.067–0.078] |
| BP4 Strong | ≤ 0.0298 | 0.010 [0.007–0.014] | 0.014 [0.012–0.018] |

For both predictors, every evidence interval met the required likelihood ratio, including at
the conservative end of its 95% interval. Test-set counts per interval (P/LP / B/LB): CADD PP3
Moderate 8,588 / 1,235; AVI PP3 Moderate 14,755 / 1,305 (full table:
`calibration/results/v17_avi/temporal_holdout/interval_LR_test_set.tsv`).

Threshold drift between releases (seed 20260925; 2024 → 2026): CADD PP3 Supporting 26.0 → 25.6,
Moderate 28.6 → 28.0; BP4 Supporting 22.8 → 22.7, Moderate 19.50 → 19.59, Strong 1.05 → 8.34.
AVI PP3 Supporting 25.32 → 25.18, Moderate 26.59 → 26.44; BP4 Supporting 23.00 → 22.95,
Moderate 20.60 → 20.52, Strong 8.94 → 9.15. All thresholds moved by ≤ 0.6 units except CADD
BP4 Strong, which lies in a sparse score region (compare sensitivity analysis S1, 1.74).

## R18. VUS of January 2024 resolved by September 2026 (Aim 7b; amendment A15)

Script `42_vus_resolved.R` (log `logs/42_vus_resolved.log`); Table 2 (manuscript v5), Figure 3B,
Tables S18–S19. ClinVar 2024-01-27: 1,081,487 SNV VUS with >= 1 star; 657,215 MANE missense in
disease genes (657,184 with AF_max < 0.01). By 2026-09-23: 1,790 resolved to P/LP and 3,989 to B/LB
(0.9% of the rare VUS); 651,405 not resolved.

Primary analysis (r2024 thresholds, AF < 0.01, resolved variants):

| Predictor | P/LP with PP3 % (mean points) | B/LB with BP4 % (mean points) | Final direction / opposite % |
|---|---|---|---|
| CADD v1.7 | 65.5 (+0.93) | 58.3 (−0.82) | 60.5 / 13.6 |
| AVI | 77.1 (+1.29) | 70.5 (−1.21) | 72.5 / 9.8 |
| REVEL | 77.1 (+2.04) | 75.9 (−1.54) | 76.3 / 7.3 |
| BayesDel | 76.4 (+1.90) | 73.4 (−1.57) | 74.3 / 8.1 |
| AlphaMissense | 76.6 (+1.88) | 70.7 (−1.40) | 72.5 / 8.2 |

Interval likelihood ratios (resolved P/LP vs B/LB): supporting and moderate intervals meet the
targets at the conservative bound for all predictors except CADD PP3 moderate (5.54, 95% CI
4.90–6.35; 5.79 required). AVI PP3 moderate 9.71 (8.64–11.0). BP4 strong not confirmed for CADD,
AVI, REVEL and BayesDel (2–7 P/LP variants per interval); AlphaMissense PP3 strong 32.0 (23.8–46.6)
short of 33.53; REVEL (46.9) and BayesDel (50.9) PP3 strong met. 2026 thresholds give the same pattern.
Prespecified examples (expert-panel records): three MYOC variants (likely benign) with BP4 from all
five predictors; RPE65 p.His527Arg (pathogenic) with PP3 from all five; LDLR p.Gly461Cys and PIK3CD
p.Asn334Lys (likely pathogenic) with BP4 from AVI and CADD (PIK3CD: gain-of-function gene).
Interpretation: calibrated evidence points in the final direction for about three quarters of resolved
VUS (AVI and missense predictors), contributing one to two points on average; AVI (negative control)
is close to the missense predictors that 2026 classifications could cite.

## R19. Current ClinVar missense VUS with concordant calibrated evidence (Aim 7c; amendment A16)

Script `46_vus_reclassification_candidates.R` (log `executed_pipeline/logs/46_vus_reclassification_candidates.log`),
2026-10-29. 1,208,461 VUS (MANE Select missense, >= 1 star, disease genes; ClinVar 2026-09-23), 1,208,430 with
AF_max < 0.01 (counts below).

| Concordant evidence (REVEL, BayesDel, AlphaMissense published thresholds + AVI) | All rare VUS | Expert panel/guideline |
|---|---|---|
| PP3 strong (+4) | 10,553 (0.87%) | 186 |
| PP3 moderate (+2) | 46,931 (3.88%) | 341 |
| PP3 supporting (+1) | 48,013 (3.97%) | 280 |
| BP4 supporting (−1) | 210,699 (17.44%) | 175 |
| BP4 moderate (−2) | 167,479 (13.86%) | 112 |
| None | 724,755 (59.97%) | 1,854 |
| Total | 1,208,430 | 2,948 |

Concordant VUS: 483,690 (483,675 rare); expert-panel concordant 1,094 in 97 genes; 23 genes with >= 10
(856 variants; Table 4): RUNX1 217 (139 PP3, 78 BP4), RYR1 101 (100 PP3), HNF1A 58, LDLR 58, HNF4A 41, PAH 40
(all PP3), MYOC 39 (30 BP4), GCK 38 (all PP3), ITGA2B 30 (19 BP4). In the gene-class annotation RUNX1, RYR1,
HNF1A, HNF4A and CDH1 are LOEUF 0–2.
Transferred thresholds (CADD v1.7 scores, published v1.6 thresholds) vs calibrated CADD v1.7: 216,454 (17.91%)
change level (92,776 BP4 supporting -> moderate; 79,467 BP4 moderate -> strong; 30,239 PP3 supporting ->
indeterminate; 8,900 BP4 supporting -> indeterminate; 5,072 PP3 supporting -> moderate), none change direction.
Transferred CADD v1.6 vs calibrated AVI: 55,003 (4.55%) opposite directions (43,919 PP3 by CADD and BP4 by AVI;
11,084 the reverse).
Interpretation: calibrated predictors agree on evidence for 40% of current missense VUS; concordant PP3 strong
completes a likely pathogenic classification with two further points. Expert panels may already apply their
own computational criteria, so the list identifies candidates for reassessment, not reclassifications.
Transferring v1.6 thresholds to v1.7 mainly weakens benign evidence.

## R20. Main thresholds outside missense variants on functional data (Aim 9 extension; amendment A17)

Script `65_functional_noncoding_main_thresholds.R` (log `65_functional_noncoding.log`), 2026-11-01; outputs in
`calibration/results/v17_avi/functional_noncoding/`. Main thresholds (Table 1) and, in the intronic splice region,
the splice-region thresholds of R13. LR = P(interval | abnormal or functional) / P(interval | normal or
non-functional); 95% CI from 1,000 bootstrap resamples (seed 20261101). Functional abnormality is not pathogenicity,
so comparisons with the ClinGen LRs are descriptive.

**Data.** SGE: 60,382 SNVs of the 12 saturation genome editing sets (all scored by CADD v1.7 and AVI); evaluable
classes (>= 30 abnormal and >= 30 normal): missense 5,436 / 28,122, other synonymous 287 / 9,628, stop gained 2,048 /
108, synonymous exon edge 56 / 337, canonical splice 1,071 / 117, intronic splice region 619 / 2,610, near intronic
386 / 7,277. MPRA: 39,170 SNV-element records (6,684 functional), the set of R13.

**Missense (control).** All five evaluable intervals of each predictor have the expected direction (CADD PP3 LR
2.04–3.29, BP4 0.12–0.38; AVI PP3 1.54–3.35, BP4 0.09–0.33), as in R16.

**Missense thresholds outside missense variants.** The BP4 Supporting interval has LR > 1, with the 95% CI above 1,
in every noncoding or synonymous class where it is evaluable, except AVI in the synonymous exon edge (0.65) and
intronic splice region (1.11, CI including 1) classes:

| Class | CADD BP4 Sup | CADD BP4 Mod | CADD BP4 Strong | AVI BP4 Sup | AVI BP4 Mod | AVI BP4 Strong |
|---|---|---|---|---|---|---|
| Other synonymous | 22.4 (10.8–41.3) | 1.25 (1.10–1.41) | 0.76 (0.66–0.87) | 5.39 (3.63–7.51) | 0.91 (0.81–1.01) | 0.65 (0.52–0.79) |
| Synonymous exon edge | 10.7 (6.4–20.7) | 0.48 (0.31–0.66) | 0.22 (0.04–0.48) | 0.65 (0.14–1.44) | 0.21 (0.08–0.36) | 0.33 (0.07–0.86) |
| Intronic splice region | 3.15 (2.55–3.85) | 1.12 (1.02–1.22) | 0.17 (0.13–0.23) | 1.11 (0.82–1.46) | 0.36 (0.30–0.43) | 0.12 (0.07–0.16) |
| Near intronic | 3.26 (1.94–4.87) | 1.28 (1.15–1.41) | 0.62 (0.53–0.70) | 5.50 (3.86–7.80) | 1.14 (0.99–1.30) | 0.42 (0.35–0.49) |
| MPRA (all elements) | 1.69 (1.53–1.88) | 1.28 (1.24–1.31) | 0.75 (0.73–0.78) | 1.83 (1.71–1.97) | 1.36 (1.33–1.40) | 0.45 (0.43–0.47) |

BP4 Strong has LR < 1 in every class,
but no BP4 interval reaches the ClinGen LR outside missense variants. Because noncoding SNVs score low, most
functionally abnormal SNVs receive BP4 under the missense thresholds: near intronic CADD 93.5%, AVI 75.1%; intronic
splice region CADD 82.1%, AVI 29.9%; other synonymous CADD 99.7%, AVI 91.3%; synonymous exon edge CADD 82.1%, AVI
25.0%; MPRA CADD 98.5%, AVI 92.2% (normal SNVs 90–100%). Per MPRA element (intervals with >= 10 records): the AVI BP4
Supporting interval has LR < 1 in 0 of 26 element conditions and BP4 Moderate in 8 of 28; BP4 Strong in 28 of 28;
CADD BP4 Supporting 2 of 15, Moderate 8 of 27, Strong 21 of 28.
The PP3 intervals of the missense thresholds have LR > 1 wherever evaluable outside the canonical splice and stop
gained classes (AVI: intronic splice region 13.3 and 19.5, near intronic 64.9 and 94.3, synonymous exon edge 5.3 and
16.0, MPRA 5.8 and 6.4), but capture few abnormal SNVs (AVI PP3 3–55% of abnormal SNVs in these classes, 2.5% in
the MPRA; CADD at most 4.6%).
Canonical splice and stop gained SNVs are almost all abnormal; their intervals have LR close to 1 (AVI PP3 Moderate
1.21 and 1.02) and are not informative.

**Splice-region thresholds in the intronic splice region.** All intervals of both predictors have the expected
direction: CADD PP3 2.18–3.65, BP4 Supporting 0.46, Moderate 0.13; AVI PP3 1.25–11.0, BP4 Supporting 0.42, Moderate
0.14, Strong 0.06. The share of abnormal SNVs receiving BP4 falls from 82.1% to 14.5% (CADD) and from 29.9% to 14.4%
(AVI); PP3 is given to 71.6% (CADD) and 76.7% (AVI) of abnormal SNVs, and to 24.5% and 14.3% of normal SNVs.

**Discrimination (AUROC, AVI − CADD with 95% CI).** Missense 0.802 vs 0.771 (+0.032, 0.028 to 0.035); other
synonymous 0.657 vs 0.600 (+0.056, 0.032 to 0.080); synonymous exon edge 0.855 vs 0.836 (+0.019, −0.016 to 0.056);
canonical splice 0.637 vs 0.529 (+0.107, 0.041 to 0.172); intronic splice region 0.884 vs 0.788 (+0.095, 0.081 to
0.111); near intronic 0.763 vs 0.665 (+0.098, 0.076 to 0.125); stop gained 0.501 vs 0.581 (−0.080, −0.126 to −0.035);
MPRA 0.679 vs 0.565 (+0.114, 0.109 to 0.120, as in R13).

**Sensitivity, without SNVs of the ClinVar reference sets** (SGE 53,209 SNVs; MPRA 39,039 records): the same pattern.
BP4 Supporting LR > 1 for CADD in all four noncoding or synonymous classes and the MPRA (1.68–17.5) and for AVI in
other synonymous (2.22), near intronic (5.87) and the MPRA (1.84); intronic splice region AVI 1.12 (0.79–1.51).
Splice-region thresholds: all intervals in the expected direction.

**Interpretation.** On data independent of ClinVar, the missense-calibrated BP4 Supporting and Moderate intervals of
CADD v1.7 and AVI do not indicate normal function outside missense variants: they contain more functionally abnormal
than normal SNVs in proportion, while BP4 Strong keeps LR < 1 but well above the ClinGen requirement. This confirms
R13 with functional data, and the class-specific calibration of the intronic splice region gives intervals that
behave as expected for both predictors. AVI has a higher AUROC than CADD in every class except stop
gained SNVs (the interval includes 0 for the synonymous exon edge class).

## R-P2-1. Paper 2: reference-set sizes by mechanism stratum and region class (amendment A-P2-1, first step)

Script `59_paper2_reference_sets.R`; log and tables in `calibration/results/paper2/refsets/`. Counts only; no score
was summarised. Rules of plan v1 (>= 1 star, exact labels, MANE Select class gene among the disease genes of the
release, gnomAD group-maximum AF < 0.01, intergenic excluded). Reference SNVs: 1,227,728 (2026-09-23) and 712,489
(2024-01-27); test set (2026 variants absent from the 2024 set) 538,416; 28 label changes excluded. Without Atlas
attributions and excluded: 76 SNVs in the 2026 set (4 P/LP) and 23 in the 2024 set (0 P/LP), all with a MANE class
gene that is a disease gene while their ClinVar GENEINFO gene is not (the query list used GENEINFO, script 53).

| Stratum | 2026 P/LP | 2026 B/LB | 2024 P/LP | 2024 B/LB | Test P/LP | Test B/LB |
|---|---|---|---|---|---|---|
| Splicing | 46,924 | 199,458 | 28,634 | 110,659 | 19,712 | 91,711 |
| Protein | 126,762 | 64,657 | 83,011 | 42,542 | 47,920 | 29,339 |
| Regulatory | 1,385 | 386,233 | 838 | 220,328 | 617 | 169,272 |
| Mixed | 1,046 | 18,428 | 708 | 11,116 | 421 | 8,879 |
| Conservation only | 531 | 169,556 | 325 | 93,811 | 221 | 76,792 |
| No predicted effect | 83 | 212,665 | 42 | 120,475 | 46 | 93,486 |

Size criterion (>= 1,000 P/LP and >= 1,000 B/LB): met in 2026 by the splicing, protein, regulatory and mixed
strata, and in 2024 by the splicing and protein strata only. Among region classes it is met in both releases only
by intronic splice region and missense (sizes by class in `sizes_by_class.tsv`).
Composition of the regulatory stratum in 2026: of 1,385 P/LP, 840 are coding (558 stop gained, 192 missense,
79 start lost, 11 stop lost), 102 canonical splice, 48 synonymous and 395 noncoding (119 upstream, 107 deep
intronic, 73 5'UTR, 35 3'UTR, 28 intronic splice region, 17 downstream, 16 near intronic)
(`class_by_stratum_2026.tsv`).
Consequences under the plan: the stopping rule of A-P2-1 is not triggered (four strata meet the criterion in
2026). The temporal hold-out, which uses thresholds calibrated on 2024, can evaluate only the splicing and protein
strata and, for region stratification, intronic splice region and missense. Noncoding regulatory variants (395
P/LP) are reported descriptively (secondary question 1), as are the conservation-only (531) and no-predicted-effect
(83) strata.

## R-P2-2. Paper 2: calibrated AVI Phred thresholds by stratum and class (script 60)

Script `60_paper2_calibrate.R`; outputs in `calibration/results/paper2/calibration/` (thresholds, curves and bootstrap
thresholds per analysis, set and seed; `set_sizes.tsv`; `60_calibrate.log`); summary tables
`thresholds_all_calibrations.tsv` and `seed_range_2026.tsv` in `calibration/results/paper2/evaluation/`. 62
calibrations; sets that meet the size criterion: in 2026 the splicing, protein, regulatory and mixed strata, the
intronic splice region and missense classes and all classes pooled; in 2024 the same except the regulatory and mixed
strata. Discounted thresholds (AVI Phred, seed 20260925), PP3 / BP4; a dash means the level is not reached.

| Set | 2026 Supporting | 2026 Moderate | 2026 Strong | 2024 Supporting | 2024 Moderate | 2024 Strong |
|---|---|---|---|---|---|---|
| Splicing stratum | 23.60 / 21.23 | 24.26 / 19.81 | 25.64 / 13.38 | 23.84 / 21.25 | 24.51 / 19.85 | 26.09 / 11.25 |
| Protein stratum | 26.14 / 23.90 | 27.54 / 22.88 | 33.76 / 16.26 | 26.21 / 23.90 | 27.58 / 22.84 | 34.14 / 15.25 |
| Regulatory stratum | 18.75 / 2.11 | 20.39 / 0.11 | – / – | not calibrated | | |
| Mixed stratum | 23.22 / 9.00 | 24.29 / 4.44 | 32.88 / – | not calibrated | | |
| Intronic splice region | 21.30 / 17.49 | 22.45 / 14.09 | 24.06 / 0.11 | 21.34 / 16.91 | 22.76 / 12.91 | 24.51 / – |
| Missense | 25.19 / 22.95 | 26.44 / 20.52 | – / 9.11 | 25.32 / 23.00 | 26.58 / 20.67 | – / 8.92 |
| All classes pooled | 23.40 / 20.65 | 24.26 / 19.02 | 26.37 / 13.96 | 23.56 / 20.69 | 24.46 / 18.99 | 26.71 / 14.01 |

PP3 Very Strong is reached only by the pooled set (2026 56.45; 2024 67.48). Five seeds: the largest range of any
threshold across seeds is 0.32 Phred (mixed stratum, BP4 Moderate, 4.22–4.54); all others are within 0.12. H1 replication: in 2026 the
splicing stratum reaches PP3 up to Strong (25.64) and BP4 up to Strong (13.38).
Prior sensitivity (0.01 and 0.10 instead of 0.0441): thresholds agree with the main ones within 0.07 Phred. This
follows from the method: with c fixed, every evidence level corresponds to a fixed likelihood ratio and the prior
cancels between the local posterior and its threshold, so the analysis as planned cannot show prior dependence.
Population control (2026 P/LP against the rare gnomAD sample): PP3 thresholds are higher and fewer levels are
reached (protein stratum PP3 Supporting 29.66, Moderate 32.50; splicing stratum PP3 Supporting 27.53, Moderate not
reached; pooled set PP3 Supporting 24.91, Moderate 26.94), in line with the dependence on the benign set reported by
Babazade et al.

## R-P2-3. Paper 2: evaluation of H-P on held-out data (script 61; amendment A-P2-3)

Script `61_paper2_evaluate.R`; outputs and log (`61_evaluate.log`) in `calibration/results/paper2/evaluation/`.
Thresholds of the 2024 calibration (seed 20260925, discounted); 1,000 bootstrap resamples (seed 20261101).

**Criterion (a), temporal hold-out** (68,937 P/LP and 469,479 B/LB test SNVs). Evaluable intervals that meet their
required LR (point estimate):

| Strategy | Evaluable | Meet | Proportion | Test SNVs with evidence | Share in intervals that meet |
|---|---|---|---|---|---|
| M, mechanism (with fallback) | 18 | 15 | 0.83 | 522,087 | 0.34 |
| R, region (with fallback) | 16 | 16 | 1.00 | 521,323 | 1.00 |
| G, genome-wide | 6 | 6 | 1.00 | 521,575 | 1.00 |
| X, missense thresholds | 5 | 4 | 0.80 | 527,644 | 0.97 |
| M without fallback | 12 | 12 | 1.00 | 177,536 | 1.00 |
| R without fallback | 10 | 10 | 1.00 | 81,415 | 1.00 |

Differences (95% bootstrap interval): M − R −0.17 (−0.17 to −0.11); M − G −0.17 (−0.17 to −0.02); M − X 0.03
(0.03 to 0.03). With the 95% bound instead of the point LR the counts are the same. All intervals of the splicing
and protein strata meet. The three failures of M are the BP4 intervals of its fallback group (strata without a 2024
calibration: regulatory, mixed, conservation only, no predicted effect; 1,305 test P/LP and 348,429 B/LB), which
receive the pooled thresholds: LR 6.28 (Supporting; required at most 0.42), 1.14 (Moderate; 0.17) and 0.26 (Strong;
0.03). 646 of the 1,305 P/LP of that group (50%) fall in BP4 intervals. The failure of X is BP4 Supporting (LR 0.86).
Under strategies R and G the same variants are pooled with stop-gained, canonical splice and coding variants, and
their BP4 intervals meet. Chromosome split (test variants outside the AVI training chromosomes): same evaluable and
meeting intervals for every strategy, same differences.

**Criterion (b), resolved VUS** (14,513 VUS of 2024 with an exact P/LP or B/LB label in 2026). (Final direction −
opposite direction) / resolved:

| Strategy | Final direction | Opposite direction | Net (95% CI) |
|---|---|---|---|
| M | 12,487 | 653 | 0.815 (0.808–0.824) |
| R | 12,362 | 705 | 0.803 (0.795–0.812) |
| G | 11,813 | 1,107 | 0.738 (0.728–0.747) |
| X | 12,651 | 668 | 0.826 (0.818–0.834) |
| M without fallback | 7,658 | 588 | 0.487 (0.477–0.496) |
| R without fallback | 5,529 | 622 | 0.338 (0.328–0.347) |

Differences: M − R 0.012 (0.008 to 0.017); M − G 0.078 (0.071 to 0.084); M − X −0.010 (−0.015 to −0.005).
Chromosome split (10,304 resolved VUS): M − R 0.013 (0.008 to 0.018); M − G 0.080 (0.072 to 0.088); M − X −0.008
(−0.013 to −0.002). Counts by strength and direction in `b_levels_*.tsv`.

**H-P.** The two criteria disagree. On criterion (a), mechanism stratification is worse than region stratification
and than the genome-wide set, entirely because of BP4 evidence given by pooled thresholds to strata that could not be
calibrated in 2024; the calibrated strata meet every requirement. On criterion (b), mechanism stratification gives a
higher net score than region stratification (+0.012) and the genome-wide set (+0.078), and a lower one than the
missense thresholds applied to every variant (−0.010). Neither result depends on the AVI training chromosomes.

**H2 replication** (strategy M, hold-out LRs by region class; `h2_strategyM_intervals_by_class.tsv`). Classes with
evaluable intervals: deep intronic, missense, near intronic, intronic splice region, stop gained, synonymous exon edge,
synonymous other. BP4 intervals meet in 5 of 21 (missense 2 of 3, intronic splice region 2 of 3, synonymous exon edge
1 of 3; none in deep intronic, near intronic, stop gained or synonymous other). PP3 intervals meet in 17 of 21 (all
except missense PP3 Strong, LR 25.8, and the three of stop gained, where LR is close to 1 with 125 test B/LB). H2 is
not supported.

**Secondary question 2, conservation only** (pooled thresholds; 221 test P/LP, 76,792 B/LB): BP4 Supporting LR 7.19,
Moderate 1.25, Strong 0.30; none meets. No PP3 interval of this stratum is reported in the output.

**Secondary question 1, descriptive bins** (2026 set; `secondary1_descriptive_bins.tsv`). In the regulatory stratum
(1,385 P/LP, 386,233 B/LB) LR rises from 0.25 in [0,10) to 12.5 in [20,25) and 2,634 at 30 or more; LR crosses 1
in the [15,20) bin (1.69). In the 5'UTR, 3'UTR, upstream 5 kb and downstream 5 kb classes LR is 4.8 to 11.9 in [20,25), 5.6 to 306 in
[25,30) and 78 to 3,445 at 30 or more, with wide intervals (54 to 165 P/LP per class).

## R-P2-4. Paper 2: functional validation, clinical outcomes and exploratory hybrid strategy (scripts 62–64; amendment A-P2-4)

Scripts `62_paper2_functional_sets.R`, `63_paper2_functional_evaluate.R` and `64_atlas_splice_tracks.py`; outputs and
logs in `calibration/results/paper2/functional/`. Thresholds of the 2026 calibration (seed 20260925, discounted);
1,000 bootstrap resamples (seed 20261101). Agreement means PP3 for an abnormal (SGE) or functional (MPRA) SNV and BP4
for a normal or non-functional one; net = (agreeing − opposite) / all SNVs.

**Sets.** SGE: 12 MaveDB score sets; 12,620 intronic SNVs placed from exonic anchors (deviation from A-P2-4), all with
a reference allele matching GRCh38 (none in CARD11). Analysis set 60,382 SNVs, of which 60,229 have Atlas attributions
(9,963 abnormal, 50,266 normal). MPRA: 29,444 SNVs (29,432 with attributions; 4,571 functional, 24,861 non-functional);
897 SNVs tested in more than one element with discordant labels were counted as functional (implementation rule, not
stated in A-P2-4). The Atlas query of the 72,828 new SNVs had one transient error, recovered on repair. Strata of the
SGE set: protein 30,375, splicing 10,376, regulatory 8,824, conservation only 3,675, no predicted effect 3,581, mixed
3,398. Strata of the MPRA set: regulatory 21,139, no predicted effect 5,165, conservation only 2,732, splicing 367,
protein 18, mixed 11.

**SGE, net agreement (95% CI)**

| Strategy | Agree | Opposite | Net |
|---|---|---|---|
| M | 35,969 | 8,814 | 0.451 (0.445–0.457) |
| R | 43,892 | 9,528 | 0.571 (0.565–0.576) |
| G | 40,773 | 13,291 | 0.456 (0.450–0.463) |
| X | 44,452 | 9,302 | 0.584 (0.578–0.589) |
| M without fallback | 28,928 | 8,649 | 0.337 (0.331–0.342) |
| R without fallback | 21,905 | 8,565 | 0.221 (0.217–0.227) |

Differences: M − R −0.120 (−0.124 to −0.115); M − G −0.005 (−0.011 to −0.0004); M − X −0.133 (−0.137 to −0.129).
Without SNVs of the 2026 ClinVar reference set (52,222 SNVs): M − R −0.103 (−0.108 to −0.099); M − G 0.027 (0.022 to
0.032); M − X −0.117 (−0.121 to −0.112). Without the AVI training and validation chromosomes (26,075 SNVs): M − R
−0.129 (−0.136 to −0.123); M − G −0.028 (−0.036 to −0.021); M − X −0.142 (−0.148 to −0.136).
By stratum (`sge_by_stratum.tsv`), M is higher than R in the protein stratum (0.402 vs 0.308), similar in the splicing
stratum (0.718 vs 0.712), and much lower in the regulatory (0.030 vs 0.911) and mixed (0.105 vs 0.798) strata. The
regulatory and mixed stratum thresholds give PP3 from AVI Phred 18.75 and 23.22 and BP4 only below 2.11 and 9.00,
whereas the pooled thresholds that R applies to these SNVs give BP4 Strong below 13.96. Most regulatory and mixed SNVs
of the SGE genes are functionally normal, so the stratum thresholds turn pooled BP4 calls into indeterminate or PP3
calls. By gene (`sge_by_gene.tsv`), M is below R in 11 of 12 genes (the exception is SFPQ, 0.138 vs 0.136).
Interval LRs (`sge_intervals.tsv`): M has 24 evaluable intervals, 22 in the expected direction; the two exceptions are
BP4 Supporting (LR 3.64) and BP4 Moderate (1.70) of the fallback group. All 21 intervals of the calibrated strata have
the expected direction (PP3 LR 1.46–37.2; BP4 LR 0.07–0.88). R 17 of 17, G 5 of 6, X 5 of 5. Intervals reaching the
ClinGen LR (descriptive): M 3, R 5, G 1, X 1.

**MPRA, net agreement (95% CI)**

| Strategy | Agree | Opposite | Net |
|---|---|---|---|
| M | 10,434 | 3,618 | 0.232 (0.224–0.239) |
| R | 23,148 | 4,007 | 0.650 (0.642–0.658) |
| G | 23,146 | 4,006 | 0.650 (0.642–0.658) |
| X | 24,534 | 4,333 | 0.686 (0.678–0.694) |
| M without fallback | 3,122 | 3,255 | −0.005 (−0.010 to 0.001) |
| R without fallback | 13 | 15 | −0.0001 (−0.0004 to 0.0003) |

Differences: M − R −0.419 (−0.428 to −0.410); M − G −0.419 (−0.428 to −0.410); M − X −0.455 (−0.465 to −0.445).
Without the 2026 ClinVar reference SNVs the differences are the same to the third decimal; without the AVI training
and validation chromosomes (5,840 SNVs): M − R and M − G −0.171 (−0.192 to −0.151), M − X −0.159 (−0.180 to −0.138).
Under R almost every MPRA SNV falls in the fallback group, so R and G are nearly identical.
Interval LRs (`mpra_intervals.tsv`): the four intervals of the M regulatory stratum all have the expected direction
(PP3 Supporting 1.41, PP3 Moderate 2.08; BP4 Supporting 0.48, BP4 Moderate 0.59). The pooled thresholds, as used by R
and G, have BP4 Supporting and Moderate in the wrong direction (LR 1.75 and 1.57), BP4 Strong 0.65, and PP3 LR 3.8 to
76 on 15 to 251 SNVs. Direction counts: M 6 of 11, R 5 of 7, G 5 of 7, X 3 of 5 (the M failures are in the fallback
group and in the small splicing stratum, 43 functional SNVs). No interval of M reaches the ClinGen LR; R, G and X one
each. The higher net score of R and G comes from BP4 Strong assigned by the pooled thresholds to 17,444 non-functional
SNVs (and 2,080 functional ones); the regulatory-stratum thresholds give BP4 to fewer SNVs and PP3 to 4,347, with LR
1.4 to 2.1. Net agreement therefore rewards BP4 volume in a set where 84% of SNVs are non-functional, while the
direction check favours the regulatory-stratum thresholds. By element (`mpra_by_element.tsv`), M is below R in 21 of
24 element groups.

**Clinical outcomes, current VUS** (1,336,827 rare VUS with at least one star in 2026 disease genes; 3,495 reviewed by
an expert panel or practice guideline). Under M: PP3 292,685 (Supporting 92,349, Moderate 172,342, Strong 27,994),
BP4 695,854 (Supporting 135,516, Moderate 360,407, Strong 199,931), indeterminate 348,288; no PP3 Very Strong. M
without fallback: the same PP3 counts, BP4 670,166, no thresholds for 26,022 VUS. For comparison, R gives PP3 318,837
and BP4 767,043; G gives PP3 499,525 and BP4 548,035 (`vus_current_levels.tsv`; by stratum and class in
`vus_current_levels_by_stratum_class.tsv`). Expert-panel subset under M: PP3 1,337 (Supporting 325, Moderate 847,
Strong 165), BP4 1,306, indeterminate 852.
RNA-study candidates (splicing-stratum VUS with at least PP3 Moderate): 25,847 VUS in 3,571 genes
(`rna_candidates_by_gene.tsv`; most in TTN 226, ATM 209, NF1 169, RYR2 168); 142 in the expert-panel subset (106 PP3
Strong, 36 Moderate; canonical splice 41, intronic splice region 49, missense 37, synonymous exon edge 11, other 4).
Atlas splicing tracks for these 142 (script 64, 0 errors): the maximum absolute SPLICE_SITE_USAGE score has median
0.87 (range 0.02–0.99) and is at least 0.2 for 133 (descriptive cut-off, not prespecified), including 31 of the 37
missense VUS; the maximum SPLICE_JUNCTIONS score has median 6.19 (1.18–10.38). The track with the maximum
SPLICE_SITE_USAGE is most often from hepatocyte (20), testis (10), foreskin keratinocyte and heart left ventricle (6
each). Per-variant track scores are not versioned (AlphaGenome terms).

**Exploratory hybrid strategy (defined after the H-P results; no test).** Criterion (a), 2024 thresholds: H 14 of 17
evaluable intervals meet (test SNVs with evidence 188,168); H_fb 17 of 23 (503,051); M 15 of 18, R 16 of 16, G 6 of 6.
Criterion (b) (resolved VUS rebuilt in script 63; 14,513, identical to script 61): H 0.577 (0.568–0.587), H_fb 0.827
(0.818–0.835), against M 0.815, R 0.803, G 0.738, X 0.826. SGE: H 0.345, H_fb 0.451. MPRA: H −0.005, H_fb 0.232.
With the 2026 thresholds the strata regulatory and mixed are calibrated, so H and H_fb differ from M without fallback
and M only in the conservation-only and no-predicted-effect strata, which are small in the functional sets.

**Summary of A-P2-4.** On both functional data sets mechanism stratification gives a lower net agreement than region
stratification and the missense thresholds, and on SGE a net agreement similar to the genome-wide set; the deficit
comes from the regulatory and mixed strata (SGE) and from fewer BP4 calls in the regulatory stratum (MPRA). In the
protein and splicing strata M is equal to or better than R. The interval LRs of the calibrated strata have the
expected direction in all 21 SGE intervals and in all four MPRA regulatory intervals, which is not true of the pooled
BP4 intervals on the MPRA. The differences between strategies do not depend on the 2026 ClinVar reference SNVs or on
the AVI training and validation chromosomes.
