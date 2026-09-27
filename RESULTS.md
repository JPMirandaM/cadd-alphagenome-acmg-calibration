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
