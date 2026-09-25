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
