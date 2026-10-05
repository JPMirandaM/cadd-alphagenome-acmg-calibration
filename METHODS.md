# Methods (draft for publication)

> Living document. Sections marked **[PENDING]** describe steps not yet executed and must be
> completed with the actual values once run. Items marked **[VERIFY]** were written from memory
> and must be checked against the original source before submission.
> Exact scripts and their SHA-256 hashes are in `executed_pipeline/`; the ordered list of
> commands is in `executed_pipeline/execution_log.ps1`. Figures and their source data are in
> `figures/` (see the Figures section).

## 1. Data source: precomputed CADD scores

We used the precomputed scores of CADD (Combined Annotation Dependent Depletion) version 1.7
for the GRCh38 assembly (Schubach et al., 2024), which score every possible single-nucleotide
variant (SNV) of the reference genome. The file `whole_genome_SNVs.tsv.gz` (87,473,403,655 bytes;
server modification date 2024-01-29; ETag `"145dd23707-61014c1cb9940"`) and its tabix index were
obtained from `https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/`. For each SNV the
file provides chromosome, 1-based position, reference allele, alternative allele, raw score
(*RawScore*) and scaled score (*PHRED*).

The file was downloaded with curl 8.21.0 as 8 byte ranges retrieved in parallel. Only HTTP 206
(partial content) responses were accepted, the segments were concatenated in order, and
integrity was checked by requiring every segment and the final file to have exactly the expected
size, and by comparing MD5 checksums with those published by the CADD developers: both the score
file (`88577a55f1cd519d44e0f415ba248eb9`) and its index (`347df8fac17ea374c4598f4f44c7ce8b`)
matched (download completed 2026-09-25; SHA-256 of the score file:
`11582f442d3db5bc529a3786e9ddfe27860642c2447d64e468387ad9c050645a`). The full scan in section 2
read 8,812,917,339 SNV records.

## 2. Variant selection

The CADD PHRED score is defined as −10·log₁₀(r/N), where r is the rank of the variant among the
N possible SNVs of the genome. Thus PHRED ≥ 10, ≥ 20 and ≥ 30 correspond to the top 10%, 1% and
0.1% of variants. The selection threshold was derived from our own calibration of CADD v1.7
against current ClinVar and gnomAD data (section 3), rather than adopting a conventional value.
The final threshold is the calibrated PP3 Supporting threshold of CADD v1.7, PHRED ≥ 25.6
(median of five seeds; section 3.4), which selected 24,926,994 SNVs (the provisional value
had been PHRED ≥ 25); the genome-wide catalogue (section 3.12) adds the SNVs with AVI Phred
≥ 25.18. No prior
filters on consequence, population frequency or ClinVar status were applied, so as not to exclude
potentially deleterious variants that have not been reported to clinical databases.

Because the file contains on the order of 9 × 10⁹ rows and does not fit in memory, it was scanned
in full with tabix queries (Li, 2011) over non-overlapping 2-Mb genomic windows, using
Rsamtools 2.22.0 (htslib via Rhtslib 3.2.0) in R 4.4.2 with 10 parallel processes (package
`parallel`). In each window only rows with PHRED at or above the threshold were retained. Since
all rows above the threshold are kept, the resulting set is exact and independent of the window
size. Selected variants were sorted with data.table 1.16.4 by: (1) PHRED, descending; (2)
RawScore, descending, to break ties between variants with equal (rounded) PHRED; (3) genomic
position, ascending (chromosomes 1–22, X, Y, MT, then position); and (4) reference and
alternative alleles, alphabetically, so that the order is deterministic among alleles at the same
site. Each variant was assigned a consecutive rank.

### 2.1 Validation of the selection procedure

Before applying it to the real data, the procedure was validated on a synthetic set of
13,470,003 SNVs across three chromosomes (1, 2 and X; 4.49 Mb each, three alternative alleles per
position), compressed with bgzip and tabix-indexed in the same format as the CADD file. The set
included forced exact ties (three variants with equal PHRED and RawScore at different positions
on chromosome 1) and positions spanning window boundaries. Scripts: `auxiliary/A1`, `A6`; log:
`logs/A1_A6_synthetic_validation.log`. Both the top 20,000 variants and the full
set with PHRED ≥ 20 (1,423,500 variants) were identical, row by row, to the result of loading and
sorting the complete set in memory with the same criteria. Ties were resolved by genomic
position, as intended.

## 3. Calibration of the CADD v1.7 PHRED threshold

### 3.1 Method

We followed the local posterior probability approach of the ClinGen Sequence Variant
Interpretation (SVI) computational working group (Pejaver et al., 2022), which maps predictor
scores to ACMG/AMP evidence strengths (PP3 for pathogenicity, BP4 for benignity) within the
Bayesian framework of Tavtigian et al. (2018). The original MATLAB implementation (GitHub
`vpejaver/clingen-svi-comp_calibration`, commit `c56d32d`, MIT license) was ported to R
(`executed_pipeline/scripts/04_pejaver_port.R`), preserving its exact semantics. Briefly:

- **Evidence levels.** With prior probability of pathogenicity α = 0.0441 and the constant
  c = 1124 (so that the likelihood ratios for VeryStrong, Strong, Moderate and Supporting are
  c^(1/2^(j−1)) for j = 1…4: 1124, 33.5, 5.79 and 2.41), the posterior probability required for
  each pathogenic level is LR·α / ((LR − 1)·α + 1), and analogously for benign levels with 1 − α.
- **Local posterior.** For each distinct score value t, a window [t − h, t + h] is widened in
  steps of 0.1 PHRED units (the increment the original code uses for CADD PHRED) until it contains
  at least 100 labelled ClinVar variants and at least 3% of the unlabelled gnomAD variants, both
  scaled proportionally at the edges of the score range. The posterior is P / (P + w·B), where P
  and B are the pathogenic and benign counts in the window and w = (1 − α)·n_P / (n_B·α)
  re-weights the benign class to match the prior.
- **Thresholds.** For each level, the threshold is the lowest score in the contiguous region,
  starting from the extreme of the scale, in which the local posterior exceeds the required
  value. Uncertainty was estimated with 10,000 bootstrap replicates (resampling the labelled
  variants and the gnomAD set independently). The final threshold is the one-sided 95% confidence
  bound across replicates. A level is not reached if more than 5% of replicates fail to reach it.
- **Reproducibility.** The original code seeds the random number generator from the clock
  (`rng('shuffle')`). We used a fixed L'Ecuyer-CMRG seed (20260925) so that the parallel bootstrap
  is reproducible. Seeds are arbitrary integers fixed before running: the main-run seed is the
  run date (2026-09-25) written as a number, and the additional seeds of the sensitivity analysis
  (section 3.2.1) are 101, 202, 303 and 404.

The only implementation difference from the original code is computational: window counts are
obtained by binary search on sorted vectors and vectorised across all thresholds. The arithmetic
of the window growth (successive addition of the increment), the inclusive window bounds, and the
threshold selection rules are identical.

### 3.2 Validation of the R port

The port was validated by re-running the published CADD v1.6 PHRED calibration on the authors'
own data (Zenodo 10.5281/zenodo.8347415, `data.zip`, MD5 `a25359e5bbb3de9ea47dd21f03950ffb`):
5,704 P/LP and 6,059 B/LB rare missense variants from ClinVar (2019-12-02) and 363,894 rare
missense variants from gnomAD not present in ClinVar. With B = 10,000 bootstrap replicates
(22.7 min on 10 processes), the R port reproduced the published thresholds (Pejaver et al.,
2022, Table 2): PP3 Supporting 25.3 and Moderate 28.1; BP4 Supporting 22.70 (published 22.7) and
Moderate 17.33 (published 17.3). The local posterior curves (Figure S1) match the CADD panel of
Figure 3 of the original article. The full comparison is given in Table S1.
For BP4 Strong the published threshold is ≤ 0.15, whereas the port gave 0.31 with seed
20260925. This level is estimated in the sparse lower tail of the score distribution (320 B/LB
and 4 P/LP ClinVar variants at PHRED ≤ 0.31): its bootstrap thresholds are widely dispersed
(1st percentile 0.15, median 1.84, 95th percentile 10.2) and 3.4% of replicates do not reach the
level, so the one-sided bound falls at about the 1.7th percentile of valid replicates. The
published value equals the 1st percentile of our replicates. Because the original code seeds the
generator from the clock, both values are single Monte Carlo draws; their seed dependence is
quantified in section 3.2.1. **[VERIFY in Table 2 the published PP3 Strong/VeryStrong and BP4
VeryStrong entries; the port gives all three not reached.]**

#### 3.2.1 Seed (Monte Carlo) sensitivity

The full calibration (B = 10,000) was repeated with four additional seeds (101, 202, 303, 404;
20.6–24.1 min each). Thresholds estimated in data-rich regions were stable: PP3 Supporting
25.2–25.3, PP3 Moderate 28.1 (all seeds), BP4 Supporting 22.6–22.7, BP4 Moderate 17.33–17.34.
BP4 Strong alternated between two values: 0.152–0.155 in three seeds (matching the published
0.15) and 0.310 in two. The switch is explained by the number of replicates that do not reach
the level (3.3–3.7%): because they count toward the 5% bound, small differences move the bound
between two plateaus of the bootstrap distribution separated by a gap between PHRED ~0.16 and
0.31 (Figure S2, Table S2). We conclude that the port reproduces the published calibration, and
that thresholds located in sparse score regions should be reported with their seed variability.
For the CADD v1.7 calibration, thresholds will therefore be reported as the median and range
over at least five seeds. Script: `scripts/05b_seed_sensitivity_cadd16.R`; log:
`logs/05b_seed_sensitivity_cadd16.log`.
Script: `scripts/05_reproduce_pejaver_cadd16.R`; log: `logs/05_reproduce_pejaver_cadd16_B10000.log`;
full results: `calibration/results/cadd16_reproduction_B10000/`.

### 3.3 Coding SNV universe and predicted consequences

All possible SNVs in the coding sequence of the 19,229 GENCODE v50 MANE Select transcripts (one
per protein-coding gene; CDS plus stop codon; chrY PAR copies excluded; 33,463,877 coding
nucleotides at 33,438,724 genomic positions) were enumerated, and their CADD v1.7 scores were
extracted by tabix from the verified local file: 100,316,172 SNVs, with exactly three
alternative alleles at every coding position and no missing or duplicated records
(`scripts/08`, `09`). Reference codons were rebuilt from the reference alleles of the CADD file
itself, on the coding strand, and translated with the standard genetic code (Biostrings);
consequences were assigned per transcript as missense, synonymous, stop gained, stop lost, stop
retained or start lost (`scripts/10`).

As a quality control, the full reference CDS of every transcript was translated and compared with
the GENCODE v50 protein sequence: 19,194 of 19,229 proteins (99.82%) were identical. All 35
differences were explained (`auxiliary/A12`): 35 selenocysteine UGA codons in 25 selenoproteins
(e.g. SELENOP with 10), 10 transcripts with non-AUG initiation codons (ACG, GTG, ATT, ATA, CTG;
annotated as Met by GENCODE), and one transcript with an incomplete terminal codon. The
consequences at these codons were corrected (`scripts/10b`): selenocysteine codons were treated
as Sec (SNVs are missense, or stop gained if they create a stop codon), every SNV in a start codon
was classified as start lost (Sequence Ontology/Ensembl VEP definition), and SNVs in the
incomplete terminal codon as incomplete terminal codon variants. In total 654 SNV × transcript
records were corrected: 280 stop lost → missense and 35 stop retained → stop gained
(selenocysteine), 333 → start lost (37 transcripts with non-ATG start codons: CTG 26, ACG 3,
GTG 2, ATT 2, TTG 2, ATA 2; only 10 had been flagged by the protein comparison because the
translation function reads CTG/TTG at the initiation position as Met), and 6 → incomplete
terminal codon variant (`logs/10b_fix_special_codons.log`, `data/coding/consequence_corrections.tsv`).

The resulting universe comprises 100,316,172 SNVs (100,391,631 SNV × transcript records, owing to
overlapping genes): 72.64% missense, 22.90% synonymous, 4.11% stop gained, 0.17% start lost,
0.15% stop lost and 0.02% stop retained. As an independent check, our consequence agreed with the
ClinVar molecular consequence (MC) field for 99.99% of the 3,405,151 ClinVar SNVs with MC located
in MANE Select coding sequence (missense 99.99%, synonymous 100.00%, stop gained 99.99%, start
lost 99.92%, stop lost 99.83%); most of the ~200 discordances were annotated by ClinVar on a
different transcript (`scripts/11`; `logs/11_consequence_summary_and_clinvar_concordance.log`).

AlphaGenome Variant Impact (AVI) scores (AlphaGenome Atlas, released 2026-09-08; SNV file with
raw and Phred-scaled AVI, GRCh38, tabix) were obtained from the permissive-use download section of
the AlphaGenome website (a single 88,473,344,811-byte ZIP archive with two stored members,
downloaded through a signed-in browser session; SHA-256
`a237c198bc1c033da127fe0257129fdf246e7aa311f8fb584bb053ffa0c7966b`). Both members passed CRC-32
verification before and after extraction, which was done in place because of disk space
(`auxiliary/A11`). The file has the same layout as CADD (`CHROM, POS, REF, ALT, raw_score, PHRED`,
with "chr" sequence names). AVI scores were extracted for the same 100,316,172 coding SNVs
(`scripts/12`): every coding position had exactly three alleles, none was missing or duplicated,
and the reference allele agreed with CADD at every position.

### 3.4 Reference sets and calibration of CADD v1.7 and AVI (main analysis)

The analysis plan was frozen before any CADD v1.7 or AVI calibration result (Git tag
`analysis-plan-v1`, 2026-09-25; `STUDY_PLAN.md`); every later addition or change was recorded as a
dated amendment (A1–A14) committed before the corresponding analysis was run, and results under
the original plan are always reported.

**Reference sets** (following Pejaver et al. 2022 as closely as possible; scripts `13`, `13a`).
ClinVar GRCh38 release 2026-09-23: SNVs that are missense in their GENCODE v50 MANE Select
transcript (our consequence annotation, section 3.3), excluding SNVs with a stop-gained, stop-lost
or start-lost consequence in another overlapping MANE transcript. Labels were taken from exact
CLNSIG values: P/LP = Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic; B/LB = Benign,
Likely_benign, Benign/Likely_benign (uncertain, conflicting and combined terms excluded). Review
status ≥ 1 star. Rare: maximum of the gnomAD v4.1 exome and genome overall allele frequency
< 0.01 (absent = rare). Disease genes: genes (ClinVar GENEINFO) with ≥ 1 P/LP record of any
variant type and ≥ 1 star. The unlabelled gnomAD set, used only to size the local windows, comprises
MANE missense SNVs observed in gnomAD v4.1 (exomes or genomes) with AF < 0.01, absent from ClinVar
(any record) and located in the same disease genes. The final sets contained 59,979 P/LP and
77,128 B/LB variants (4,826 genes) and 1,955,995 gnomAD variants. Declared deviations from the
original study: no exclusion of HGMD variants (licensed resource; see the HGMD sensitivity analysis,
section 3.8), no exclusion of predictor-training variants (neither CADD nor AVI is trained on
ClinVar), and no gnomAD genotype-quality or low-complexity filters (not available in the VEP cache).

**Scores.** CADD v1.7 PHRED as distributed (one decimal); AVI Phred rounded to two decimals.

**Calibration.** Validated R port (section 3.2; the faster implementation `04b_pejaver_port_fast.R`
counts window members with `tabulate` instead of sorting and was verified to give bit-identical
results, `auxiliary/A14`); α = 0.0441, c = 1124, B = 10,000, minimum 100 ClinVar variants and 3%
of the gnomAD set per window, window increment 0.1 for both Phred scales, one-sided 95% bootstrap
bound. Main analysis: five seeds (20260925, 101, 202, 303, 404); thresholds reported as the median
and range across seeds (script `14_calibrate_cadd17_avi.R`; summary `15`; Figure 1, script `16`).

### 3.5 Sensitivity analyses and discrimination

Four sensitivity analyses were pre-specified (one seed each): S1, review status ≥ 2 stars; S2,
rare defined as AF < 0.001 (AVI is trained with gnomAD frequency proxy labels, so frequency-based
B/LB classifications could favour it); S3, gene-level (cluster) bootstrap, resampling genes instead
of variants in both sets; S4, exclusion of missense SNVs within 3 nt of an internal exon boundary.
Discrimination: AUROC and AUPRC on the main reference set with 95% bootstrap intervals (1,000
variant-level replicates); the paired AVI − CADD difference (1,000 paired replicates) was added by
amendment A1 (script `17`).

### 3.6 Temporal hold-out validation (Aim 2c; amendment A2)

Both predictors were recalibrated (seed 20260925, B = 10,000) on reference sets built with the same
criteria from the ClinVar release of 2024-01-27, the last release before the CADD v1.7 score files.
The test set comprised variants of the 2026 main set absent from the 2024 set (variants whose label
changed between releases were excluded and counted). For each evidence interval, the positive
likelihood ratio on the test set (95% bootstrap interval, 1,000 replicates) was compared with the
likelihood ratio required for that evidence strength. Because AVI was released in 2026, this tests
the stability of its calibration on newly classified variants rather than a post-release hold-out
(scripts `18`, `19`; Figure S3).

### 3.7 Comparison with already-recommended predictors (Aim 6; amendment A3)

AlphaMissense (hg38 predictions, Zenodo 10813168), REVEL v1.3 (official file, GRCh38 positions)
and BayesDel without allele frequency (UCSC hg19 bigWig tracks, one per alternative allele, lifted
over to GRCh38 with the UCSC chain file and checked against the hg19 reference base; minus-strand
mappings complemented) were extracted for all reference variants (script `20`); where a predictor
had several values for an SNV, the maximum was used. Precision and window increments matched the
original study (BayesDel two decimals, increment 0.01; REVEL three decimals and AlphaMissense four
decimals, increment 0.001). The main comparison used the reference variants scored by all five
predictors (intersection), five seeds per predictor, with AUROC/AUPRC and paired bootstrap
differences versus AVI (script `21`).

### 3.8 Threshold drift, circularity and HGMD exclusion (amendments A6, A7)

*Release series (A6a).* Reference sets were rebuilt with the same criteria from ClinVar releases of
2019-12-02 (used by Pejaver et al.), 2022-11-19 (last release before the ClinGen recommendations),
2024-01-27 and 2026-09-23, and all five predictors were calibrated on the variants of each release
scored by all five (one seed; scripts `24a`–`24c`). *Before/after split (A6b).* Within the 2026 main
set, predictors were calibrated separately on variants already present in the November 2022 set and
on variants added afterwards. Because AVI (released 2026-09-08) cannot have been used as evidence in
any ClinVar classification, it served as a negative control: a threshold relaxation on post-2022
variants present for REVEL/BayesDel/CADD but absent for AVI was pre-specified as supporting
circularity. A post hoc common-scale comparison is reported separately (script `24d`).
*HGMD exclusion (A7).* HGMD Professional is licensed; the Ensembl VEP 116 cache contains HGMD Public
2020.4 records (position only). Every SNV at the position of a single-nucleotide HGMD Public record
was excluded from the ClinVar and gnomAD sets of the 2019 and 2026 releases, and the five predictors
were recalibrated (scripts `25a`, `25b`). The HGMD-included analysis remains primary; the excluded
analysis is a sensitivity analysis.

### 3.9 Calibration by gene class (Aim 1; amendment A8)

Using the 2026 intersection sets, genes (MANE Select gene; SNVs assigned to more than one gene
excluded) were classified by (1) constraint: gnomAD v4.1 LOEUF decile of the MANE transcript —
highly constrained (deciles 0–2), intermediate (3–6), tolerant (7–9); (2) inheritance: ClinGen
gene–disease validity (Definitive, Strong or Moderate), or GenCC with the same classifications when
ClinGen had none — autosomal dominant only vs autosomal recessive only; (3) mechanism:
Gene2Phenotype molecular mechanism (confidence definitive, strong or moderate) — loss of function vs
gain of function or dominant negative. For each class the gnomAD set was restricted to the class's
genes; classes with < 1,000 P/LP or < 1,000 B/LB variants were not calibrated. All five predictors
were calibrated per class (one seed, global prior). Two classes were reported as different for a
threshold when the 2.5th–97.5th percentile intervals of their bootstrap thresholds did not overlap;
the percentage of variants changing evidence level under class-specific thresholds is also reported
(scripts `26`, `27`; Figure 2).

### 3.10 Evidence yield for variants of uncertain significance (Aim 7; amendment A5)

ClinVar 2026-09-23 SNVs with CLNSIG exactly "Uncertain_significance", ≥ 1 star, MANE Select
missense (same annotation and exclusions) in disease genes were assigned the PP3/BP4 level of each
predictor (main-analysis thresholds for CADD and AVI; intersection thresholds for all five), with
the corresponding ACMG/AMP points (Tavtigian et al. 2020). No frequency filter was applied in the
main count (a secondary count is restricted to AF < 0.01). This quantifies computational evidence
only (script `23`).

### 3.10b VUS resolved between ClinVar releases (Aim 7b; amendment A15)

SNVs classified Uncertain significance (>= 1 star) in the ClinVar release of 27 January 2024 that are
MANE Select missense in a disease gene (disease genes from the 2024 release; same exclusions as the
reference sets) were followed to the release of 23 September 2026 and grouped as resolved to P/LP,
resolved to B/LB or not resolved. Primary analysis: gnomAD v4.1 AF_max < 0.01 and thresholds
calibrated on the January 2024 release only (amendment A6 release series, set r2024, five-predictor
intersection, one seed). Secondary: 2026 main and intersection thresholds, and all VUS. For each
predictor: percentage of resolved variants with PP3 or BP4 evidence and mean ACMG/AMP points by
outcome; concordance with the final classification; likelihood ratio of each interval among resolved
variants (P/LP vs B/LB, 2,000 bootstrap replicates) against the ClinGen targets. AVI (not citable by
ClinVar submitters before 2026) is the negative control for circularity. Illustrative examples follow
a rule fixed before inspection (2026 record reviewed by an expert panel or practice guideline; up to
three per direction with the largest absolute AVI points). Script `42_vus_resolved.R`; Table 2,
Figure 3B, Tables S18–S19.

### 3.11 Population frequency and ancestry (Aim 8; amendments A4, A12)

*gnomAD (A4).* Universe: MANE Select missense SNVs observed in gnomAD v4.1. Per-group AF = maximum
of exome and genome AF for AFR, AMR, ASJ, EAS, FIN, MID, NFE and SAS, and genome AF for AMI. For
PP3-level SNVs (main median thresholds) we report the number with a group AF > 0.01 (proxy for a
gene-agnostic BS1 level) and > 0.05 (BA1), the group with the highest AF, per-group proportions
with Wilson intervals, and group-specific common variants (script `22`). The VEP cache has no
allele numbers, so the ClinGen minimum-AN requirement for BA1 could not be applied.
*Mexico City Prospective Study (A12a).* The official per-chromosome allele-frequency VCFs of the
RGC-MCPS browser (release 20230130; ~140,000 adults; frequencies overall and by local ancestry:
African, European and Indigenous Mexican) were downloaded (script `31`; sizes verified) and queried
by tabix at the positions of the PP3-level MANE missense SNVs of the Aim 4 catalogue, matching
position and alleles exactly (script `32`). For each predictor and PP3 level we report the number
observed, with MCPS AF > 0.01 and > 0.05, and with MCPS AF > 0.01 while ≤ 0.01 in every gnomAD
group ("MCPS-only frequency conflicts"), overall and for the Indigenous Mexican component; ClinVar
P/LP variants were counted as a quality check, and the MCPS AF was correlated with the gnomAD AMR AF
(Spearman). *TOPMed BRAVO (A12b).* The public Freeze 10 VCFs (overall AC/AF only; downloaded by the
author from signed links; sizes verified) were recompressed with bgzip, indexed, and analysed in
the same way using FILTER = PASS sites (script `34`). gnomAD v4.1 genomes include many TOPMed
samples, so BRAVO is not independent of gnomAD; MCPS participants are not expected to overlap
with either. The resources are reported separately and never pooled; no threshold is recalibrated
with these data.

### 3.12 Genome-wide catalogue (Aim 4; amendment A9)

Every SNV of the CADD v1.7 and AVI files with CADD PHRED ≥ 25.6 or AVI Phred ≥ 25.18 (main median
PP3 Supporting thresholds) was retrieved by full tabix scans (scripts `02`, `28`) and annotated with
both scores, the PP3 level of each predictor, the MANE Select consequence (coding SNVs) or overlap
with GENCODE v50 gene bodies (other SNVs; "intergenic" otherwise), ClinVar 2026-09-23 (CLNSIG,
stars, variation ID; exact allele match) and gnomAD v4.1 exome, genome and per-group AF and grpmax
(VEP 116 cache) (script `29`). A column `calibrated_context` is TRUE only for MANE Select missense
SNVs; for all other SNVs the PP3 level is provided for prioritisation only. A per-gene summary
(MANE missense) reports catalogued SNVs, SNVs at each PP3 level, gnomAD observation, common variants
(group AF > 1% and > 5%; overall AF > 1% for comparison; script `29b`) and ClinVar P/LP and B/LB
counts.

### 3.13 Noncoding variants (Aim 9; amendments A10, A11)

Every ClinVar 2026-09-23 SNV outside MANE Select coding sequence was assigned, relative to the MANE
Select transcripts, to one region class (canonical splice ±1–2, intronic splice region ±3–8,
near-intronic 9–50, deep intronic > 50 nt, 5′UTR, 3′UTR, upstream and downstream ≤ 5 kb,
intergenic), and coding synonymous SNVs were added as a class split by distance to an internal exon
boundary (≤ 3 nt vs other). Reference variants follow the plan-v1 rules. Classes with ≥ 1,000 P/LP
and ≥ 1,000 B/LB variants would be calibrated; otherwise the missense-calibrated thresholds are
evaluated with interval-based likelihood ratios and AUROC. In addition, CADD and AVI were evaluated
against saturation-mutagenesis MPRA data of 20 regulatory elements (Kircher et al. 2019; functional
= p < 10⁻⁵ and |coefficient| ≥ 0.1) (script `30`). SpliceVarDB-based splicing calibration was
deferred to a second publication (A11). **[PENDING: results.]**

### 3.14 Validation with multiplexed assays of variant effect (Aim 3; amendment A13)

MaveDB (accessed 2026-09-27) was searched for published human score sets with a score calibration
defining both "normal" and "abnormal" functional classes (metadata inventory, script `33a`; 92
score sets in ~40 genes). For each set the investigator-provided calibration was used when
available (otherwise the first calibration not flagged research-use-only). Range-based classes
were applied to the scores; class-based (categorical) classes were retrieved as variant lists from
the MaveDB calibration endpoint. Variants were mapped to GRCh38 with MaveDB's own mappings:
nucleotide-level sets through genomic HGVS (SNVs only; current mapping, else the most recent), and
protein-level sets through the mapped RefSeq protein change, which was matched to all MANE Select
missense SNVs encoding the same amino-acid change; protein-level sets were kept only if ≥ 90% of
their mapped (position, reference residue) pairs agreed with the MANE Select protein (script
`33b`). SNVs classified discordantly by different sets were excluded. Evaluation (script `33c`):
for each evidence interval of the calibrated thresholds, the proportion of abnormal variants and
LR = P(interval | abnormal) / P(interval | normal) with 95% bootstrap intervals (B = 2,000);
AUROC pooled and per gene (≥ 30 abnormal and ≥ 30 normal); paired AVI − CADD difference;
secondary analysis on the SNVs scored by all five predictors; sensitivity analysis excluding
variants of the ClinVar calibration set. Functional abnormality is not equivalent to
pathogenicity, so the LRs test direction and consistency rather than the ClinGen LR targets.
**[PENDING: results.]**

### 3.15 CADD v1.6 versus v1.7 (Aim 2b; amendment A14)

CADD v1.6 GRCh38 scores were obtained by partial download of the official whole-genome file (86.6
GB): the tabix linear index was parsed to find the BGZF blocks covering the positions of the 2026
reference sets (merged into 6,419 byte ranges, ~9.4 GB; script `35a`), the ranges were downloaded
(script `35b`) and assembled into a valid BGZF file of complete blocks, which was indexed and queried
with exact allele matching (script `35c`). CADD v1.6 is then calibrated on the same reference
variants as the main analysis (five seeds), and thresholds, discrimination and evidence levels are
compared with v1.7. **[PENDING: results.]**

## 4. Variant annotation

### 4.1 Gene
Coding SNVs are assigned to the gene of their GENCODE v50 MANE Select transcript (section 3.3).
Other SNVs of the catalogue are annotated with every GENCODE v50 *basic* gene body they overlap on
GRCh38 (released 2026-06-09; Frankish et al., 2023 **[VERIFY most recent citation]**), with gene
identifiers and biotypes, using GenomicRanges 1.58.0; SNVs overlapping no gene are labelled
intergenic.

### 4.2 ClinVar
Clinical classifications are taken from the ClinVar GRCh38 VCF (`fileDate=2026-09-23`; archived
releases 2019-12-02, 2022-11-19 and 2024-01-27 for Aims 2 and 6; Landrum et al. **[VERIFY
citation]**) by exact match of chromosome, position, reference and alternative allele (CLNSIG,
review stars derived from CLNREVSTAT, variation ID). ClinVar labels define the calibration
reference sets; in the genome-wide catalogue ClinVar is only a descriptive annotation and is
**not a selection criterion**.

### 4.3 Population frequencies (gnomAD)
Frequencies are obtained in two stages:

1. **All selected variants.** gnomAD allele frequencies for exomes and genomes, overall and per
   population, from the Ensembl VEP release 116 indexed cache for GRCh38
   (`homo_sapiens_vep_116_GRCh38.tar.gz`, released 2026-04-09), read directly with tabix, with
   exact matching on chromosome, position and allele. The cache's `info.txt` confirms gnomAD
   v4.1 for both exomes and genomes (`source_gnomADe v4.1`, `source_gnomADg v4.1`) and the same
   gene build used for gene annotation (GENCODE 50). Only the dbSNP (`rs`) records carry
   frequencies, stored per allele (`ALT:AF`, four significant digits). As a check, the cache
   values for TP53 c.524G>A (chr17:7675088 C>T; exomes 4.104 × 10⁻⁶) and C>G (6.841 × 10⁻⁷)
   agreed with the gnomAD v4.1 exomes VCF (4.10437 × 10⁻⁶ and 6.84052 × 10⁻⁷)
   (`auxiliary/A10`, `A2`). Integrity of the 27,644,657,162-byte archive was checked by exact
   size and by full decompression (gzip CRC).
2. **Variants prioritised at the end of the analysis.** Allele count (AC), allele number (AN),
   number of homozygotes (nhomalt), combined exome + genome frequency (AF_joint), population with
   the highest frequency (grpmax) and filter status (FILTER), from the gnomAD v4.1 *joint* VCF
   (`gs://gcp-public-data--gnomad/release/4.1/vcf/joint/`) **[VERIFY citation]** through remote
   tabix queries. **[PENDING: not yet run; planned for the variants with frequency conflicts.]**

Throughout, "overall AF" denotes the maximum of the gnomAD exome and genome overall AF (used for
the plan-v1 rarity rule), and "group AF" the maximum over genetic ancestry groups of the per-group
exome and genome AF (used for common-variant counts). An early version of the per-gene catalogue
summary and a first run of the MCPS analysis used the overall AF where the group AF was specified;
both were corrected (scripts `29b`, `32`) and the corrections are logged.

The two-stage design was necessary because querying the *joint* VCF for all variants was not
feasible. In a test over the coding regions of chromosome 21 (GENCODE CDS ± 10 bp; 300 regions,
0.052 Mb), the remote query returned 41,328 records and 494 MB of data (~12 KB per record, 479
INFO fields) in 539 s. Extrapolated to the coding regions of the whole genome alone (41.4 Mb),
the volume would exceed 390 GB. Variants absent from gnomAD will be reported as "not observed",
not as zero frequency.

## 5. Software

| Component | Version |
|---|---|
| Operating system | Windows 11 Enterprise 10.0.26200 |
| R | 4.4.2 (2024-10-31 ucrt) |
| Bioconductor | 3.20 |
| data.table | 1.16.4 |
| Rsamtools | 2.22.0 |
| Rhtslib (htslib) | 3.2.0 |
| GenomicRanges | 1.58.0 |
| R.utils | 2.13.0 |
| ggplot2 | 3.5.2 |
| patchwork | 1.3.2 |
| ragg | 1.4.0 |
| curl | 8.21.0 (Schannel) |
| PowerShell | 5.1.26100 |

An `renv` lockfile of the packages loaded by the executed scripts and their dependencies (86
packages; R 4.4.2, Bioconductor 3.20) is provided in `executed_pipeline/renv.lock` and can be
restored with `renv::restore(lockfile = "executed_pipeline/renv.lock")`.
Full `sessionInfo()` output: `executed_pipeline/R_sessionInfo.txt`; every installed R package
with version and repository: `executed_pipeline/package_inventory.tsv`. Hardware (Intel Core Ultra
7 155H, 16 cores / 22 threads, 31.5 GB RAM) and run times: `executed_pipeline/compute_environment.txt`.

## 6. Code and data availability

Scripts exactly as executed, with SHA-256 hashes, are in `executed_pipeline/scripts/` (main
pipeline) and `executed_pipeline/scripts/auxiliary/` (validation and benchmarking); console
output of each run is in `executed_pipeline/logs/`. The file `execution_log.ps1` lists the
commands in the order they were run. Every external data resource is listed in
`executed_pipeline/data_manifest.tsv` with URL, version, release and download dates, size,
checksum, licence and citation **[VERIFY licence terms, in particular CADD's terms for
non-commercial use, before deposit]**.

## Figures

Each figure is stored as PNG and PDF in `figures/`, together with the data needed to redraw it
(`figures/data/<figure>/`) and the command that generates it, so it can be modified without
repeating the computations.

| Figure | Content | Files | Status |
|---|---|---|---|
| S1 | Local posterior probability of CADD v1.6 PHRED, reproduced from the data of Pejaver et al. (2022); validation of the R port | `figures/FigS1_cadd16_reproduction.{png,pdf}`; data `figures/data/FigS1_cadd16_reproduction/` | Done 2026-09-25 |
| Table S1 | Published vs reproduced CADD v1.6 thresholds | `figures/TableS1_port_validation.tsv` | Done 2026-09-25 (PP3 Strong/VeryStrong and BP4 VeryStrong rows to verify) |
| S2 | Seed sensitivity of the one-sided BP4 Strong bound (5 seeds × 10,000 replicates) | `figures/FigS2_seed_sensitivity.{png,pdf}`; data `figures/data/FigS2_seed_sensitivity/`; script `scripts/06b_plot_seed_sensitivity.R` | Done 2026-09-25 |
| Table S2 | All thresholds per seed, with fraction of failed replicates | `figures/TableS2_seed_sensitivity.tsv` | Done 2026-09-25 |
| S3 | Temporal validation: interval likelihood ratios on variants classified after ClinVar 2024-01-27 | `figures/FigS3_temporal_validation.{png,pdf}`; data `figures/data/FigS3_temporal_validation/`; script `scripts/19_plot_temporal_lr.R` | Done 2026-09-25 |
| 5 (Figure 2 until manuscript v4) | Calibrated PP3/BP4 thresholds by gene class (constraint, inheritance, mechanism), five predictors | `figures/Fig5_gene_classes.{png,pdf}`; data `figures/data/Fig5_gene_classes/`; script `scripts/27_plot_gene_classes.R` | Done 2026-09-26; renumbered 2026-10-04 |
| 1 | Local posterior probability of CADD v1.7 PHRED and AVI Phred (main analysis, seed 20260925) | `figures/Fig1_cadd17_avi_calibration.{png,pdf}`; data `figures/data/Fig1_cadd17_avi_calibration/`; script `scripts/16_plot_figure1.R` | Done 2026-09-25 |
| 2 (Figure 3 until manuscript v4) | Head-to-head comparison on the intersection set: AUROC and evidence levels of P/LP and B/LB variants, five predictors | `figures/Fig2_comparators.{png,pdf}`; data `figures/data/Fig2_comparators/`; script `scripts/36_plot_fig3_comparators.R` | Done 2026-09-27; renumbered 2026-10-04 |
| Not in manuscript v6 (Figure 4 until v4, Figure S5 in v5; population analyses reserved for a separate study) | Population frequency of PP3-level missense SNVs: gnomAD group with highest AF; MCPS overall and Indigenous Mexican frequencies | `figures/FigS5_population.{png,pdf}`; data `figures/data/FigS5_population/`; script `scripts/37_plot_fig4_population.R` | Done 2026-09-27; moved to the supplement 2026-10-04 |
| 3 (Figure 5 until manuscript v4, now panel A) | Computational evidence for ClinVar missense VUS, five predictors: (A) 2026 VUS, intersection thresholds; (B) 2024 VUS resolved by 2026, thresholds calibrated on the 2024 release (amendment A15) | `figures/Fig3_vus.{png,pdf}`; data `figures/data/Fig3_vus/`; script `scripts/44_plot_fig3_vus.R` (replaces `38_plot_fig5_vus.R`) | Done 2026-10-04 |
| 4 | Missense-calibrated thresholds outside missense variation: likelihood ratios of the BP4 intervals and AUROC of CADD v1.7 and AVI per ClinVar region class (grayscale) | `figures/Fig4_noncoding.{png,pdf}`; data `figures/data/Fig4_noncoding/`; script `scripts/43_plot_fig_noncoding.R` | Done 2026-10-04 |
| Table 2 (manuscript v5) | Calibrated evidence for VUS of January 2024 resolved by September 2026 (amendment A15); the former Table 2 (MAVE) is now Table S17 | `figures/Table2_vus_resolved.{tsv,md}`; data `figures/data/Table2_vus_resolved/`; script `scripts/45_make_table2_vus_resolved.R` | Done 2026-10-04 |
| Tables 1–2 | Table 1: calibrated thresholds of CADD v1.7 and AVI with evidence levels and published CADD v1.6 thresholds; Table 2: likelihood ratios of the calibrated intervals on MAVE functional classes, five predictors | `figures/Table1_thresholds.tsv`, `figures/Table2_mave.tsv`, `figures/Tables_1_2.docx`; data `figures/data/Table1_thresholds/`, `figures/data/Table2_mave/`; script `scripts/39_make_tables_1_2.R` | Done 2026-09-27 (Table 1 gains a CADD v1.6 2026-data column when A14 finishes) |
| S4 | AUROC of five predictors for MAVE functional classes per gene (36 genes) | `figures/FigS4_mave_genes.{png,pdf}`; data `figures/data/FigS4_mave_genes/`; script `scripts/40_plot_figS4_mave_genes.R` | Done 2026-09-27 |
| Tables S1–S18 (manuscript v6) | Supplementary tables (reproduction, seeds, predictor settings, sensitivity, HGMD, temporal hold-out, comparators, release series, gene classes, VUS, splice region, noncoding and MPRA, MaveDB score sets (S13), CADD v1.6 vs v1.7 (S14), plan and amendments (S15), MAVE intervals (S16), resolved VUS (S17), illustrative resolved VUS (S18); the MCPS/BRAVO table was removed in v6) | `figures/Supplementary_Tables_S1-S18.xlsx` (Contents sheet with titles and legends); TSV copies in `figures/data/SupplementaryTables/`; script `scripts/41_make_supplementary_tables.R` | Done 2026-09-27; S16–S18 added 2026-10-04; renumbered 2026-10-05 |

## References [VERIFY all before submission]

- Schubach M, Maass T, Nazaretyan L, Röner S, Kircher M. CADD v1.7: using protein language models, regulatory CNNs and other nucleotide-level scores to improve genome-wide variant predictions. *Nucleic Acids Res.* 2024;52(D1):D1143–D1154.
- Kircher M, et al. A general framework for estimating the relative pathogenicity of human genetic variants. *Nat Genet.* 2014;46:310–315.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet.* 2022;109(12):2163–2177.
- Tavtigian SV, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med.* 2018;20(9):1054–1060.
- Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. *Bioinformatics.* 2011;27(5):718–719.
- Frankish A, et al. GENCODE: reference annotation for the human and mouse genomes in 2023. *Nucleic Acids Res.* 2023;51(D1):D942–D949.
- ClinVar and gnomAD v4.1: citations to be confirmed.

---

## Appendix: process log

| Date | Event |
|---|---|
| 2026-09-24 | WSL not available on the machine (not installed; firmware virtualisation disabled; no administrator rights). A native Windows implementation with R + Rsamtools was adopted. |
| 2026-09-24 | Rsamtools and dependencies installed in the user R library. |
| 2026-09-24 | Validation on synthetic data (section 2.1): a defect was found and fixed (parallel workers did not load data.table); the corrected version passed all comparisons. |
| 2026-09-24 | Single-connection download failed a certificate revocation check on the local network; `--ssl-no-revoke` was added to curl (the certificate is still validated). At 0.25–0.68 MB/s it was not feasible, so 8 parallel byte ranges were used instead (~1.4–2 MB/s). |
| 2026-09-24 | Selection criterion changed from "top 500,000 PHRED" (initial, arbitrary) to a PHRED threshold (provisionally 25). Decided not to filter by ClinVar. |
| 2026-09-24 | GENCODE v50 basic and ClinVar (fileDate 2026-09-23) downloaded. The Ensembl VEP REST API was unavailable (HTTP 500/503). MyVariant.info was discarded because of outdated versions (gnomAD 2.1.1, ClinVar 2025-05). Remote tabix queries to gnomAD v4.1 were verified to work from Windows. |
| 2026-09-24 | Performance test of gnomAD *joint* (section 4.3) showed that querying it for all variants is not feasible; a two-stage strategy was adopted and the VEP cache (26 GB) download started. |
| 2026-09-25 | Decided to calibrate the threshold with the method of Pejaver et al. (2022) instead of using PHRED 25 or the published CADD v1.6 value (25.3), which predates CADD v1.7 and current ClinVar/gnomAD. Code and data of the original study obtained (GitHub, Zenodo). |
| 2026-09-25 | Count of current ClinVar SNVs (≥ 1 star): missense 62,231 P/LP and 147,638 B/LB, enough for calibration. Stop-gained and canonical splice variants have almost no B/LB, and synonymous/UTR variants almost no P/LP, so a genome-wide threshold cannot be calibrated for every consequence class; the missense-calibrated threshold will be applied genome-wide and this will be stated as a limitation. |
| 2026-09-25 | MATLAB code ported to R. A 20-replicate test on the CADD v1.6 data gave thresholds close to the published ones. |
| 2026-09-25 | All deliverables (this document, execution log, figures, scripts) switched to English. The download scripts already executed (01, 03) were kept verbatim with Spanish comments. |
| 2026-09-25 | R port validated: B = 10,000 reproduction of the CADD v1.6 calibration matched the published thresholds (Table S1, Figure S1). |
| 2026-09-25 | Reproducibility records added: auxiliary scripts A1–A7 archived verbatim (A1–A5 were run with Spanish comments/labels), console logs, data manifest, package inventory, compute environment, README. Synthetic validation and ClinVar counts re-run from the archived scripts to produce logs; results unchanged. Correction: the forced ties of the synthetic test are all on chromosome 1 (an earlier draft said "different chromosomes"). |
| 2026-09-25 | Study scope agreed with the author: Aims 1–4 (stratified calibration by gene class, temporal/version stability, MAVE validation, genome-wide resource); Aim 5 (local population cohort) recorded as a future proposal because no cohort data are currently available. Literature scan performed to position each aim. See `STUDY_PLAN.md`. |
| 2026-09-25 | Aim 1 gene annotations downloaded (gnomAD v4.1 constraint; GenCC; ClinGen gene validity; ClinGen dosage; Gene2Phenotype), with provenance and hashes. Caveat for the mechanism stratum: 3,642 of 3,851 Gene2Phenotype records have an *inferred* (not evidence-based) molecular mechanism, and only ~330 records are gain of function or dominant negative. |
| 2026-09-25 | Repository prepared for GitHub (private until publication): local user-profile paths redacted in versioned text files (A9; only path strings changed, A1 re-hashed); the single-connection download script that produced the first 84 MB of the CADD file archived as script 00; Git initialised with byte-for-byte storage so that published hashes remain verifiable; third-party data excluded and documented in the data manifest. |
| 2026-09-25 | Scope extended to increase clinical impact (STUDY_PLAN.md): head-to-head calibration of AVI (AlphaGenome Atlas, released 2026-09-08), CADD v1.7 and already-recommended predictors (Aim 6); VUS evidence yield (Aim 7); temporal hold-out validation (added to Aim 2); ancestry analysis (Aim 8). Noncoding/splicing calibration deferred (Aim 9). Two outputs planned (main article, resource article); target journals recorded. |
| 2026-09-25 | The reproduced BP4 Strong threshold (0.31) was questioned. Diagnosis: point estimate 1.90; bootstrap BP4 Strong thresholds are widely dispersed (1st percentile 0.15, median 1.84, 95th percentile 10.2; 3.4% of replicates fail the level, with ~450 ClinVar variants below PHRED 1). Because failed replicates count toward the 5% bound, the reported value is the ~1.7th percentile of valid replicates, a tail region sensitive to the random seed. A seed-sensitivity analysis (4 additional seeds, B = 10,000) was started (step 5c). Author confirmed non-commercial use of CADD. |
| 2026-09-25 | Seed sensitivity completed: BP4 Strong = 0.152–0.155 in 3/5 seeds (published 0.15) and 0.31 in 2/5; all other levels stable (section 3.2.1). A first draft of Figure S2 excluded the failed replicates from the cumulative curve, making the 5% line cross at ~0.5; it was corrected so that failed replicates are counted first, as in the original rule. Decision: CADD v1.7 thresholds will be reported as median and range over ≥ 5 seeds. |
| 2026-09-25 | Both parallel downloads (VEP cache, CADD) completed all parts but stopped before merging because a curl retry warning was treated as a terminating error; fixed in script v3 (03b), which also checks every part's exact length before merging. VEP cache verified (size, SHA-256, gzip CRC); gnomAD v4.1 confirmed in its info.txt. |
| 2026-09-25/26 | Main calibration, sensitivity analyses, temporal hold-out, comparators, ancestry, VUS, drift/circularity (A6), HGMD exclusion (A7), gene classes (A8) and the genome-wide catalogue (A9) run as registered (RESULTS R1–R12). Two jobs were killed by Claude Code under low memory (scripts 13 and 18); both were rewritten to be memory-lean and rerun with the author's approval. The AVI scan was stopped preventively at 8 workers and resumed with 6. The original bootstrap implementation would have taken ~30 h for the comparator runs; a faster port (window counts by `tabulate`) was verified bit-identical (A14 auxiliary) and used thereafter. |
| 2026-09-26 | Decision with the author: the HGMD-included analysis stays primary; HGMD-excluded is a sensitivity analysis (R10 shows that HGMD inclusion explains the discrepancy with the published REVEL/BayesDel thresholds). Noncoding aim activated (A10); SpliceVarDB deferred to a second paper (A11) because splicevardb.org is not reachable from the local network. |
| 2026-09-27 | Population resources assessed (A12): MCPS per-chromosome frequency VCFs found on the RGC-MCPS Resources page (pointed out by the author) and downloaded; RGC Million Exome and All of Us offer no bulk download usable here (the RGC browsers' internal API is not scraped); TOPMed BRAVO VCFs downloaded from signed links obtained by the author. Sample overlap between BRAVO and gnomAD genomes noted (not independent). |
| 2026-09-27 | Two definition errors found and corrected before interpretation: the per-gene catalogue summary and a first MCPS run used the gnomAD overall AF where the group maximum was specified (scripts 29b, 32 rerun). An MCPS conflict (ALDH7A1 p.Asn195Ser, AF 0.28, array-only call) was identified as a probable array artefact; variant source is reported post hoc. |
| 2026-09-27 | MAVE validation registered (A13) from MaveDB metadata before any score was downloaded. A first run of the label builder had two implementation errors (class-based calibrations not applied; the protein-level QC compared whole amino-acid changes, of which only ~1/3 are SNV-reachable, instead of reference residues); both corrected before any evaluation. CADD v1.6 comparison registered (A14); v1.6 scores obtained by partial BGZF download planned from the tabix index (~9.4 GB of 86.6 GB). |
| 2026-09-27 | The Claude Code session ended while three jobs were running (noncoding lookups after >4 h, MAVE labels, CADD v1.6 download); all were restarted, the download resuming from its completed parts, and the noncoding lookups rewritten as cached per-chromosome parallel jobs. |
