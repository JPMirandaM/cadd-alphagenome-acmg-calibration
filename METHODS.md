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
integrity was checked by requiring the final size to match the size reported by the server
exactly. **[PENDING: completion date, SHA-256 of the final file, and total number of rows
scanned in section 2.]**

## 2. Variant selection

The CADD PHRED score is defined as −10·log₁₀(r/N), where r is the rank of the variant among the
N possible SNVs of the genome. Thus PHRED ≥ 10, ≥ 20 and ≥ 30 correspond to the top 10%, 1% and
0.1% of variants. The selection threshold was derived from our own calibration of CADD v1.7
against current ClinVar and gnomAD data (section 3), rather than adopting a conventional value.
**[PENDING: final threshold. Provisional value: PHRED ≥ 25, i.e. the top ~0.32%.]** No prior
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
site. Each variant was assigned a consecutive rank. **[PENDING: number of selected variants and
their distribution across PHRED bins.]**

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

### 3.3 Application to CADD v1.7 with current data [PENDING]

The reference sets will be rebuilt with current data, following the inclusion criteria of the
original study as closely as possible. **[PENDING: detailed definition of the ClinVar
(2026-09-23 release) P/LP and B/LB missense sets and of the gnomAD set; documented deviations
from the original criteria (e.g., HGMD exclusion is not reproducible without a licence); number
of variants; resulting thresholds; Figure 1.]**

## 4. Variant annotation [PENDING]

### 4.1 Gene
Variants will be assigned to genes by coordinate overlap with the GENCODE v50 *basic* annotation
on GRCh38 (released 2026-06-09; Frankish et al., 2023 **[VERIFY most recent citation]**), using
GenomicRanges 1.58.0. **[PENDING: exact criterion, e.g., the gene containing the variant, its
biotype (protein-coding, lncRNA, etc.), and handling of variants in more than one gene or in
intergenic regions.]**

### 4.2 ClinVar
Clinical classifications will be taken from the ClinVar GRCh38 VCF (`fileDate=2026-09-23`;
Landrum et al. **[VERIFY citation]**) by exact match of chromosome, position, reference and
alternative allele. We will report clinical significance (CLNSIG), review status (CLNREVSTAT),
associated condition (CLNDN) and the ClinVar variation ID. ClinVar is used only as a descriptive
annotation, **not as a selection criterion**.

### 4.3 Population frequencies (gnomAD)
Frequencies will be obtained in two stages:

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
   tabix queries.

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
| 1 | Local posterior probability of CADD v1.7 PHRED with current ClinVar/gnomAD data | `figures/Fig1_cadd17_calibration.*` | [PENDING] |

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
| 2026-09-25 | The reproduced BP4 Strong threshold (0.31) was questioned. Diagnosis: point estimate 1.90; bootstrap BP4 Strong thresholds are widely dispersed (1st percentile 0.15, median 1.84, 95th percentile 10.2; 3.4% of replicates fail the level, with ~450 ClinVar variants below PHRED 1). Because failed replicates count toward the 5% bound, the reported value is the ~1.7th percentile of valid replicates, a tail region sensitive to the random seed. A seed-sensitivity analysis (4 additional seeds, B = 10,000) was started (step 5c). Author confirmed non-commercial use of CADD. |
| 2026-09-25 | Seed sensitivity completed: BP4 Strong = 0.152–0.155 in 3/5 seeds (published 0.15) and 0.31 in 2/5; all other levels stable (section 3.2.1). A first draft of Figure S2 excluded the failed replicates from the cumulative curve, making the 5% line cross at ~0.5; it was corrected so that failed replicates are counted first, as in the original rule. Decision: CADD v1.7 thresholds will be reported as median and range over ≥ 5 seeds. |
| 2026-09-25 | Both parallel downloads (VEP cache, CADD) completed all parts but stopped before merging because a curl retry warning was treated as a terminating error; fixed in script v3 (03b), which also checks every part's exact length before merging. VEP cache verified (size, SHA-256, gzip CRC); gnomAD v4.1 confirmed in its info.txt. |
