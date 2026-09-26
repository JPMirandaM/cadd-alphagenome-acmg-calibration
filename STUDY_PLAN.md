# Study plan

> Agreed with the author on 2026-09-25. In scope: Aims 1–4 and 6–8 (Aims 6–8 and the temporal
> hold-out of Aim 2 added on 2026-09-25 to increase clinical impact). Aim 5 is a future proposal;
> Aim 9 (noncoding/splicing calibration) is deferred to a later phase. The literature scan
> (2026-09-25) is summarised under "Positioning"; every claim of novelty must be re-checked
> against the literature immediately before submission.

## Core analysis (prerequisite for all aims)

Calibration of **CADD v1.7 PHRED** and of the **AlphaGenome Variant Impact (AVI) score**
(AlphaGenome Atlas, released 2026-09-08) for missense variants with the ClinGen local posterior
probability method (Pejaver et al., 2022), using current ClinVar and gnomAD v4.1 reference sets,
followed by genome-wide selection of SNVs meeting each PP3 evidence level and their annotation.
Thresholds are reported as median and range over ≥ 5 random seeds.

## Planned outputs

1. **Main article**: comparative calibration (Aim 6) with gene-class stratification (Aim 1),
   temporal stability and hold-out validation (Aim 2), MAVE validation (Aim 3), VUS evidence yield
   (Aim 7) and ancestry analysis (Aim 8). Preprint (medRxiv) as soon as the core results are solid.
2. **Resource article**: genome-wide catalogue of SNVs meeting each calibrated evidence level
   (Aim 4), with ready-to-use threshold tables and VEP-compatible files.

## Pre-specified analysis plan v1 (frozen before any CADD v1.7 / AVI calibration result)

Frozen on 2026-09-25 with the Git tag `analysis-plan-v1`. Any later change must be recorded as
an amendment with date and reason, and results under the original plan must remain reported.

**Reference sets** (following Pejaver et al. 2022 as closely as possible):
1. ClinVar GRCh38 release 2026-09-23. SNVs that are missense in their GENCODE v50 MANE Select
   transcript (our validated consequence annotation); SNVs with a non-missense protein-altering
   consequence (stop gained/lost, start lost) in another overlapping MANE transcript are excluded.
2. Labels: P/LP = Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic; B/LB = Benign,
   Likely_benign, Benign/Likely_benign (exact CLNSIG values; VUS, conflicting and any combination
   with other terms excluded).
3. Review status ≥ 1 star (criteria provided).
4. Rare: maximum of gnomAD v4.1 exome and genome AF (overall) < 0.01; absent from gnomAD = rare.
5. Disease genes: genes (ClinVar GENEINFO) with ≥ 1 P/LP record of any variant type and ≥ 1 star.
6. gnomAD set (window sizing only): MANE missense SNVs observed in gnomAD v4.1 (exomes or genomes)
   with AF < 0.01, absent from ClinVar (any record), in the same disease genes.
7. Declared deviations: no HGMD exclusion (licensed resource); no training-set exclusion (neither
   CADD nor AVI was trained on ClinVar; verified); gnomAD genotype-quality (DP/GQ) and
   low-complexity-region filters not reproducible from the VEP cache.

**Scores**: CADD v1.7 PHRED (1 decimal as distributed); AVI Phred rounded to 2 decimals (the
method evaluates every distinct score value; 5-decimal AVI would multiply computation without
changing thresholds reported at 0.1 resolution).

**Calibration**: R port of Pejaver et al. (validated), α = 0.0441, c = 1124, B = 10,000,
minimum 100 ClinVar variants and 3% of the gnomAD set per window, window increment 0.1 for both
Phred scales, one-sided 95% bootstrap bound. Main analysis: 5 seeds (20260925, 101, 202, 303,
404); thresholds reported as median and range. Each sensitivity analysis: 1 seed (20260925).

**Pre-specified sensitivity analyses** (the only ones to be run):
- S1: review status ≥ 2 stars.
- S2: rare defined as AF < 0.001 (addresses the overlap between frequency-based B/LB
  classifications and AVI's frequency proxy training labels).
- S3: gene-level (cluster) bootstrap: genes are resampled with replacement instead of variants,
  for both the ClinVar and the gnomAD sets.
- S4: excluding missense SNVs within 3 nt of an internal exon boundary (possible splicing effect).

**Discrimination**: AUROC and AUPRC on the main reference set, with 95% bootstrap intervals
(1,000 replicates, variant-level).

**Amendments to plan v1** (dated; the original plan above remains the primary specification)
- A1 (2026-09-25, after the main results): paired bootstrap comparison of AUROC/AUPRC between
  AVI and CADD v1.7 (1,000 paired replicates). Reason: plan v1 specified per-predictor intervals
  only, which do not test the difference.
- A2 (2026-09-25, before running it): temporal hold-out design for Aim 2c. Calibration set =
  ClinVar release 2024-01-27 (last release before the CADD v1.7 score files, 2024-01-29/30),
  built with the plan-v1 criteria; test set = variants of the 2026-09-23 main set absent from the
  2024 set (variants whose label changed between releases are excluded and counted). Evaluation:
  interval-based positive likelihood ratios of each evidence interval on the test set (95%
  bootstrap intervals, 1,000 replicates) compared with the required LRs. AVI is evaluated with the
  same split; because AVI was released in 2026, this tests the stability of its calibration on
  newly classified variants rather than a post-release hold-out. Calibration: 1 seed (20260925),
  B = 10,000. gnomAD frequencies for all MANE coding SNVs are extracted once (script 13a) and the
  2026 reference sets are rebuilt from them as a consistency check before building the 2024 sets.

- A3 (2026-09-25, before running it): Aim 6 comparator calibration. Sources: AlphaMissense
  (Zenodo 10813168, hg38), REVEL v1.3 (official file, GRCh38 positions), BayesDel noAF (UCSC hg19
  track, lifted to GRCh38 with the UCSC chain and reference-base check; the developer site and
  Google Drive were not reachable from the local network). Where a predictor has several values
  for an SNV (transcripts/genes), the maximum is used. Score precision and window increment, chosen
  to match the original study: BayesDel rounded to 2 decimals, increment 0.01 (the BayesDel data of
  Pejaver et al. have 2 decimals); REVEL as distributed (3 decimals), increment 0.001; AlphaMissense
  as distributed (4 decimals), increment 0.001 (default of the original code); CADD and AVI as in
  plan v1 (increment 0.1). Main comparison: ClinVar and gnomAD reference variants scored by all
  five predictors (intersection), five seeds per predictor (25 runs). Discrimination (AUROC,
  AUPRC) on the intersection with paired bootstrap differences versus AVI.

- A4 (2026-09-25, before running it): Aim 8 ancestry analysis. Universe: MANE Select missense
  SNVs observed in gnomAD v4.1 (AF > 0 in exomes or genomes, VEP 116 cache). Genetic ancestry
  groups: AFR, AMR, ASJ, EAS, FIN, MID, NFE, SAS (per-group AF = maximum of exome and genome AF)
  and AMI (genomes only); "remaining" excluded. PP3 level assigned with the main-analysis median
  thresholds of CADD v1.7 and AVI (Supporting, Moderate). Outcomes, per predictor: (1) number and
  proportion of PP3-level variants with a group AF > 0.01 (proxy for a gene-agnostic BS1 level)
  and > 0.05 (BA1 level); (2) the group with the highest AF (grpmax) among those variants;
  (3) per group, the proportion of PP3-level variants observed in that group whose AF in that group
  exceeds 0.01 / 0.05, with Wilson 95% intervals, and the number of PP3-level variants that exceed
  0.01 only in that group ("group-specific common"). Descriptive analysis; no hypothesis test.
  Limitation declared: the VEP cache has no allele numbers (AN), so the ClinGen requirement of a
  minimum AN for BA1 cannot be applied; AN will be retrieved from the gnomAD joint VCF for the
  variants in conflict (stage 2 of section 4.3).

- A5 (2026-09-25, before running it): Aim 7 VUS evidence yield. Set: ClinVar 2026-09-23 SNVs
  with CLNSIG exactly "Uncertain_significance", review status ≥ 1 star, MANE Select missense
  (same annotation and exclusions as the reference sets), in disease genes; no frequency filter in
  the main count (a secondary count restricted to AF < 0.01 is also reported). Each VUS is assigned
  the PP3/BP4 level of each predictor with the median thresholds of (a) the main analysis for CADD
  v1.7 and AVI and (b) the Aim 6 intersection analysis for all five predictors. Reported: number
  and proportion of VUS per evidence level and the corresponding ACMG/AMP points (Tavtigian et
  al. 2020: Supporting ±1, Moderate ±2, Strong ±4); agreement between CADD and AVI (cross-table).
  This quantifies computational evidence only; it does not reclassify variants.

- A6 (2026-09-25, approved by the author, before running it): threshold drift and circularity
  (extends Aim 2). Motivation: with 2026 data, REVEL, BayesDel and AlphaMissense thresholds were
  less stringent than published (R6). Design:
  (a) Release series: reference sets built with the plan-v1 criteria from ClinVar GRCh38 releases
      of December 2019 (the release used by Pejaver et al.), the last release before the ClinGen
      recommendations (November 2022), 2024-01-27 and 2026-09-23; each predictor (REVEL,
      BayesDel, AlphaMissense, CADD v1.7, AVI) calibrated on the variants of each release scored
      by all five (1 seed, 20260925, B = 10,000, same gnomAD set construction per release).
  (b) Before/after split within the 2026 main set: "pre" = variants already in the November 2022
      reference set; "post" = variants absent from it. All five predictors calibrated separately
      on pre and post variants (same 2026 gnomAD set; 1 seed).
  Interpretation rule, fixed in advance: AVI (released 2026-09-08) cannot have been used as
  evidence in any ClinVar classification and serves as a negative control. A relaxation of
  thresholds on post-2022 variants that is present for REVEL/BayesDel (and CADD) but absent or
  clearly smaller for AVI supports circularity; a similar relaxation for AVI points to changes in
  ClinVar composition instead. No formal test; thresholds and differences reported descriptively.

- A7 (2026-09-26, approved by the author, before running it): HGMD exclusion sensitivity
  analysis. HGMD Professional (licensed) cannot be used; the Ensembl VEP 116 cache contains
  HGMD Public release 2020.4 records (allele_string "HGMD_MUTATION"; position only, no allele).
  (a) Count HGMD Public 2020.4 records in the cache (all chromosomes; total, unique accessions,
      unique positions, accession-prefix classes) to document its size (for reference: HGMD
      Professional 2026.4 has 584,565 variants and the current public release 354,396; author).
  (b) Exclude from the ClinVar and gnomAD reference sets every SNV located at a position of a
      single-nucleotide HGMD Public record (start = end), i.e. a conservative position-level
      exclusion, as an approximation of the original "not in HGMD 2019" filter.
  (c) Recalibrate REVEL, BayesDel, AlphaMissense, CADD v1.7 and AVI on the HGMD-excluded
      intersection sets of the ClinVar 2019-12-02 and 2026-09-23 releases (1 seed, 20260925,
      B = 10,000; precision/increments as in A3), and compare with the non-excluded results and
      with the published thresholds. Expectation fixed in advance: if HGMD inclusion explains the
      discrepancy, the 2019 HGMD-excluded REVEL/BayesDel BP4 thresholds should move towards the
      published values. Declared limitations: HGMD Public is a subset of HGMD Professional;
      position-level exclusion also removes other alleles at the same site; predictor-training
      variants are still not excluded.

- A8 (2026-09-26, before running it): Aim 1 gene-class stratification. Reference sets: the 2026
  main intersection sets (HGMD included, as in the primary analysis). Gene assignment uses the
  gene of the MANE Select transcript; SNVs assigned to more than one gene are excluded.
  Classes: (1) Constraint: gnomAD v4.1 LOEUF decile (lof.oe_ci.upper_bin_decile) of the MANE
  Select transcript: highly constrained (deciles 0–2), intermediate (3–6), tolerant (7–9).
  (2) Inheritance: gene-level mode from ClinGen gene–disease validity (Definitive, Strong or
  Moderate), and from GenCC (same classifications) when ClinGen has no such curation; "AD only"
  vs "AR only" (genes with both modes, X-linked, semidominant or other modes are excluded).
  (3) Mechanism: Gene2Phenotype molecular mechanism for records with confidence definitive,
  strong or moderate: "loss of function" vs "gain of function or dominant negative"; genes with
  both categories or only "undetermined" are excluded. For every class, the gnomAD set is
  restricted to the genes of the class. Classes with < 1,000 P/LP or < 1,000 B/LB variants are
  not calibrated. Calibration: REVEL, BayesDel, AlphaMissense, CADD v1.7 and AVI, 1 seed
  (20260925), B = 10,000, precision/increments as in A3; global prior α = 0.0441 in all classes
  (class-specific priors are a declared limitation). Pre-specified comparison: for each threshold,
  the 2.5th–97.5th percentile interval of its 10,000 bootstrap values; two classes are reported
  as different for a threshold when these intervals do not overlap. Also reported: the percentage
  of each class's reference variants whose evidence level changes when class-specific instead of
  global (main intersection, Aim 6) thresholds are used.

- A9 (2026-09-26, before running it): Aim 4 genome-wide catalogue. Scope: every SNV in the CADD
  v1.7 and AVI files (all chromosomes) with CADD PHRED ≥ 25.6 or AVI Phred ≥ 25.18 (main-analysis
  median PP3 Supporting thresholds), found by full tabix scans (script 02 for CADD; an equivalent
  scan for AVI). Each SNV carries both scores and the PP3 level of each predictor (main median
  thresholds). Annotation: MANE Select consequence (coding SNVs) or overlap with GENCODE v50 gene
  bodies (other SNVs; "intergenic" if none); ClinVar 2026-09-23 (CLNSIG, review stars, variation
  ID; exact allele match); gnomAD v4.1 exome/genome/per-group AF and grpmax from the VEP 116 cache
  (streamed over the whole cache; "not observed" when absent). A column "calibrated_context" is
  TRUE only for MANE Select missense SNVs; for all other SNVs the PP3 level is provided for
  prioritisation only and is not calibrated evidence. Per-gene summary (MANE missense): possible
  SNVs, SNVs at each PP3 level per predictor, observed in gnomAD, common (group AF > 1%), and
  ClinVar P/LP and B/LB counts. Outputs are kept outside the Git repository (third-party score
  licences; CADD non-commercial); distribution terms to be verified before any release.

**Fair comparison across predictors** (when REVEL, BayesDel and AlphaMissense are added): the main
comparison uses the intersection of variants scored by all predictors; per-predictor full-set
results are reported as secondary.

## Aims

### Aim 1 — Calibration stratified by gene class
Calibrate CADD v1.7 separately by biologically defined gene classes instead of individual genes:
(a) constraint (gnomAD v4.1 LOEUF deciles/bins), (b) mode of inheritance (autosomal dominant vs
recessive), (c) disease mechanism (loss of function / haploinsufficiency vs gain of function /
dominant negative). Test whether thresholds differ between classes beyond bootstrap and seed
variability, and quantify how many variants change evidence level when class-specific thresholds
are used.
Data: gnomAD v4.1 constraint metrics; GenCC and/or ClinGen gene–disease validity (inheritance);
ClinGen dosage sensitivity and Gene2Phenotype (mechanism). All open access [VERIFY terms].

### Aim 2 — Stability of calibrated thresholds over time and across CADD versions
Separate the effect of the predictor version from that of the reference data:
(a) same CADD version (v1.7), ClinVar snapshots from 2019 to 2026; (b) same data (current
ClinVar), CADD v1.6 vs v1.7. (c) **Temporal hold-out validation** (addresses circularity and
overfitting): thresholds derived on the ClinVar release preceding each predictor's publication
(CADD v1.7: January 2024; AVI: September 2026 [PENDING: define AVI hold-out, few post-release
variants will exist]) and interval-based likelihood ratios evaluated on variants first
classified afterwards, as done with ClinVar 2020 in the original study.
Data: archived ClinVar VCFs (NCBI); CADD v1.6 GRCh38 scores for the reference variants
[PENDING: obtain via remote tabix of the v1.6 file or full download, ~80 GB].

### Aim 3 — Independent validation with multiplexed assays of variant effect (MAVE)
Evaluate the calibrated CADD v1.7 thresholds against functional classes from saturation
mutagenesis / deep mutational scanning data that do not depend on ClinVar (e.g., BRCA1, TP53,
PTEN, MSH2 [PENDING: final gene list]), addressing the circularity of ClinVar-based calibration.
Data: MaveDB (open), ClinMAVE and published MAVE calibrations [VERIFY availability and licences].

### Aim 4 — Genome-wide resource of SNVs meeting calibrated evidence levels
Catalogue every possible SNV reaching each PP3 level of CADD v1.7, annotated with gene, ClinVar
and gnomAD v4.1, with per-gene summaries (e.g., proportion of PP3-level SNVs observed in gnomAD,
genes where many PP3-level variants are common in the population, which may indicate reduced
penetrance or miscalibration). Released as downloadable tables.
Data: outputs of the core pipeline.

### Aim 6 — Head-to-head calibration of AVI, CADD v1.7 and already-recommended predictors
Calibrate AVI and CADD v1.7 on identical current reference sets, together with predictors
already recommended by ClinGen (REVEL, BayesDel, AlphaMissense) as benchmarks, and compare the
fraction of variants reaching each evidence level and interval-based likelihood ratios. The
AVI decomposition into feature attributions (protein change vs splicing vs expression) will be
used to interpret discordances with CADD.
Circularity check (medRxiv preprint of AlphaGenome Atlas, 2026-09-16): AVI is a supervised
combiner of AlphaGenome, AlphaMissense, VEP loss-of-function annotations and two conservation
scores, trained with gnomAD v4.1 frequency proxy labels (AF > 0.001 = proxy benign; AF < 0.001 =
proxy impactful), not with ClinVar; its Phred scale is the rank among all SNVs (hg38). Because
many ClinVar B/LB variants are classified by frequency, B/LB reference variants with
0.001 ≤ AF < 0.01 overlap AVI's proxy-benign training labels and could favour AVI on the BP4 side.
Planned sensitivity analysis: repeat the calibration of all predictors restricting the ClinVar
reference set to AF < 0.001.
Data: AVI SNV scores (AlphaGenome Atlas, permissive-use download, 88.5 GB, downloaded
2026-09-25); REVEL/BayesDel/AlphaMissense scores [PENDING: source and licence — dbNSFP requires
academic registration; AlphaMissense predictions licence to be verified].

### Aim 7 — Evidence yield for variants of uncertain significance (VUS)
Apply the calibrated thresholds to ClinVar missense VUS and report how many would receive PP3 or
BP4 evidence at each strength, overall and by gene class, as a clinically interpretable measure
of potential reclassification impact (in the ACMG/AMP points framework).
Data: ClinVar (2026-09-23), core pipeline outputs.

### Aim 8 — Ancestry and population frequency of variants with computational evidence
Using gnomAD v4.1 per-population frequencies and grpmax, quantify how often SNVs reaching PP3
levels are observed, and at which frequencies, in each genetic ancestry group; identify
variants where population frequency criteria (BA1/BS1) would conflict with PP3, and assess
whether this is more frequent in under-represented groups. Prepares Aim 5.
Data: gnomAD v4.1 (VEP 116 cache; joint VCF for prioritised variants).

### Aim 9 — DEFERRED (later phase): noncoding and splicing calibration
ClinGen-style calibration of AVI and CADD v1.7 for splice-region, intronic and UTR variants.
Deferred because P/LP reference variants in these classes are scarce (see ClinVar counts,
`executed_pipeline/logs/A5_clinvar_reference_counts.log`) and levels may not be reachable; to be
considered for a second article, contrasting with the Youden-based noncoding thresholds of
Tenywa et al. (2025).

### Aim 5 — FUTURE PROPOSAL (not in current scope)
Allele frequencies of PP3-level variants in a local cohort from a population under-represented
in gnomAD, to estimate how many variants with computational evidence of pathogenicity are
population-specific and how often BA1/BS1 frequency criteria would override PP3. Requires access
to a local sequencing cohort, which is not currently available (author, 2026-09-25). The pipeline
(Aim 4 catalogue + exact allele matching) is designed so that a cohort VCF can be added later.

## Dissemination and uptake (agreed 2026-09-25)

- Preprint on medRxiv/bioRxiv as soon as core results are solid (AVI was released 2026-09-08;
  other groups are likely calibrating it).
- Ready-to-use threshold tables and VEP-compatible files; Zenodo DOI for code and results;
  `renv.lock` for the R environment.
- Optional: contact the ClinGen group (Pejaver/Radivojac; Variant Classification Working Group);
  collaboration with a clinical laboratory for diagnostic-yield evaluation; conference
  presentation (ASHG/ESHG); visual abstract.

## Target journals (to be re-evaluated with the results)

| Output | First choice | Alternatives |
|---|---|---|
| Main article | *Genetics in Medicine* | *AJHG* (if the methodological findings, e.g. threshold drift, are strong); *Genome Medicine*; *HGG Advances* |
| Resource article | *NAR Genomics and Bioinformatics* | *Human Mutation*; *NAR* Database issue (only with a maintained web portal) |

## Positioning (literature scan, 2026-09-25)

| Work | What it covers | Gap left for this study |
|---|---|---|
| Pejaver et al. 2022, AJHG | ClinGen calibration of 13 tools incl. **CADD v1.6**, ClinVar 2019 | No CADD v1.7; data from 2019 |
| Bergquist et al. 2025, Genet Med | Same method for AlphaMissense, ESM1b, VARITY (EVE considered) | CADD v1.7 not calibrated |
| Stenton et al. 2024, Genet Med | Evidence yield of calibrated PP3/BP4 per proband (Rare Genomes Project), by inheritance | Per-proband yield, not a genome-wide catalogue; not CADD v1.7 |
| Tejura et al. 2024, AJHG | Genome-wide calibration masks per-gene heterogeneity (REVEL, BayesDel) | No stratification by constraint, inheritance or mechanism; no CADD (Aim 1) |
| Chen et al. 2026, preprint (Research Square) | Gene- and domain-aware calibration (REVEL, MutPred2, AlphaMissense) | Notes but does not test constraint/inheritance/mechanism; no CADD (Aim 1) |
| Tenywa et al. 2025, NAR Genom Bioinform | Region-aware **CADD v1.7** thresholds for **noncoding** variants (Youden's J, ClinVar 2024-11) | Not mapped to ACMG/AMP evidence strengths; no coding/missense thresholds |
| MAVE–predictor comparisons (e.g., bioRxiv 2025 "Why variant effect predictors and multiplexed assays agree and disagree"; bioRxiv 2026 "Clinical evidence yield as a framework...") | Agreement between predictors and MAVEs; evidence-yield metrics | [VERIFY whether CADD v1.7 and calibrated thresholds are included] (Aim 3) |

| AlphaGenome Atlas (Google DeepMind, 2026-09-08; technical report and medRxiv 2026-09-16) | AVI score for all ~9 billion SNVs, combining AlphaGenome and AlphaMissense; benchmarks for rare-disease prioritisation | Not calibrated to ACMG/AMP evidence strengths [VERIFY in the technical report and medRxiv preprint] (Aim 6) |

No published ClinGen-style (local posterior) calibration of **CADD v1.7** or of **AVI** was found, and no study
was found that quantifies the drift of calibrated thresholds across ClinVar releases and predictor
versions (Aim 2) [VERIFY before submission].

## Sources of the literature scan
- Pejaver et al. 2022: https://www.cell.com/ajhg/fulltext/S0002-9297(22)00461-X
- Bergquist et al. 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC11429929/
- Stenton et al. 2024: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10942508/
- Tejura et al. 2024: https://pmc.ncbi.nlm.nih.gov/articles/PMC11393694/
- Chen et al. 2026 (preprint): https://pmc.ncbi.nlm.nih.gov/articles/PMC13174790/
- Tenywa et al. 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC12629842/
- MAVE/predictors: https://www.biorxiv.org/content/10.1101/2025.07.31.667868.full.pdf ; https://www.biorxiv.org/content/10.64898/2026.03.27.714777v1
- AlphaGenome Atlas: https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/ ; technical report https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/alphagenome-atlas.pdf ; medRxiv https://www.medrxiv.org/content/10.64898/2026.09.16.26363192v1.full.pdf ; Ensembl VEP integration https://www.ensembl.info/2026/09/08/alphagenome-variant-impact-scores-integrated-into-ensembl-vep/
