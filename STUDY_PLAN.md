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
