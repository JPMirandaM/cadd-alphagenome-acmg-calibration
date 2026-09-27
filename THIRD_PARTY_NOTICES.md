# Third-party notices

## ClinGen SVI computational calibration code (Pejaver et al., 2022)

The following files are an R translation (derivative work) of the MATLAB code in
https://github.com/vpejaver/clingen-svi-comp_calibration, commit
`c56d32d4db45389bbd386c909cd49a625e3e2dc7`, and are distributed under the original MIT license:

- `calibration/R/pejaver_port.R` and its frozen copy `executed_pipeline/scripts/04_pejaver_port.R`
- `calibration/R/plot_posteriors.R` and its frozen copy `executed_pipeline/scripts/06_plot_posteriors.R`
  (reproduces the figure layout of `plotting/plot_both_posteriors_pub.m`)

```
MIT License

Copyright (c) 2022 Vikas Pejaver

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

The data associated with that study (Zenodo 10.5281/zenodo.8347415, CC-BY-4.0) are used for
validation but are not redistributed here.

Please also cite: Pejaver V, et al. Calibration of computational tools for missense variant
pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.
*Am J Hum Genet.* 2022;109(12):2163–2177.

## Third-party data used in the analyses (not redistributed in this repository)

None of the following data files are included in the repository; scripts download or read them
from local copies. Terms must be re-verified before any data release (for example, the Aim 4
catalogue contains CADD and AVI scores). **[VERIFY each licence before deposit.]**

| Resource | Use here | Terms (as understood; verify) | Required attribution |
|---|---|---|---|
| CADD v1.7 and v1.6 GRCh38 scores | main predictor; Aim 2b | free for non-commercial use; commercial use requires a licence from the University of Washington | Schubach et al. 2024 (v1.7); Rentzsch et al. 2019/2021 (v1.6) |
| AlphaGenome Variant Impact (AVI) SNV file | predictor | AlphaGenome Atlas download terms (permissive-use section) [VERIFY] | AlphaGenome Atlas preprint (medRxiv, 2026) |
| AlphaMissense hg38 predictions (Zenodo 10813168) | comparator | CC BY-NC-SA 4.0 | Cheng et al. 2023 |
| REVEL v1.3 | comparator | free for non-commercial use [VERIFY] | Ioannidis et al. 2016 |
| BayesDel noAF (UCSC track) | comparator | per developer / UCSC terms [VERIFY] | Feng 2017 |
| ClinVar (GRCh38 VCF; 2019–2026 releases) | reference labels | public domain (NCBI) | Landrum et al. |
| gnomAD v4.1 (via Ensembl VEP 116 cache) | frequencies | ODC Open Database License / CC0 for summary data [VERIFY] | Chen et al. 2024 |
| HGMD Public 2020.4 positions (VEP cache) | sensitivity analysis | HGMD Public terms (non-commercial) [VERIFY] | Stenson et al. |
| GENCODE v50 / MANE | annotation | open | Frankish et al.; Morales et al. 2022 |
| Mexico City Prospective Study frequencies (RGC-MCPS release 20230130) | Aim 8 extension | RGC licence agreement: research use; no redistribution; no re-identification; any public use must cite the publication and state that RGC and its collaborators (University of Oxford, UNAM, INMEGEN) were the source of the data | Ziyatdinov et al. *Nature* 2023;622:784–793 |
| TOPMed BRAVO Freeze 10 public VCFs | Aim 8 extension | BRAVO data use terms accepted by the author at download [VERIFY acknowledgement text] | TOPMed / BRAVO (Taliun et al. 2021) |
| MaveDB score sets and calibrations | Aim 3 | per score set: CC0, CC BY 4.0, CC BY-SA 4.0 or CC BY-NC-SA 4.0 (listed in `calibration/results/v17_avi/mave/mavedb_human_scoreset_inventory.tsv`) | MaveDB (Rubin et al. 2025) and the primary publication of each score set |
| Kircher et al. 2019 saturation mutagenesis MPRA | Aim 9 | per repository licence [VERIFY] | Kircher et al. *Nat Commun* 2019 |
| gnomAD constraint, GenCC, ClinGen validity/dosage, Gene2Phenotype | Aim 1 | open (see `annot/genes/download_provenance.tsv`) | respective resources |

**Before making the repository public [TO DO]:** the following versioned result files contain
per-variant third-party scores (CADD, AVI, and in some files REVEL/BayesDel/AlphaMissense) and must
be reviewed against the licences above (removed from history or replaced by identifiers only if
redistribution is not allowed): `calibration/results/cadd16_reproduction_B10000/full_result.rds`,
`calibration/results/v17_avi/ancestry/pp3_level_missense_observed_in_gnomad.rds`,
`calibration/results/v17_avi/release_series/refsets_r2019.rds`, `refsets_r2022.rds`,
`calibration/results/v17_avi/temporal_holdout/refsets_2024.rds`, `test_set.rds`,
`calibration/results/v17_avi/vus/vus_scores.rds`.
