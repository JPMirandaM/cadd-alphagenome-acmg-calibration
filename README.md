# Beyond a single threshold: calibration of CADD v1.7 and AlphaGenome for ACMG/AMP computational evidence across reference sets, genes and populations

Pre-registered ClinGen (Pejaver et al. 2022) calibration of CADD v1.7 and the AlphaGenome Variant
Impact (AVI) score for PP3/BP4 evidence, with comparison against REVEL, BayesDel and AlphaMissense,
temporal hold-out, reference-set dependence (ClinVar release series, HGMD exclusion), MAVE
validation (MaveDB), gene-class calibration, VUS evidence yield, population frequencies (gnomAD
v4.1, Mexico City Prospective Study, TOPMed BRAVO) and a genome-wide catalogue of SNVs at PP3 levels
(GRCh38). Target journal for the main article: *The American Journal of Human Genetics*. Runs
natively on Windows (R 4.4.2 + Rsamtools); no WSL/Linux required.

## Directory map

| Path | Content |
|---|---|
| `METHODS.md` | Publication-style methods, figure index, references and process log (appendix) |
| `STUDY_PLAN.md` | Aims 1–4, future proposal (Aim 5), positioning against the literature |
| `manuscript/` | Manuscript drafts (v3: style-revised AJHG draft with AI-use declaration; v2: first full draft; Markdown + Word) and abstract drafts |
| `LICENSE`, `THIRD_PARTY_NOTICES.md`, `CITATION.cff` | MIT license (code), notice for code derived from Pejaver et al. (MIT), how to cite |
| `executed_pipeline/execution_log.ps1` | Every command, in the order run, with status and date |
| `executed_pipeline/scripts/` | Frozen copies of the scripts exactly as executed (`01`–`06` main pipeline) |
| `executed_pipeline/scripts/auxiliary/` | Validation and benchmarking scripts cited in the methods (`A1`–`A7`) |
| `executed_pipeline/logs/` | Console output of each run |
| `executed_pipeline/SHA256SUMS_scripts.txt` | Hashes of all frozen scripts |
| `executed_pipeline/SHA256SUMS_data.txt`, `data_manifest.tsv` | Provenance of every external data resource (URL, version, dates, size, checksum, licence, citation) |
| `executed_pipeline/R_sessionInfo.txt`, `package_inventory.tsv`, `compute_environment.txt` | Software and hardware environment, run times |
| `figures/` | Figures (PNG + PDF) and tables |
| `figures/data/<figure>/` | Data needed to redraw each figure without recomputing |
| `calibration/` | Pejaver et al. code (pinned commit), Zenodo data, R port (`R/`) and full results (`results/`) |
| `pj/` | Zenodo `data.zip` extracted to a short path (Windows 260-character path limit) |
| `annot/` | GENCODE and ClinVar reference files |
| `data/` | CADD v1.7 scores and pipeline outputs |
| `vep/` | Ensembl VEP 116 cache (gnomAD frequencies) |

Working copies of the scripts (`C:\cadd\*.R`, `*.ps1`, `calibration\R\`) may evolve; the versions
that produced reported results are always the frozen copies in `executed_pipeline/scripts/`.

## Redrawing a figure

```powershell
& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' C:\cadd\calibration\R\plot_posteriors.R `
  C:/cadd/figures/data/FigS1_cadd16_reproduction C:/cadd/figures/FigS1_cadd16_reproduction `
  "CADD v1.6 PHRED (data of Pejaver et al. 2022)"
```
Edit `plot_posteriors.R` (labels, colours, sizes) and rerun; no recomputation is needed.
