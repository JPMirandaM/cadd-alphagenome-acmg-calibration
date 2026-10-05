# =====================================================================================
# EXECUTION LOG - commands exactly as they were run, in order.
# Project: GRCh38 SNVs with CADD v1.7 PHRED above a calibrated threshold, annotated with gene,
#          ClinVar and gnomAD; threshold calibrated following Pejaver et al. (2022).
# Machine: Windows 11 Enterprise 10.0.26200, 22 logical CPUs, 31.5 GB RAM, PowerShell 5.1.
# Scripts invoked are frozen copies in .\scripts (hashes in SHA256SUMS_scripts.txt).
# Step status: [DONE], [RUNNING] or [PENDING]. Times are the machine's local time.
# Note: scripts 01 and 03 were executed with Spanish comments/messages and are kept verbatim;
#       all later scripts are in English.
# =====================================================================================

$R = 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe'

# -------------------------------------------------------------------------------------
# Step 1 [DONE 2026-09-24] R packages (user library, no administrator rights).
# data.table and BiocManager were already installed.
# -------------------------------------------------------------------------------------
& $R -e "options(repos=c(CRAN='https://cloud.r-project.org')); BiocManager::install('Rsamtools', update=FALSE, ask=FALSE, type='binary')"
# The previous call left two dependencies uninstalled (GenomeInfoDbData: source only; BH: transient SSL error):
& $R -e "BiocManager::install(c('GenomeInfoDbData','BH'), update=FALSE, ask=FALSE)"
& $R -e "install.packages('R.utils', repos='https://cloud.r-project.org', type='binary')"
# [DONE 2026-09-25] plotting
& $R -e "install.packages('patchwork', repos='https://cloud.r-project.org', type='binary')"

# -------------------------------------------------------------------------------------
# Step 2 [DONE 2026-09-25] CADD v1.7 GRCh38 whole_genome_SNVs.tsv.gz  (started 2026-09-24 ~23:20)
# All 8 parts complete at 02:31, but script 01 stopped before merging (same curl-stderr issue as the
# VEP cache, see step 3b). Completed with script v3 (03b; verifies part lengths, then merges):
#   powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\03b_download_parallel_v3.ps1 -Parts 8   -> logs\03b_cadd_merge.log
# Verification: MD5 matches the official checksum files (whole_genome_SNVs.tsv.gz.md5, MD5SUMs,
# downloaded to C:\cadd\data\checksums\); index MD5 also matches -> logs\cadd_file_checksums.log
# 87,473,403,655 bytes; server Last-Modified 2024-01-29 12:26:37 GMT, ETag "145dd23707-61014c1cb9940".
# The first 84,180,992 bytes (part 0) come from an earlier single-connection attempt (curl -C -,
# same server and file); the rest from 8 parallel byte ranges. The script accepts only HTTP 206
# responses and checks the exact final size.
# Script 00 (single connection, run 2026-09-24 ~23:10, stopped after 84,180,992 bytes because of
# its speed) was archived on 2026-09-25 from the working copy that produced those bytes.
# -------------------------------------------------------------------------------------
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\00_download_cadd_single_connection.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\01_download_parallel.ps1 -Parts 8

# -------------------------------------------------------------------------------------
# Step 3 [DONE 2026-09-24] Annotation references (hashes in SHA256SUMS_data.txt)
# GENCODE v50 basic (released 2026-06-09); ClinVar GRCh38 VCF (fileDate=2026-09-23).
# -------------------------------------------------------------------------------------
New-Item -ItemType Directory -Force C:\cadd\annot | Out-Null
curl.exe -sS --ssl-no-revoke --fail -L -o C:\cadd\annot\gencode.v50.basic.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz
curl.exe -sS --ssl-no-revoke --fail -L -o C:\cadd\annot\clinvar.vcf.gz     https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl.exe -sS --ssl-no-revoke --fail -L -o C:\cadd\annot\clinvar.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# -------------------------------------------------------------------------------------
# Step 3b [DONE 2026-09-25 01:47] Ensembl VEP 116 GRCh38 indexed cache (27,644,657,162 bytes),
# source of gnomAD allele frequencies for all selected variants (METHODS.md, section 4.3).
# Generalised version of the download script (-Url/-Dest); same byte-range logic and validation.
# -------------------------------------------------------------------------------------
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\03_download_parallel_v2.ps1 -Parts 8 -Url https://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz -Dest C:\cadd\vep\homo_sapiens_vep_116_GRCh38.tar.gz -MinFreeGB 110
# The run above downloaded all 8 parts completely but stopped before merging: a curl retry
# warning on stderr ("curl: (56) schannel: server closed abruptly") was turned into a terminating
# error by Receive-Job under $ErrorActionPreference = 'Stop'. Script v3 (03b) tolerates stderr
# from the jobs and checks every part's exact length before merging. Rerun (no new bytes needed):
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\03b_download_parallel_v3.ps1 -Parts 8 -Url https://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz -Dest C:\cadd\vep\homo_sapiens_vep_116_GRCh38.tar.gz -MinFreeGB 110
# -> logs\03b_vep_cache_merge.log. Integrity: SHA-256 and full gzip CRC check (tar listing, 3.3 min):
tar.exe -tzf C:\cadd\vep\homo_sapiens_vep_116_GRCh38.tar.gz > C:\cadd\vep\tar_listing.txt
# Extraction of info.txt and all_vars.gz(+.csi) for chromosomes 1-22, X, Y, MT (21.37 GB, 7.4 min):
tar.exe -xzf C:\cadd\vep\homo_sapiens_vep_116_GRCh38.tar.gz homo_sapiens/116_GRCh38/info.txt
tar.exe -xzf C:\cadd\vep\homo_sapiens_vep_116_GRCh38.tar.gz homo_sapiens/116_GRCh38/<chr>/all_vars.gz homo_sapiens/116_GRCh38/<chr>/all_vars.gz.csi   # for each chr
# Lookup check (TP53 R175H vs gnomAD VCF): auxiliary\A10 -> logs\A10_test_vep_cache_lookup.log
& $R .\scripts\auxiliary\A10_test_vep_cache_lookup.R

# -------------------------------------------------------------------------------------
# Step 4 [DONE 2026-09-25] Pejaver et al. (2022) code and data
# Code: GitHub vpejaver/clingen-svi-comp_calibration, commit c56d32d4db45389bbd386c909cd49a625e3e2dc7 (MIT).
# Data: Zenodo 10.5281/zenodo.8347415, data.zip (MD5 a25359e5bbb3de9ea47dd21f03950ffb, verified).
# data.zip is extracted to the short path C:\cadd\pj because one file name exceeds the
# 260-character Windows path limit.
# -------------------------------------------------------------------------------------
curl.exe -sS --ssl-no-revoke -L --fail -o C:\cadd\calibration\repo.zip https://github.com/vpejaver/clingen-svi-comp_calibration/archive/c56d32d4db45389bbd386c909cd49a625e3e2dc7.zip
curl.exe -sS --ssl-no-revoke -L --fail -o C:\cadd\calibration\pejaver2022_data.zip "https://zenodo.org/api/records/8347415/files/data.zip/content"
tar.exe -xf C:\cadd\calibration\pejaver2022_data.zip -C C:\cadd\pj --exclude "__MACOSX" --exclude "*.DS_Store"

# -------------------------------------------------------------------------------------
# Step 5 [DONE 2026-09-25, 22.7 min] Validation of the R port: reproduce the published CADD v1.6
# calibration with the original data (B = 10,000 bootstrap replicates, 10 workers, seed 20260925).
# Outputs: C:\cadd\calibration\results\cadd16_reproduction_B10000\ ; log: logs\05_reproduce_pejaver_cadd16_B10000.log
# Result: PP3 Supporting 25.3, Moderate 28.1; BP4 Supporting 22.70, Moderate 17.33 (= published).
# -------------------------------------------------------------------------------------
& $R .\scripts\05_reproduce_pejaver_cadd16.R 10000 10

# Step 5b [DONE 2026-09-25] Figure S1 (data copied to figures\data\FigS1_cadd16_reproduction\)
Copy-Item C:\cadd\calibration\results\cadd16_reproduction_B10000\posterior_curves.tsv, C:\cadd\calibration\results\cadd16_reproduction_B10000\posterior_levels.tsv, C:\cadd\calibration\results\cadd16_reproduction_B10000\thresholds.tsv, C:\cadd\calibration\results\cadd16_reproduction_B10000\bootstrap_thresholds.tsv C:\cadd\figures\data\FigS1_cadd16_reproduction
& $R .\scripts\06_plot_posteriors.R C:/cadd/figures/data/FigS1_cadd16_reproduction C:/cadd/figures/FigS1_cadd16_reproduction "CADD v1.6 PHRED (data of Pejaver et al. 2022)"

# Step 5c [DONE 2026-09-25 ~02:40; 20.6-24.1 min per seed] -> Table S2; Figure S2 via:
#   & $R .\scripts\06b_plot_seed_sensitivity.R C:/cadd/figures/data/FigS2_seed_sensitivity C:/cadd/figures/FigS2_seed_sensitivity
#   (data copied from calibration\results\cadd16_seed_sensitivity_B10000\; figure script revised twice
#    before the final render: failed replicates counted first; label position; main-run seed
#    listed first and labelled, caption explaining the seeds)
# Originally: [RUNNING since 2026-09-25 ~01:15] Seed (Monte Carlo) sensitivity of the thresholds:
# 4 additional seeds (101, 202, 303, 404), same data and parameters, B = 10,000, 10 workers.
# The seed-20260925 results of step 5 were copied into the same folder in the same format first.
# Outputs: C:\cadd\calibration\results\cadd16_seed_sensitivity_B10000\ ; log: logs\05b_seed_sensitivity_cadd16.log
# Motivation: BP4 Strong = 0.31 in step 5 sits at ~the 1.7th percentile of valid bootstrap
# thresholds (3.4% of replicates fail the level), a region expected to vary between seeds.
& $R .\scripts\05b_seed_sensitivity_cadd16.R 101,202,303,404 10000 10

# -------------------------------------------------------------------------------------
# Auxiliary runs (scripts in .\scripts\auxiliary, logs in .\logs)
# A1-A5 were first run on 2026-09-24/25 from a temporary folder and archived verbatim;
# A1+02+A6 and A5 were re-run from the archived copies on 2026-09-25 to produce the logs.
# -------------------------------------------------------------------------------------
# [DONE] A1 build synthetic set, 02 in both modes, A6 compare  -> logs\A1_A6_synthetic_validation.log
$T = '<SCRATCH>\test'   # path hard-coded in A1 (<SCRATCH> = session temporary folder, redacted by A9)
& $R .\scripts\auxiliary\A1_make_synthetic_test_set.R
& $R .\scripts\02_rank_top_windows.R "$T\t.tsv.bgz" 20000 20 4
& $R .\scripts\02_rank_top_windows.R "$T\t.tsv.bgz" 0 20 4
& $R .\scripts\auxiliary\A6_compare_synthetic_ranking.R $T
# [DONE] A2 remote tabix test (TP53 R175H in gnomAD v4.1 exomes); A3 gnomAD joint INFO fields
& $R .\scripts\auxiliary\A2_test_remote_gnomad_tabix.R
& $R .\scripts\auxiliary\A3_inspect_gnomad_joint_fields.R
# [DONE] A4 remote gnomAD benchmark                            -> logs\A4_benchmark_gnomad_remote.log
& $R .\scripts\auxiliary\A4_benchmark_gnomad_remote.R
# [DONE] A5 ClinVar reference counts                           -> logs\A5_clinvar_reference_counts.log
& $R .\scripts\auxiliary\A5_clinvar_reference_counts.R
# [DONE 2026-09-25] A11 space-saving extractor for STORED zip archives (for the 88.5 GB AlphaGenome
# AVI zip, which requires a signed-in browser session and was downloaded manually by the author).
# Tested on a synthetic ZIP64 archive built with bsdtar in store mode (4.35 GB, 3 members incl. one
# > 4 GB): Verify passed CRC-32 for all members; Extract produced files MD5-identical to the
# originals, including the first member shifted in place. The test also confirmed that archives
# written by .NET with NoCompression (method 8) are rejected by the safety check.
# powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\auxiliary\A11_unzip_stored_inplace.ps1 -Zip <zip> -Mode List|Verify|Extract
# [DONE] A7 package inventory                                  -> package_inventory.tsv
& $R .\scripts\auxiliary\A7_record_environment.R C:\cadd\executed_pipeline

# -------------------------------------------------------------------------------------
# Step 5d [DONE 2026-09-25 04:21-04:22 UTC] Aim 1 gene-level annotations (constraint, inheritance,
# mechanism) -> C:\cadd\annot\genes\ (+ download_provenance.tsv); log: logs\07_download_aim1_gene_annotations.log
# Inspection of contents: auxiliary\A8 -> logs\A8_inspect_aim1_gene_annotations.log
# (A8 was re-run once after fixing the tab separator for the ClinGen dosage file.)
# -------------------------------------------------------------------------------------
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\07_download_aim1_gene_annotations.ps1
& $R .\scripts\auxiliary\A8_inspect_aim1_gene_annotations.R

# -------------------------------------------------------------------------------------
# Step 5e [DONE 2026-09-25 05:48 UTC] GENCODE v50 protein translations (QC of the consequence engine)
curl.exe -sS --ssl-no-revoke --fail -L -o C:\cadd\annot\gencode.v50.pc_translations.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.pc_translations.fa.gz
# Step 5f [DONE 2026-09-25] MANE Select coding positions (19,229 transcripts; 33,463,877 coding nt)
& $R .\scripts\08_build_mane_cds.R                     # -> logs\08_build_mane_cds.log
# Step 5g [DONE 2026-09-25, 6.0 min] CADD v1.7 scores of all possible SNVs at those positions
# (100,316,172 SNVs; 3 alleles at every position; none missing or duplicated)
& $R .\scripts\09_extract_cadd_coding.R 10             # -> logs\09_extract_cadd_coding.log
# Step 5h [DONE 2026-09-25, ~1.5 h] Consequence per SNV x MANE transcript + translation QC
& $R .\scripts\10_annotate_consequences.R              # -> logs\10_annotate_consequences.log
# QC: 19,194/19,229 proteins identical to GENCODE; mismatches classified by auxiliary A12
& $R .\scripts\auxiliary\A12_inspect_translation_mismatches.R   # -> logs\A12_inspect_translation_mismatches.log
# Step 5i [DONE 2026-09-25] Correction of consequences at selenocysteine, start and incomplete codons (654 records)
& $R .\scripts\10b_fix_special_codons.R                # -> logs\10b_fix_special_codons.log
# Step 5j [DONE 2026-09-25] Consequence counts; 37 non-ATG starts; ClinVar MC concordance 99.99%
& $R .\scripts\11_consequence_summary_and_clinvar_concordance.R   # -> logs\11_consequence_summary_and_clinvar_concordance.log
# Step 5k [DONE 2026-09-25] AlphaGenome AVI scores. ZIP downloaded manually by the author (signed-in
# browser session) and placed in C:\cadd\avi. SHA-256 + CRC check, then in-place extraction
# (authorised by the author; the ZIP is consumed; no backup copy kept, by the author's decision):
$Z = 'C:\cadd\avi\avi_scores_snvs_tabix.zip'
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\auxiliary\A11_unzip_stored_inplace.ps1 -Zip $Z -Mode List      # -> logs\A11_avi_zip_list.log
Get-FileHash $Z -Algorithm SHA256                                                                                               # -> logs\A11_avi_zip_sha256.log
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\auxiliary\A11_unzip_stored_inplace.ps1 -Zip $Z -Mode Verify    # -> logs\A11_avi_zip_verify.log
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\auxiliary\A11_unzip_stored_inplace.ps1 -Zip $Z -OutDir C:\cadd\avi -Mode Extract   # -> logs\A11_avi_zip_extract.log
& $R .\scripts\auxiliary\A13_inspect_avi_file.R        # format check, TP53 R175H -> logs\A13_inspect_avi_file.log
& $R .\scripts\12_extract_avi_coding.R 10              # 100,316,172 SNVs, 3.6 min -> logs\12_extract_avi_coding.log
# Analysis plan v1 frozen: commit 148d5ff, tag analysis-plan-v1 (pushed to GitHub).
# Step 5l [INTERRUPTED 2026-09-25] Reference sets, first version of script 13: stopped by Claude Code
# because the machine ran low on memory, while loading ~73 M missense SNV rows with string columns
# for all chromosomes at once. No output was written. Script revised (v2, memory-lean: per
# chromosome, disease genes filtered before aggregation, fewer columns, 6 workers); criteria unchanged.
# Rerun only after the author's approval.
# Step 5l [RERUN 2026-09-25, approved] v2 ran with peak R memory 4.5 GB but failed when parsing the
# VEP cache: (a) fread col.names given for all 42 columns while only 7 were selected; (b) the cache
# stores end = "." for SNVs, so the filter start == end removed every record. Fixed (SNVs identified
# by single-base allele_string; empty regions handled); parsing tested on TP53 before the rerun.
# Note: some alleles carry AF = 0 in the cache; per plan v1 these are not "observed" for the gnomAD
# set (AF > 0 required) and count as rare for the ClinVar set.
& $R .\scripts\13_build_reference_sets.R 6            # -> logs\13_build_reference_sets.log
# [DONE 2026-09-25, 23.1 min VEP step] Main set: 59,979 P/LP, 77,128 B/LB (4,826 genes); gnomAD set 1,955,995.

# -------------------------------------------------------------------------------------
# Step 7 CADD v1.7 and AVI calibration (analysis plan v1)
# Smoke tests (B = 20, all 18 configurations; outputs deleted) found and fixed: integer overflow in
# AUROC (n1 * n0 > 2^31) and column types in the summary. Timing with B = 200 extrapolated to
# ~30 h for 18 runs, because each replicate sorts ~2 M gnomAD scores. A faster implementation
# (04b_pejaver_port_fast.R: per-value counts instead of sorting; same sample.int draws) was
# verified bit-identical to the validated port on the Pejaver CADD v1.6 data with auxiliary A14:
# identical point posteriors, 25/25 identical replicates, and the B = 10,000 run with seed 20260925
# reproduced all 80,000 bootstrap thresholds (same NA pattern, max abs difference 0) and all final
# thresholds (6.9 vs 22.7 min). The element-wise check was first run inline (appended to the A14 log)
# and then added to the A14 script.
& $R .\scripts\auxiliary\A14_verify_fast_port.R        # -> logs\A14_verify_fast_port.log
# [DONE 2026-09-25 09:58-12:50] full calibration, 18 runs (3.6-12.0 min each), B = 10,000, 10 workers:
& $R .\scripts\14_calibrate_cadd17_avi.R 10            # -> logs\14_calibrate_cadd17_avi.log
# Note: the end-of-script summary in 14 matched "thresholds.tsv" also in bootstrap_thresholds.tsv
# (spurious rows with level = NA); the per-run outputs are unaffected. Corrected summary:
& $R .\scripts\15_summarize_calibration.R              # -> logs\15_summarize_calibration.log; tables in calibration\results\v17_avi\summary\
# Figure 1 (data copied from main\<tool>\seed20260925 to figures\data\Fig1_cadd17_avi_calibration\<tool>\):
& $R .\scripts\16_plot_figure1.R
# Paired AVI - CADD discrimination (amendment A1):
& $R .\scripts\17_paired_discrimination.R              # -> logs\17_paired_discrimination.log

# -------------------------------------------------------------------------------------
# Step 8 Aim 2c temporal hold-out (amendment A2)
# ClinVar 2024-01-27 archive (MD5 matches official 6959856d8b7a071c0942149e0e31c46f):
#   https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240127.vcf.gz (+ .md5, .tbi) -> C:\cadd\annot\clinvar_archive\
# gnomAD v4.1 AF (overall + per population) for all MANE coding SNVs (13,032,719 alleles, 13.6 min):
& $R .\scripts\13a_extract_gnomad_af_all_mane.R 6      # -> logs\13a_extract_gnomad_af_all_mane.log
# [INTERRUPTED 2026-09-25] first run of script 18 stopped by Claude Code (system low on memory; other
# applications were using ~15 GB). No output written. Script revised (v2, memory-lean: disease-gene
# flags per chromosome instead of splitting ~27 M gene strings; update joins; 8 workers). Rerun only
# after the author's approval.
# [DONE 2026-09-25, rerun approved after the author closed other applications; peak R 4.7 GB]
& $R .\scripts\18_temporal_holdout.R 8                 # -> logs\18_temporal_holdout.log; results in calibration\results\v17_avi\temporal_holdout\
# -------------------------------------------------------------------------------------
# Step 9 Aim 6 comparators (amendment A3)
# AlphaMissense (MD5 matches Zenodo) and REVEL v1.3 downloaded to C:\cadd\comparators\ (see manifest);
# BayesDel noAF bigWigs, hg38ToHg19 chain and hg19.2bit from UCSC -> C:\cadd\comparators\bayesdel\ (download_provenance.tsv)
# rtracklayer 1.66.0 installed (BiocManager, binary) to read bigWig files and lift over coordinates.
# Script 20 failed twice before completing: (1) a GNU tar earlier in R's PATH read "C:" as a remote
# host -> Windows bsdtar called by full path; (2) import.chain() does not read .gz -> chain gunzipped.
# Intermediate AlphaMissense/REVEL tables are cached between attempts. Peak R memory 14.4 GB (REVEL read).
& $R .\scripts\20_extract_comparators.R                # -> logs\20_extract_comparators.log
# [DONE 2026-09-25 18:28-20:50] 25 calibration runs (5 predictors x 5 seeds) + discrimination on the intersection:
& $R .\scripts\21_calibrate_comparators.R 10           # -> logs\21_calibrate_comparators.log; results in calibration\results\v17_avi\comparators\
# Step 10 Aim 8 ancestry (amendment A4) [DONE 2026-09-25 ~21:05]:
& $R .\scripts\22_ancestry.R                           # -> logs\22_ancestry.log; results in calibration\results\v17_avi\ancestry\
# Step 12 Amendment A6 (drift and circularity) [DONE 2026-09-26]
# ClinVar archives 2019-12-02 and 2022-11-19 (MD5 match official; SHA-256 22423f28... / f3a6cccf...) -> annot\clinvar_archive\
& $R .\scripts\24a_build_release_sets.R                # -> logs\24a_build_release_sets.log
& $R .\scripts\24b_extract_comparators_releases.R      # 393,927 new SNVs; peak R 12 GB -> logs\24b_extract_comparators_releases.log
& $R .\scripts\24c_calibrate_release_series.R 10       # 30 runs -> logs\24c_calibrate_release_series.log
& $R .\scripts\24d_threshold_shift_common_scale.R      # POST HOC common-scale comparison -> logs\24d_threshold_shift_common_scale.log
# Step 15 Aim 4 genome-wide catalogue (amendment A9)
# CADD scan [DONE 2026-09-26, ~3 h, 8 workers, peak R 11.7 GB]: 8,812,917,339 SNVs read; 24,926,994 with PHRED >= 25.6
& $R .\scripts\02_rank_top_windows.R C:/cadd/data/whole_genome_SNVs.tsv.gz 0 25.6 8   # -> logs\02_rank_cadd_ge25.6.log
# AVI scan: started with 8 workers; stopped preventively by the assistant after 21 windows because free
# memory fell to ~5 GB (no error occurred); the 21 saved windows were verified readable; resumed with 6 workers:
& $R .\scripts\28_scan_avi.R 25.18 6                   # -> logs\28_scan_avi.log  [DONE; 26,743,913 SNVs]
# Catalogue assembly: tested first on chromosome 21 with 4 workers (7.3 min, output kept); the script
# was then set to 3 workers and run for all chromosomes (296.6 min; 32,539,174 SNVs):
& $R .\scripts\29_build_catalogue.R 21                 # test run (chr21)
& $R .\scripts\29_build_catalogue.R                    # -> logs\29_build_catalogue.log; data\catalogue\
# Correction 2026-09-27: the per-gene "common" count used the overall gnomAD AF instead of the group AF
# specified in A9; recomputed (previous file kept as data\catalogue\gene_summary_v1_overallAF.tsv):
& $R .\scripts\29b_gene_summary_group_af.R             # -> logs\29b_gene_summary_group_af.log
# Step 16 Aim 9 noncoding (amendments A10, A11) [2026-09-26/27]; third attempt after fixing NA
# joins (inner joins for exon table) and skipping sequences absent from the tabix indexes (MT) ran >4 h in
# its sequential tabix lookups and was lost when the Claude Code session ended (log kept as
# logs\30_noncoding_clinvar_mpra_run3_interrupted.log). The lookups were rewritten to run per chromosome
# in 4 workers with per-chromosome caching (results unchanged by design) and the script was rerun:
& $R .\scripts\30_noncoding_clinvar_mpra.R             # -> logs\30_noncoding_clinvar_mpra.log
# Run 4 failed at the MPRA join (an SNV can be tested in more than one element; log kept as
# logs\30_noncoding_clinvar_mpra_run4_mpra_join_error.log); lookups deduplicated and rerun (run 5, completed).
# The intronic splice region was the only class with >= 1,000 P/LP and >= 1,000 B/LB (A10 rule); calibrated:
& $R .\scripts\30b_calibrate_splice_region.R           # -> logs\30b_calibrate_splice_region.log
# The MAVE label builder (33b) and the CADD v1.6 range download (35b) were also interrupted by the
# session end; 33b was rerun from scratch (downloads cached) and 35b resumed (completed parts kept).
# Step 17 Aim 8 extension, population frequencies (amendment A12) [2026-09-27]
# (a) MCPS: official per-chromosome frequency VCFs (release 20230130), 23 files + .tbi, sizes verified:
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\31_download_mcps_frequencies.ps1   # -> logs\31_download_mcps_frequencies.log; data\mcps\SHA256SUMS.txt
# First run compared with the gnomAD overall AF instead of the pre-specified group maximum (bug; log kept as
# logs\32_mcps_frequencies_run1_overallAF_bug.log); corrected and rerun:
& $R .\scripts\32_mcps_frequencies.R                   # 4 workers, ~8 min, peak R ~10 GB -> logs\32_mcps_frequencies.log
# Post hoc: MCPS SOURCE (array / exome / genome) added to mcps_only_frequency_conflicts.tsv by a one-line
# tabix lookup (Rscript -e; logged in RESULTS R14).
# (b) BRAVO (TOPMed Freeze 10): 23 public VCFs downloaded manually from signed links obtained by the author
# (links expire after 15 min and contain temporary credentials; not stored). Sizes verified against the
# server Content-Length (headers kept in data\bravo\chr*.headers.txt). Files are plain gzip; script 34
# recompresses them with bgzip and indexes them (originals kept):
& $R .\scripts\34_bravo_frequencies.R                  # 3 workers -> logs\34_bravo_frequencies.log
# Figures 3-5 [2026-09-27] (source tables copied into figures\data\<figure>\; Figure 3 legend set to two rows
# after the first render truncated it):
& $R .\scripts\36_plot_fig3_comparators.R C:/cadd/figures/data/Fig3_comparators C:/cadd/figures/Fig3_comparators
& $R .\scripts\37_plot_fig4_population.R C:/cadd/figures/data/Fig4_population C:/cadd/figures/Fig4_population
& $R .\scripts\38_plot_fig5_vus.R C:/cadd/figures/data/Fig5_vus C:/cadd/figures/Fig5_vus
# Figure revision 2026-09-27 (author's comments): titles and explanatory notes removed from the images and
# moved to the manuscript legends; panel tags A/B; Figure 3B percentages printed in segments >= 6% (a first
# attempt labelled only those segments, which misplaced the stacked labels; fixed by labelling all segments
# with empty text below 6%); Figure 4 long labels wrapped and panel spacing widened. plot_posteriors.R gained
# show_caption (default FALSE), so Figure S1 was re-rendered as well. Re-run commands:
& $R .\scripts\16_plot_figure1.R
& $R .\scripts\06_plot_posteriors.R C:/cadd/figures/data/FigS1_cadd16_reproduction C:/cadd/figures/FigS1_cadd16_reproduction "CADD v1.6 PHRED (data of Pejaver et al. 2022)"
& $R .\scripts\27_plot_gene_classes.R C:/cadd/figures/data/Fig2_gene_classes C:/cadd/figures/Fig2_gene_classes
& $R .\scripts\36_plot_fig3_comparators.R C:/cadd/figures/data/Fig3_comparators C:/cadd/figures/Fig3_comparators
& $R .\scripts\37_plot_fig4_population.R C:/cadd/figures/data/Fig4_population C:/cadd/figures/Fig4_population
& $R .\scripts\38_plot_fig5_vus.R C:/cadd/figures/data/Fig5_vus C:/cadd/figures/Fig5_vus
& $R .\scripts\40_plot_figS4_mave_genes.R C:/cadd/figures/data/FigS4_mave_genes C:/cadd/figures/FigS4_mave_genes
# Tables 1-2 (TSV + Word, landscape); source tables copied into figures\data\Table1_thresholds and Table2_mave:
& $R .\scripts\39_make_tables_1_2.R C:/cadd/figures
& $R .\scripts\40_plot_figS4_mave_genes.R C:/cadd/figures/data/FigS4_mave_genes C:/cadd/figures/FigS4_mave_genes   # Figure S4
& $R .\scripts\41_make_supplementary_tables.R            # Supplementary Tables S1-S15 (xlsx + TSV); first run failed on a fread argument, fixed
# Step 19 Aim 2b CADD v1.6 vs v1.7 (amendment A14) [2026-09-27]
& $R .\scripts\35a_cadd16_byte_ranges.R 20000          # tabix linear index -> 6,419 byte ranges (9.38 GB); also run with 2000 and 200000 to compare volumes
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\35b_download_cadd16_ranges.ps1   # 8 curl processes; interrupted by the session end, resumed; all ranges length-verified
# 35c run 1 concatenated the BGZF blocks of all ranges for tabix indexing and failed (records cut at range
# junctions; logs\35c_cadd16_assemble_extract_run1_error.log); rewritten to decompress each range separately:
& $R .\scripts\35c_cadd16_assemble_extract.R           # 938,158,116 records read; 2,093,102/2,093,102 reference SNVs scored; 28.6 min
& $R .\scripts\35d_calibrate_cadd16_vs_17.R 6          # 2 versions x 5 seeds, B = 10,000 [DONE] -> calibration\results\v17_avi\cadd16_vs_17\
# Step 18 Aim 3 MAVE validation (amendment A13) [2026-09-27]
& $R .\scripts\33a_mave_inventory.R                    # MaveDB metadata only -> logs\33a_mave_inventory.log
# 33b run 1 had two implementation errors (class-based calibrations not applied -> sets with 0 labels;
# protein QC compared whole amino-acid changes -> 36 sets wrongly dropped at ~0.30 match); log kept as
# logs\33b_mave_build_labels_run1_bugs.log. Fixed (class members from the calibration endpoint; reference-
# residue check) and rerun: 89/89 sets ok, residue match 0.992-1.000, 98,413 labelled SNVs in 39 genes.
& $R .\scripts\33b_mave_build_labels.R                 # scores + GRCh38 mappings -> logs\33b_mave_build_labels.log
& $R .\scripts\33c_mave_evaluation.R                   # B = 2,000 -> logs\33c_mave_evaluation.log; calibration\results\v17_avi\mave\
# Step 14 Aim 1 gene-class stratification (amendment A8) [DONE 2026-09-26] - 35 runs
& $R .\scripts\26_gene_class_stratification.R 10       # -> logs\26_gene_class_stratification.log; calibration\results\v17_avi\gene_classes\
& $R .\scripts\27_plot_gene_classes.R C:/cadd/figures/data/Fig2_gene_classes C:/cadd/figures/Fig2_gene_classes   # Figure 2
# Step 13 Amendment A7 (HGMD Public 2020.4 exclusion) [DONE 2026-09-26]
& $R .\scripts\25a_hgmd_public_positions.R             # 1,429,070,556 cache records scanned, 11.1 min -> logs\25a_hgmd_public_positions.log
& $R .\scripts\25b_calibrate_hgmd_excluded.R 10        # 10 runs -> logs\25b_calibrate_hgmd_excluded.log; calibration\results\v17_avi\hgmd_excluded\
# Step 11 Aim 7 VUS evidence (amendment A5) [DONE 2026-09-25]:
& $R .\scripts\23_vus_yield.R                          # -> logs\23_vus_yield.log
# Figure S3 (data copied to figures\data\FigS3_temporal_validation\; caption wrapped once after first render):
& $R .\scripts\19_plot_temporal_lr.R C:/cadd/figures/data/FigS3_temporal_validation C:/cadd/figures/FigS3_temporal_validation

# -------------------------------------------------------------------------------------
# Step 6 [PENDING: requires step 2] Selection and ranking of SNVs above the PHRED threshold
# Arguments: file, N = 0 (all above threshold), PHRED threshold, 10 workers.
# The threshold will be the one obtained from the CADD v1.7 calibration (step 7); 25 is the
# provisional value.
# -------------------------------------------------------------------------------------
# & $R .\scripts\02_rank_top_windows.R C:\cadd\data\whole_genome_SNVs.tsv.gz 0 <threshold> 10

# -------------------------------------------------------------------------------------
# Step 7 [PENDING] CADD v1.7 calibration with current ClinVar and gnomAD reference sets.
# Step 8 [PENDING] Annotation: gene (GENCODE v50), ClinVar, gnomAD.
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Repository preparation [DONE 2026-09-25]
# A9 replaces local user-profile paths with <SCRATCH> / <USER_HOME> in all versioned text files
# (the user name is read at run time and not stored in the script). Only path strings change.
# Affected: execution_log.ps1, package_inventory.tsv, logs\A1_A6_synthetic_validation.log,
# scripts\auxiliary\A1_make_synthetic_test_set.R (its SHA-256 in SHA256SUMS_scripts.txt is the
# post-redaction hash). Rerun A9 before each commit.
# -------------------------------------------------------------------------------------
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\auxiliary\A9_redact_user_paths.ps1
# Git (user install, 2.55.0.windows.5): repository initialised in C:\cadd, branch main,
# core.autocrlf false and .gitattributes "* -text" so files are stored byte-for-byte.
git init -b main
git config core.autocrlf false
git config core.longpaths true
git config user.name "Jose Patricio Miranda"
git config user.email "90576224+JPMirandaM@users.noreply.github.com"   # GitHub no-reply address
git add -A; git commit    # initial commit 3a5971f
# Private GitHub repository (created by the author): JPMirandaM/cadd-v1.7-pp3-bp4-calibration
git remote add origin https://github.com/JPMirandaM/cadd-v1.7-pp3-bp4-calibration.git
git push -u origin main   # run by the author (browser authentication via Git Credential Manager)

# -------------------------------------------------------------------------------------
# Repository release preparation [2026-09-27]
# Per-variant .rds files (third-party scores, gnomAD-derived sets) are no longer versioned (git rm --cached;
# files kept locally; .gitignore rule). They remain in earlier commits of the history. A mirror backup of the
# repository was made (C:\cadd_git_backup_2026-09-27.git). git filter-branch failed (Git Bash fork error).
# UPDATE 2026-09-27: A15 was RUN locally with the author's approval. The first pass removed nothing (git rm --cached
# refused without -f in the temporary index); fixed and rerun: 0 .rds and 0 manuscript files remain in the 62
# commits, and the analysis-plan-v1 commit keeps its original hash 3fe5c22 (it predates every removed file).
# The cleaned history is pushed to a new public repository, not force-pushed over the private one.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\auxiliary\A15_rewrite_history_remove_rds.ps1
# Original note: auxiliary\A15_rewrite_history_remove_rds.ps1 (plumbing-only rewrite) was prepared first; rewriting the
# history requires the author's explicit authorization before the repository is made public.
# The manuscript/ folder was also removed from version control (kept locally; the author makes the final edits
# manually in Word) and was added to the A15 rewrite so that drafts do not appear in the public history.
# -------------------------------------------------------------------------------------
# New repository for the cleaned history [2026-09-27] (created by the author, private until the medRxiv preprint)
# The old private repository is kept as remote "old-private" until the author deletes it.
git remote rename origin old-private
git remote add origin https://github.com/JPMirandaM/cadd-alphagenome-acmg-calibration.git
git push -u origin main; git push origin analysis-plan-v1
# CITATION.cff repository-code updated to the new URL.
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Manuscript v5 after co-author review [2026-10-04]
# Amendment A15 registered in STUDY_PLAN.md before running (VUS of 2024 resolved by 2026).
& $R .\scripts\42_vus_resolved.R                         # -> logs\42_vus_resolved.log (~70 min; RESULTS R18)
# Figure 4 (noncoding; data copied from calibration/results/v17_avi/noncoding/ to figures\data\Fig4_noncoding\):
& $R .\scripts\43_plot_fig_noncoding.R C:/cadd/figures/data/Fig4_noncoding C:/cadd/figures/Fig4_noncoding
# Table 2 (resolved VUS) and Figure 3 (VUS; data in figures\data\Fig3_vus\, former Fig5_vus):
& $R .\scripts\45_make_table2_vus_resolved.R             # -> logs\45_make_table2_vus_resolved.log
& $R .\scripts\44_plot_fig3_vus.R C:/cadd/figures/data/Fig3_vus C:/cadd/figures/Fig3_vus
# Figures renumbered for v5 (git mv): Fig3_comparators -> Fig2_comparators; Fig2_gene_classes -> Fig5_gene_classes;
# Fig4_population -> FigS5_population; Fig5_vus -> Fig3_vus (now two panels). Former Table 2 (MAVE) -> Table S17.
& $R .\scripts\41_make_supplementary_tables.R            # -> Supplementary_Tables_S1-S19.xlsx (S17-S19 added)
# Manuscript v5 and supplemental material v2 built locally (manuscript/ not versioned).
# -------------------------------------------------------------------------------------
