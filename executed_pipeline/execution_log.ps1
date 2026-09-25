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
# Step 2 [RUNNING since 2026-09-24 ~23:20] CADD v1.7 GRCh38 whole_genome_SNVs.tsv.gz
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
#    before the final render: failed replicates counted first; label position)
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
