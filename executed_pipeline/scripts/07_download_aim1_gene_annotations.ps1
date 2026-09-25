# Aim 1: gene-level annotation resources (constraint, inheritance, disease mechanism).
# Downloads each file to C:\cadd\annot\genes, and writes a provenance table
# (URL, download time UTC, size, SHA-256, first line) to download_provenance.tsv.
# Usage: powershell -ExecutionPolicy Bypass -File download_aim1_gene_annotations.ps1
$ErrorActionPreference = 'Stop'
$dir = 'C:\cadd\annot\genes'
New-Item -ItemType Directory -Force $dir | Out-Null
$stamp = (Get-Date).ToUniversalTime().ToString('yyyy-MM-dd')

$sources = @(
    @{ name = 'gnomad.v4.1.constraint_metrics.tsv';
       url  = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv' },
    @{ name = "gencc_submissions_$stamp.tsv";
       url  = 'https://search.thegencc.org/download/action/submissions-export-tsv' },
    @{ name = "clingen_gene_validity_$stamp.csv";
       url  = 'https://search.clinicalgenome.org/kb/gene-validity/download' },
    @{ name = "clingen_dosage_curation_list_GRCh38_$stamp.tsv";
       url  = 'https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv' },
    @{ name = "gene2phenotype_all_panels_$stamp.csv";
       url  = 'https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download/' }
)

$prov = foreach ($s in $sources) {
    $out = Join-Path $dir $s.name
    & curl.exe -sS -L --fail --ssl-no-revoke --retry 5 --retry-delay 10 -o $out $s.url
    if ($LASTEXITCODE -ne 0) { throw "Download failed ($LASTEXITCODE): $($s.url)" }
    $f = Get-Item $out
    if ($f.Length -lt 1000) { throw "File suspiciously small ($($f.Length) bytes): $($s.name)" }
    $first = (Get-Content $out -TotalCount 1)
    [pscustomobject]@{
        file         = $s.name
        url          = $s.url
        downloaded   = (Get-Date).ToUniversalTime().ToString('yyyy-MM-dd HH:mm:ss') + ' UTC'
        size_bytes   = $f.Length
        sha256       = (Get-FileHash $out -Algorithm SHA256).Hash.ToLower()
        first_line   = $first.Substring(0, [math]::Min(200, $first.Length))
    }
}
$prov | Export-Csv -Delimiter "`t" -NoTypeInformation -Encoding UTF8 (Join-Path $dir 'download_provenance.tsv')
$prov | Format-Table file, size_bytes, sha256 -AutoSize | Out-String -Width 250
