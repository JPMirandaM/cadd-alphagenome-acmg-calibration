# Download the Mexico City Prospective Study (MCPS) allele-frequency VCFs (RGC-MCPS browser,
# Resources page, release 20230130; exome + WGS variants of ~140,000 adults) with their tabix
# indexes. Official download links of https://rgc-mcps.regeneron.com/resources. Use is subject to
# the RGC-MCPS licence agreement and terms of use (research use; cite Ziyatdinov et al. 2023,
# Nature 622:784; acknowledge RGC, University of Oxford, UNAM and INMEGEN; no redistribution).
# Files are sequentially downloaded (one connection; no crawling) and checked against the
# Content-Length reported by the server.
$ErrorActionPreference = 'Continue'
$base = 'https://rgc-mcps.regeneron.com/downloads/20230130/'
$out  = 'C:\cadd\data\mcps'
New-Item -ItemType Directory -Force $out | Out-Null
$chrs = @(1..22 | ForEach-Object { "$_" }) + 'X'
$log = @()
foreach ($c in $chrs) {
  foreach ($suffix in '.freq.vcf.gz', '.freq.vcf.gz.tbi') {
    $name = "chr$c$suffix"; $dest = Join-Path $out $name
    $hdr = C:\Windows\System32\curl.exe --ssl-no-revoke -s -D - -o NUL -r 0-0 -m 60 ($base + $name)
    $total = [int64](([regex]::Match(($hdr -join "`n"), 'Content-Range: bytes 0-0/(\d+)')).Groups[1].Value)
    for ($try = 1; $try -le 5; $try++) {
      $have = 0; if (Test-Path $dest) { $have = (Get-Item $dest).Length }
      if ($total -gt 0 -and $have -eq $total) { break }
      C:\Windows\System32\curl.exe --ssl-no-revoke -s -S -L -C - --retry 5 --retry-delay 10 -o $dest ($base + $name)
    }
    $have = 0; if (Test-Path $dest) { $have = (Get-Item $dest).Length }
    $status = if ($have -eq $total -and $total -gt 0) { 'OK' } else { 'INCOMPLETE' }
    $line = "{0}`t{1}`t{2}`t{3}" -f $name, $total, $have, $status
    Write-Output $line; $log += $line
  }
}
$log | Set-Content -Encoding utf8 (Join-Path $out 'download_status.tsv')
Get-ChildItem $out -Filter *.freq.vcf.gz* | Get-FileHash -Algorithm SHA256 |
  ForEach-Object { '{0}  {1}' -f $_.Hash.ToLower(), (Split-Path $_.Path -Leaf) } |
  Set-Content -Encoding utf8 (Join-Path $out 'SHA256SUMS.txt')
