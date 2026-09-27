# Aim 2b (A14), step b: download the planned byte ranges of the CADD v1.6 GRCh38 whole-genome SNV
# file (one file per range; 8 concurrent curl processes; each part verified against its expected
# length and re-requested up to 5 times).
param([string]$RangesTsv = 'C:\cadd\data\cadd16\byte_ranges.tsv', [int]$Parallel = 8)
$url = 'https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz'
$out = 'C:\cadd\data\cadd16\parts'
New-Item -ItemType Directory -Force $out | Out-Null
$ranges = Import-Csv $RangesTsv -Delimiter "`t"
function Need($r) { $f = Join-Path $out ("part_{0:D5}.bin" -f [int]$r.id); $len = [int64]$r.e - [int64]$r.s + 1
  -not ((Test-Path $f) -and (Get-Item $f).Length -eq $len) }
for ($pass = 1; $pass -le 5; $pass++) {
  $todo = @($ranges | Where-Object { Need $_ })
  Write-Output ("pass {0}: {1} ranges to fetch" -f $pass, $todo.Count)
  if ($todo.Count -eq 0) { break }
  $running = @()
  foreach ($r in $todo) {
    while (@($running | Where-Object { -not $_.HasExited }).Count -ge $Parallel) { Start-Sleep -Milliseconds 200 }
    $running = @($running | Where-Object { -not $_.HasExited })
    $f = Join-Path $out ("part_{0:D5}.bin" -f [int]$r.id)
    $running += Start-Process -FilePath C:\Windows\System32\curl.exe -NoNewWindow -PassThru `
      -ArgumentList @('--ssl-no-revoke', '-s', '-S', '--retry', '3', '-r', ("{0}-{1}" -f $r.s, $r.e), '-o', $f, $url)
  }
  $running | Wait-Process -ErrorAction SilentlyContinue
}
$bad = @($ranges | Where-Object { Need $_ })
Write-Output ("incomplete ranges after retries: {0}" -f $bad.Count)
$tot = ($ranges | Measure-Object -Property e -Sum).Sum - ($ranges | Measure-Object -Property s -Sum).Sum + $ranges.Count
Write-Output ("bytes planned: {0:N0}" -f $tot)
