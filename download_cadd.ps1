# Verifica espacio y descarga CADD v1.7 GRCh38 whole_genome_SNVs (+ indice tabix) en C:\cadd\data
$ErrorActionPreference = 'Stop'

$dir     = 'C:\cadd\data'
$baseUrl = 'https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38'
$file    = 'whole_genome_SNVs.tsv.gz'
$minFreeGB = 100

New-Item -ItemType Directory -Force $dir | Out-Null
$dest = Join-Path $dir $file

# Lo ya descargado cuenta para el minimo (la descarga se reanuda)
$haveGB = 0
if (Test-Path $dest) { $haveGB = [math]::Floor((Get-Item $dest).Length / 1GB) }
$needGB = $minFreeGB - $haveGB
$freeGB = [math]::Floor((Get-PSDrive C).Free / 1GB)
Write-Host "Libre en C: $freeGB GB (se necesitan $needGB GB; ya descargados $haveGB GB)"
if ($freeGB -lt $needGB) { throw "Espacio insuficiente en C:" }

foreach ($f in @("$file.tbi", $file)) {
    $out = Join-Path $dir $f
    for ($i = 1; $i -le 50; $i++) {
        # --ssl-no-revoke: la red corporativa bloquea la consulta de revocacion (el certificado si se valida)
        $curlArgs = @('-L', '--fail', '--ssl-no-revoke', '--retry', '10', '--retry-delay', '30',
                      '--continue-at', '-', '-o', $out, "$baseUrl/$f")
        & curl.exe @curlArgs
        # 33 = el servidor rechaza el rango porque el archivo ya esta completo
        if ($LASTEXITCODE -eq 0 -or $LASTEXITCODE -eq 33) { break }
        Write-Host "curl termino con codigo $LASTEXITCODE; reintento $i en 60 s..."
        Start-Sleep 60
    }
    if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne 33) { throw "No se pudo descargar $f" }
}
Write-Host "Descarga completa: $dest ($([math]::Round((Get-Item $dest).Length / 1GB, 1)) GB)"
