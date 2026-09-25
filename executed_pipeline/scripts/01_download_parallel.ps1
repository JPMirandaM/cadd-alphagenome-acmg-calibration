# Descarga en paralelo por rangos de bytes (curl.exe), reanudable, y une las partes.
# Uso: powershell -ExecutionPolicy Bypass -File C:\cadd\download_parallel.ps1 [-Parts 8]
# No cambies -Parts a mitad de una descarga: los rangos de cada parte dependen de ese numero.
param([int]$Parts = 8)
$ErrorActionPreference = 'Stop'

$dir     = 'C:\cadd\data'
$url     = 'https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz'
$dest    = Join-Path $dir 'whole_genome_SNVs.tsv.gz'
$partDir = Join-Path $dir 'parts_dl'
$minFreeGB = 100
New-Item -ItemType Directory -Force $partDir | Out-Null

# Tamano total segun el servidor
$head = & curl.exe -sIL --ssl-no-revoke $url
$total = [int64](($head | Select-String 'Content-Length:\s*(\d+)' | Select-Object -Last 1).Matches[0].Groups[1].Value)
Write-Host "Tamano total: $total bytes"

# Si hay una descarga previa de una sola conexion (incompleta), se usa como inicio de la parte 0
if ((Test-Path $dest) -and (Get-Item $dest).Length -lt $total -and -not (Test-Path "$partDir\p00.000")) {
    Move-Item $dest "$partDir\p00.000"
}
if ((Test-Path $dest) -and (Get-Item $dest).Length -eq $total) { Write-Host "Ya estaba completo."; exit 0 }

$haveGB = [math]::Floor(((Get-ChildItem $partDir -File | Measure-Object Length -Sum).Sum) / 1GB)
$freeGB = [math]::Floor((Get-PSDrive C).Free / 1GB)
# Se necesita espacio para las partes y, al unir, para una copia completa
Write-Host "Libre en C: $freeGB GB (minimo $minFreeGB GB; ya descargados $haveGB GB)"
if ($freeGB -lt ($minFreeGB - $haveGB)) { throw "Espacio insuficiente en C:" }

$chunk = [math]::Ceiling($total / $Parts)
$worker = {
    param($k, $start, $end, $url, $partDir)
    $len = $end - $start + 1
    for ($try = 1; $try -le 200; $try++) {
        $files = @(Get-ChildItem $partDir -Filter ("p{0:D2}.*" -f $k) | Sort-Object Name)
        $have = [int64]($files | Measure-Object Length -Sum).Sum
        if ($have -ge $len) { return "parte $k completa ($have bytes)" }
        $out = Join-Path $partDir ("p{0:D2}.{1:D3}" -f $k, $files.Count)
        # --speed-limit/--speed-time: si se estanca 2 min, corta y el bucle reanuda
        $code = & curl.exe -sS -L --ssl-no-revoke --fail --speed-limit 1024 --speed-time 120 -r "$($start + $have)-$end" -o $out -w '%{http_code}' $url
        # Solo sirve una respuesta parcial (206); cualquier otra cosa se descarta
        if ($code -ne '206' -and (Test-Path $out)) { Remove-Item $out -Force }
        if ((Test-Path $out) -and (Get-Item $out).Length -eq 0) { Remove-Item $out -Force }
        if ($LASTEXITCODE -ne 0) { Start-Sleep 30 }
    }
    throw "parte $k no se completo"
}

$jobs = for ($k = 0; $k -lt $Parts; $k++) {
    $start = [int64]$k * $chunk
    $end   = [math]::Min($start + $chunk - 1, $total - 1)
    Start-Job -ScriptBlock $worker -ArgumentList $k, $start, $end, $url, $partDir
}

# Progreso cada minuto
while ($jobs | Where-Object State -eq 'Running') {
    Start-Sleep 60
    # El listado de directorio no actualiza el tamano de archivos abiertos; se abre cada uno para leerlo
    $got = [int64]0
    foreach ($f in (Get-ChildItem $partDir -File)) {
        try {
            $h = [System.IO.File]::Open($f.FullName, 'Open', 'Read', 'ReadWrite')
            $got += $h.Length; $h.Dispose()
        } catch { $got += $f.Length }
    }
    Write-Host ("[{0:HH:mm}] {1:N2} / {2:N2} GB ({3:P1})" -f (Get-Date), ($got / 1GB), ($total / 1GB), ($got / $total))
}
$jobs | Receive-Job
if ($jobs | Where-Object State -ne 'Completed') { throw "Alguna parte fallo; vuelve a ejecutar el script para reanudar." }

# Unir en orden sobre la primera parte, borrando cada parte al agregarla (no duplica 81 GB)
$all = @(Get-ChildItem $partDir -File | Sort-Object Name)
$tmp = "$dest.joining"
Move-Item $all[0].FullName $tmp
$fs = [System.IO.File]::Open($tmp, [System.IO.FileMode]::Append)
try {
    foreach ($f in ($all | Select-Object -Skip 1)) {
        $in = [System.IO.File]::OpenRead($f.FullName)
        try { $in.CopyTo($fs, 16MB) } finally { $in.Dispose() }
        Remove-Item $f.FullName -Force
    }
} finally { $fs.Dispose() }
if ((Get-Item $tmp).Length -ne $total) { throw "Tamano final incorrecto: $((Get-Item $tmp).Length) != $total" }
Move-Item $tmp $dest
Remove-Item -Recurse -Force $partDir
Write-Host "Descarga completa y verificada: $dest ($total bytes)"
