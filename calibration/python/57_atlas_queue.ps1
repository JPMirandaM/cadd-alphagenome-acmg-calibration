# Paper 2: runs the remaining Atlas API retrievals one at a time (two concurrent jobs exceeded the per-minute
# quota on 2026-10-30). Waits for the running ClinVar batch-1 job (process id given as argument), then the
# gnomAD sample, ClinVar batch 2, and a repair pass over every output directory. Resumable: rerunning skips
# finished chunks. The API key is read from the user environment variable ALPHAGENOME_API_KEY.
param([int]$WaitPid = 0)
$env:ALPHAGENOME_API_KEY = [Environment]::GetEnvironmentVariable('ALPHAGENOME_API_KEY', 'User')
$py = 'C:\cadd\calibration\python'; $u = 'C:\cadd\calibration\results\paper2\universe'; $o = 'C:\cadd\calibration\results\paper2'
Set-Location $env:TEMP
if ($WaitPid -gt 0) { try { Wait-Process -Id $WaitPid -ErrorAction Stop } catch {} }
function Run-Batch($in, $out) { uv run --quiet --python 3.12 --with alphagenome python -I "$py\55_atlas_batch.py" $in $out --workers 32 --chunk 10000 }
Run-Batch "$u\variants_clinvar_batch1.tsv" "$o\atlas_clinvar"     # no-op if batch 1 finished; resumes it otherwise
Run-Batch "$u\variants_gnomad.tsv" "$o\atlas_gnomad"
Run-Batch "$u\variants_clinvar_batch2.tsv" "$o\atlas_clinvar_batch2"
uv run --quiet --python 3.12 --with alphagenome python -I "$py\56_atlas_repair.py" "$o\atlas_clinvar" "$o\atlas_gnomad" "$o\atlas_clinvar_batch2" --workers 8
"queue finished $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
