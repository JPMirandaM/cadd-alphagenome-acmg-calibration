# Replaces machine-specific paths that contain the local user name with neutral placeholders,
# in all text files of the project tree that go to the repository (data folders excluded).
# Only path strings change; no code logic is modified. Every replacement is logged.
#   <SCRATCH>    = Claude Code session scratchpad (temporary folder used for tests)
#   <USER_HOME>  = user profile folder (C:\Users\<user>)
# Usage: powershell -ExecutionPolicy Bypass -File A9_redact_user_paths.ps1 [-DryRun]
param([switch]$DryRun)
$ErrorActionPreference = 'Stop'
$root = 'C:\cadd'
$log  = Join-Path $root 'executed_pipeline\logs\A9_redact_user_paths.log'

# The user folder name (long and 8.3 short form) is read at run time, so it is not stored here.
$long  = Split-Path $env:USERPROFILE -Leaf
$short = Split-Path ((New-Object -ComObject Scripting.FileSystemObject).GetFolder($env:USERPROFILE).ShortPath) -Leaf
$user  = '(' + [regex]::Escape($long) + '|' + [regex]::Escape($short) + ')'
# Order matters: the most specific pattern first. Both slash styles are matched.
$rules = @(
    @{ re = '(?i)C:[\\/]Users[\\/]' + $user + '[\\/]AppData[\\/]Local[\\/]Temp[\\/]claude[\\/]C--cadd[\\/][0-9a-f-]{36}[\\/]scratchpad'; to = '<SCRATCH>' },
    @{ re = '(?i)C:[\\/]Users[\\/]' + $user;                                                                                          to = '<USER_HOME>' }
)
$skip = '\\(data|vep|pj|annot|\.git)\\|\\calibration\\pejaver2022_repo\\'
$exts = '\.(R|ps1|sh|md|txt|tsv|log|csv|cff|yml|yaml|gitignore|gitattributes)$'

$entries = foreach ($f in Get-ChildItem $root -Recurse -File | Where-Object { $_.FullName -notmatch $skip -and $_.Name -match $exts }) {
    # Files still being written by a running job are skipped and reported (rerun later)
    try { $text = [IO.File]::ReadAllText($f.FullName) } catch { "{0}`tSKIPPED (file in use)" -f $f.FullName.Replace("$root\", ''); continue }
    $new = $text; $n = 0
    foreach ($r in $rules) { $n += ([regex]::Matches($new, $r.re)).Count; $new = [regex]::Replace($new, $r.re, $r.to) }
    if ($n -gt 0) {
        if (-not $DryRun) { [IO.File]::WriteAllText($f.FullName, $new, (New-Object Text.UTF8Encoding($false))) }
        "{0}`t{1}" -f $f.FullName.Replace("$root\", ''), $n
    }
}
$hdr = "# $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') path redaction (DryRun=$DryRun)`nfile`treplacements"
if (-not $DryRun) { ($hdr, $entries) | Out-File -Encoding utf8 $log }
$hdr; $entries
