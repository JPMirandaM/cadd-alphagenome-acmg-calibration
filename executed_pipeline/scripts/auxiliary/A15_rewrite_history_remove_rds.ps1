# A15 (2026-09-27): rewrites the Git history of the repository removing every *.rds file and the manuscript/ folder
# (drafts kept outside the public repository, author 2026-09-27) from every
# commit, before the repository is made public. The removed files hold per-variant third-party scores
# (CADD, AVI, REVEL, BayesDel, AlphaMissense) or gnomAD-derived reference sets, which must not be
# redistributed (see THIRD_PARTY_NOTICES.md). git filter-branch could not be used because the Git-for-
# Windows shell failed to fork on this machine; this script uses only Git plumbing commands instead.
# Every commit is recreated with identical author, committer, dates and message; only the tree (minus
# *.rds) and the parent hashes change. The annotated/lightweight tag analysis-plan-v1 is moved to the
# rewritten commit. A mirror backup of the original repository must exist before running this script.
# Usage: powershell -File A15_rewrite_history_remove_rds.ps1   (run in C:\cadd with a clean working tree)
$ErrorActionPreference = 'Stop'
$git = "$env:LOCALAPPDATA\Programs\Git\cmd\git.exe"
Set-Location C:\cadd
[Console]::OutputEncoding = [Text.Encoding]::UTF8
$utf8 = New-Object Text.UTF8Encoding($false)
$tmpIndex = Join-Path $env:TEMP 'a15_index'
$tmpObj = Join-Path $env:TEMP 'a15_object.txt'
$map = @{}
$commits = & $git rev-list --reverse --topo-order main
foreach ($c in $commits) {
  $env:GIT_INDEX_FILE = $tmpIndex
  & $git read-tree $c
  & $git rm -r -f --cached --ignore-unmatch -q -- '*.rds' 'manuscript' | Out-Null   # -f: the temporary index need not match the working tree
  $tree = (& $git write-tree).Trim()
  Remove-Item Env:GIT_INDEX_FILE
  $raw = & $git cat-file commit $c
  $out = New-Object System.Collections.Generic.List[string]
  $inHeader = $true
  foreach ($line in $raw) {
    if ($inHeader -and $line -match '^tree ') { $out.Add("tree $tree"); continue }
    if ($inHeader -and $line -match '^parent ([0-9a-f]{40})$') { $out.Add("parent $($map[$Matches[1]])"); continue }
    if ($inHeader -and $line -eq '') { $inHeader = $false }
    $out.Add($line)
  }
  [IO.File]::WriteAllText($tmpObj, (($out -join "`n") + "`n"), $utf8)
  $new = (& $git hash-object -t commit -w $tmpObj).Trim()
  $map[$c] = $new
}
$oldHead = (& $git rev-parse main).Trim()
& $git update-ref refs/heads/main $map[$oldHead] $oldHead
# tag
$tagType = (& $git cat-file -t analysis-plan-v1).Trim()
$tagTarget = (& $git rev-parse 'analysis-plan-v1^{commit}').Trim()
if ($tagType -eq 'tag') {
  $raw = & $git cat-file tag analysis-plan-v1
  $raw = $raw | ForEach-Object { if ($_ -match '^object ') { "object $($map[$tagTarget])" } else { $_ } }
  [IO.File]::WriteAllText($tmpObj, (($raw -join "`n") + "`n"), $utf8)
  $newTag = (& $git hash-object -t tag -w $tmpObj).Trim()
  & $git update-ref refs/tags/analysis-plan-v1 $newTag
} else { & $git update-ref refs/tags/analysis-plan-v1 $map[$tagTarget] }
& $git read-tree HEAD   # index now matches the rewritten HEAD; *.rds files stay on disk, untracked
Remove-Item $tmpIndex, $tmpObj -ErrorAction SilentlyContinue
Write-Output ("Rewritten {0} commits; main {1} -> {2}; tag analysis-plan-v1 ({3}) -> commit {4}" -f $commits.Count, $oldHead.Substring(0,7), $map[$oldHead].Substring(0,7), $tagType, $map[$tagTarget].Substring(0,7))
