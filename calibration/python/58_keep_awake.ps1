# Paper 2: keeps Windows from sleeping on idle while the Atlas API queue (57_atlas_queue.ps1) runs.
# Uses SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED); the display may still turn off.
# Exits, releasing the request, when the queue log reports "queue finished" or no Atlas queue or python process remains.
$sig = '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);'
$k = Add-Type -MemberDefinition $sig -Name Power -Namespace KeepAwake -PassThru
$log = 'C:\cadd\calibration\results\paper2\atlas_queue.log'
while ($true) {
  [void]$k::SetThreadExecutionState([uint32]2147483649)   # 0x80000001
  $done = (Test-Path $log) -and (Select-String -Path $log -Pattern 'queue finished' -Quiet)
  $busy = Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -match 'atlas_(batch|repair|queue)' -and $_.CommandLine -notmatch 'keep_awake' }
  if ($done -or -not $busy) { break }
  Start-Sleep -Seconds 60
}
[void]$k::SetThreadExecutionState([uint32]2147483648)     # 0x80000000: clear the request
"keep-awake released $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
