# Extracts a large ZIP archive whose members are STORED (method 0, no compression, e.g. already
# bgzip-compressed tabix files) using almost no extra disk space:
#   - List   : parse the central directory (ZIP64 aware) and print members; saves <zip>.members.tsv
#   - Verify : CRC-32 of every member's bytes inside the ZIP vs the central directory (read-only)
#   - Extract: members are processed from the last to the first by position; each is copied out
#              and the ZIP is then truncated before it, so space is never duplicated; the first
#              member is shifted to the start of the file in place and the file is renamed.
#              Every extracted file is CRC-checked against the central directory.
# Extract DESTROYS the ZIP (by design). Run Verify first.
# Usage: powershell -ExecutionPolicy Bypass -File A11_unzip_stored_inplace.ps1 -Zip <file> -OutDir <dir> -Mode List|Verify|Extract
param(
    [Parameter(Mandatory)] [string]$Zip,
    [string]$OutDir = (Split-Path $Zip),
    [ValidateSet('List', 'Verify', 'Extract')] [string]$Mode = 'List',
    [int64]$MarginBytes = 2GB
)
$ErrorActionPreference = 'Stop'

Add-Type -TypeDefinition @'
using System;
using System.IO;
using System.Collections.Generic;
public class ZipMember {
    public string Name; public int Method; public uint Crc; public long CompSize; public long Size;
    public long HeaderOffset; public long DataOffset;
}
public static class StoredZip {
    static uint[] table = MakeTable();
    static uint[] MakeTable() {
        var t = new uint[256];
        for (uint i = 0; i < 256; i++) { uint c = i; for (int k = 0; k < 8; k++) c = (c & 1) != 0 ? 0xEDB88320u ^ (c >> 1) : c >> 1; t[i] = c; }
        return t;
    }
    public static uint Crc32(string path, long offset, long length) {
        uint c = 0xFFFFFFFFu; var buf = new byte[16 << 20];
        using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20)) {
            fs.Seek(offset, SeekOrigin.Begin); long left = length;
            while (left > 0) {
                int n = fs.Read(buf, 0, (int)Math.Min(buf.Length, left)); if (n <= 0) throw new IOException("unexpected EOF");
                for (int i = 0; i < n; i++) c = table[(c ^ buf[i]) & 0xFF] ^ (c >> 8);
                left -= n;
            }
        }
        return c ^ 0xFFFFFFFFu;
    }
    static long U(byte[] b, int o, int n) { long v = 0; for (int i = n - 1; i >= 0; i--) v = (v << 8) | b[o + i]; return v; }
    public static List<ZipMember> Read(string path) {
        using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read)) {
            long len = fs.Length; int tail = (int)Math.Min(len, 65557 + 20 + 56);
            var t = new byte[tail]; fs.Seek(len - tail, SeekOrigin.Begin); fs.Read(t, 0, tail);
            int e = -1; for (int i = tail - 22; i >= 0; i--) if (t[i] == 0x50 && t[i+1] == 0x4b && t[i+2] == 5 && t[i+3] == 6) { e = i; break; }
            if (e < 0) throw new Exception("End of central directory not found");
            long cdCount = U(t, e + 10, 2), cdSize = U(t, e + 12, 4), cdOff = U(t, e + 16, 4);
            int loc = e - 20;
            if (loc >= 0 && U(t, loc, 4) == 0x07064b50) {           // ZIP64 locator
                long z64 = U(t, loc + 8, 8); var r = new byte[56]; fs.Seek(z64, SeekOrigin.Begin); fs.Read(r, 0, 56);
                if (U(r, 0, 4) != 0x06064b50) throw new Exception("Bad ZIP64 EOCD");
                cdCount = U(r, 32, 8); cdSize = U(r, 40, 8); cdOff = U(r, 48, 8);
            }
            var cd = new byte[cdSize]; fs.Seek(cdOff, SeekOrigin.Begin); fs.Read(cd, 0, (int)cdSize);
            var list = new List<ZipMember>(); int p = 0;
            for (long k = 0; k < cdCount; k++) {
                if (U(cd, p, 4) != 0x02014b50) throw new Exception("Bad central directory entry");
                var m = new ZipMember();
                m.Method = (int)U(cd, p + 10, 2); m.Crc = (uint)U(cd, p + 16, 4);
                m.CompSize = U(cd, p + 20, 4); m.Size = U(cd, p + 24, 4);
                int nl = (int)U(cd, p + 28, 2), el = (int)U(cd, p + 30, 2), cl = (int)U(cd, p + 32, 2);
                m.HeaderOffset = U(cd, p + 42, 4);
                m.Name = System.Text.Encoding.UTF8.GetString(cd, p + 46, nl);
                int x = p + 46 + nl, xe = x + el;
                while (x + 4 <= xe) {                                  // ZIP64 extra field
                    int id = (int)U(cd, x, 2), sz = (int)U(cd, x + 2, 2), q = x + 4;
                    if (id == 1) {
                        if (m.Size == 0xFFFFFFFF) { m.Size = U(cd, q, 8); q += 8; }
                        if (m.CompSize == 0xFFFFFFFF) { m.CompSize = U(cd, q, 8); q += 8; }
                        if (m.HeaderOffset == 0xFFFFFFFF) { m.HeaderOffset = U(cd, q, 8); q += 8; }
                    }
                    x += 4 + sz;
                }
                var lh = new byte[30]; fs.Seek(m.HeaderOffset, SeekOrigin.Begin); fs.Read(lh, 0, 30);
                if (U(lh, 0, 4) != 0x04034b50) throw new Exception("Bad local header for " + m.Name);
                m.DataOffset = m.HeaderOffset + 30 + U(lh, 26, 2) + U(lh, 28, 2);
                list.Add(m); p += 46 + nl + el + cl;
            }
            return list;
        }
    }
    public static void CopyOut(string src, long off, long len, string dst) {
        var buf = new byte[16 << 20];
        using (var i = new FileStream(src, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20))
        using (var o = new FileStream(dst, FileMode.CreateNew, FileAccess.Write, FileShare.None, 1 << 20)) {
            i.Seek(off, SeekOrigin.Begin); long left = len;
            while (left > 0) { int n = i.Read(buf, 0, (int)Math.Min(buf.Length, left)); if (n <= 0) throw new IOException("EOF"); o.Write(buf, 0, n); left -= n; }
        }
    }
    public static void Truncate(string path, long newLen) {
        using (var fs = new FileStream(path, FileMode.Open, FileAccess.Write)) fs.SetLength(newLen);
    }
    // Moves [off, off+len) to [0, len) inside the same file (safe: destination precedes source), then truncates
    public static void ShiftToStart(string path, long off, long len) {
        var buf = new byte[16 << 20];
        using (var fs = new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.None, 1 << 20)) {
            long done = 0;
            while (done < len) {
                int n = (int)Math.Min(buf.Length, len - done);
                fs.Seek(off + done, SeekOrigin.Begin); int r = 0; while (r < n) { int k = fs.Read(buf, r, n - r); if (k <= 0) throw new IOException("EOF"); r += k; }
                fs.Seek(done, SeekOrigin.Begin); fs.Write(buf, 0, n); done += n;
            }
            fs.SetLength(len);
        }
    }
}
'@

$members = [StoredZip]::Read($Zip)
$tsv = "$Zip.members.tsv"
$members | Select-Object Name, Method, @{n='Crc';e={'{0:x8}' -f $_.Crc}}, CompSize, Size, HeaderOffset, DataOffset |
    Export-Csv -Delimiter "`t" -NoTypeInformation -Encoding UTF8 $tsv
$members | Format-Table Name, Method, @{n='GB';e={[math]::Round($_.Size/1GB,2)}}, @{n='Crc';e={'{0:x8}' -f $_.Crc}} -AutoSize | Out-String -Width 200
"Members: $($members.Count); total {0:N2} GB; methods: {1}" -f (($members | Measure-Object Size -Sum).Sum/1GB), (($members.Method | Sort-Object -Unique) -join ',')
if ($Mode -eq 'List') { return }

if ($members | Where-Object { $_.Method -ne 0 -or $_.CompSize -ne $_.Size }) { throw "Some members are compressed (method != 0); this script handles STORED members only." }

if ($Mode -eq 'Verify') {
    $bad = 0
    foreach ($m in $members) {
        $c = [StoredZip]::Crc32($Zip, $m.DataOffset, $m.Size)
        $ok = ($c -eq $m.Crc); if (-not $ok) { $bad++ }
        "{0,-60} CRC {1:x8} expected {2:x8} {3}" -f $m.Name, $c, $m.Crc, ($(if ($ok) { 'OK' } else { 'MISMATCH' }))
    }
    if ($bad) { throw "$bad member(s) failed CRC" } else { "All members passed CRC." }
    return
}

# Extract
New-Item -ItemType Directory -Force $OutDir | Out-Null
$order = $members | Sort-Object HeaderOffset -Descending
$drive = (Get-Item $OutDir).PSDrive.Name
for ($i = 0; $i -lt $order.Count; $i++) {
    $m = $order[$i]; $dst = Join-Path $OutDir ($m.Name -replace '/', '\')
    New-Item -ItemType Directory -Force (Split-Path $dst) | Out-Null
    if (Test-Path $dst) { throw "Destination exists: $dst" }
    if ($i -lt $order.Count - 1) {
        $free = (Get-PSDrive $drive).Free
        if ($free -lt $m.Size + $MarginBytes) { throw ("Not enough space for {0}: need {1:N1} GB, free {2:N1} GB" -f $m.Name, (($m.Size + $MarginBytes)/1GB), ($free/1GB)) }
        [StoredZip]::CopyOut($Zip, $m.DataOffset, $m.Size, $dst)
        [StoredZip]::Truncate($Zip, $m.HeaderOffset)
    } else {
        [StoredZip]::ShiftToStart($Zip, $m.DataOffset, $m.Size)
        Move-Item $Zip $dst
    }
    $c = [StoredZip]::Crc32($dst, 0, $m.Size)
    if ($c -ne $m.Crc) { throw ("CRC mismatch after extraction: {0} ({1:x8} vs {2:x8})" -f $m.Name, $c, $m.Crc) }
    "{0} extracted, {1:N2} GB, CRC {2:x8} OK" -f $m.Name, ($m.Size/1GB), $c
}
"Done: $($order.Count) members extracted to $OutDir; members list kept in $tsv"
