$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Port = 8787
$Prefix = "http://127.0.0.1:$Port/"
$OllamaExe = Join-Path $env:LOCALAPPDATA "Programs\Ollama\ollama.exe"
$OllamaApp = Join-Path $env:LOCALAPPDATA "Programs\Ollama\ollama app.exe"

function Test-Ollama {
    try {
        Invoke-RestMethod -Uri "http://127.0.0.1:11434/api/tags" -TimeoutSec 1 | Out-Null
        return $true
    } catch {
        return $false
    }
}

if (-not (Test-Ollama)) {
    if (Test-Path -LiteralPath $OllamaApp) {
        Start-Process -FilePath $OllamaApp -WindowStyle Hidden
    } elseif (Test-Path -LiteralPath $OllamaExe) {
        Start-Process -FilePath $OllamaExe -ArgumentList "serve" -WindowStyle Hidden
    }
    Start-Sleep -Seconds 3
}

function Get-ContentType([string]$Path) {
    switch ([IO.Path]::GetExtension($Path).ToLowerInvariant()) {
        ".html" { "text/html; charset=utf-8"; break }
        ".css" { "text/css; charset=utf-8"; break }
        ".js" { "application/javascript; charset=utf-8"; break }
        ".json" { "application/json; charset=utf-8"; break }
        ".svg" { "image/svg+xml"; break }
        default { "application/octet-stream"; break }
    }
}

$listener = [Net.Sockets.TcpListener]::new([Net.IPAddress]::Parse("127.0.0.1"), $Port)
$listener.Start()

function Write-Response($Client, [int]$StatusCode, [string]$StatusText, [string]$ContentType, [byte[]]$Body) {
    $stream = $Client.GetStream()
    $headers = @(
        "HTTP/1.1 $StatusCode $StatusText",
        "Content-Type: $ContentType",
        "Content-Length: $($Body.Length)",
        "Cache-Control: no-store",
        "Connection: close",
        "",
        ""
    ) -join "`r`n"
    $headerBytes = [Text.Encoding]::ASCII.GetBytes($headers)
    $stream.Write($headerBytes, 0, $headerBytes.Length)
    if ($Body.Length -gt 0) {
        $stream.Write($Body, 0, $Body.Length)
    }
    $stream.Flush()
}

try {
    while ($true) {
        $client = $listener.AcceptTcpClient()
        try {
            $stream = $client.GetStream()
            $reader = [IO.StreamReader]::new($stream, [Text.Encoding]::ASCII, $false, 4096, $true)
            $requestLine = $reader.ReadLine()
            if ([string]::IsNullOrWhiteSpace($requestLine)) {
                continue
            }

            while ($true) {
                $line = $reader.ReadLine()
                if ($null -eq $line -or $line.Length -eq 0) {
                    break
                }
            }

            $parts = $requestLine.Split(" ")
            $rawPath = if ($parts.Length -ge 2) { $parts[1] } else { "/" }
            $requestPath = [Uri]::UnescapeDataString($rawPath.Split("?")[0].TrimStart("/"))
            if ([string]::IsNullOrWhiteSpace($requestPath)) {
                $requestPath = "index.html"
            }

            $candidate = Join-Path $Root $requestPath
            $fullPath = [IO.Path]::GetFullPath($candidate)
            $rootFullPath = [IO.Path]::GetFullPath($Root)
            $rootPrefix = $rootFullPath.TrimEnd([IO.Path]::DirectorySeparatorChar) + [IO.Path]::DirectorySeparatorChar

            if (
                ($fullPath -ne $rootFullPath) -and
                (-not $fullPath.StartsWith($rootPrefix, [StringComparison]::OrdinalIgnoreCase) -or
                -not (Test-Path -LiteralPath $fullPath -PathType Leaf))
            ) {
                $bytes = [Text.Encoding]::UTF8.GetBytes("Not found")
                Write-Response $client 404 "Not Found" "text/plain; charset=utf-8" $bytes
            } else {
                $bytes = [IO.File]::ReadAllBytes($fullPath)
                Write-Response $client 200 "OK" (Get-ContentType $fullPath) $bytes
            }
        } catch {
            try {
                $bytes = [Text.Encoding]::UTF8.GetBytes("Server error")
                Write-Response $client 500 "Internal Server Error" "text/plain; charset=utf-8" $bytes
            } catch {
                # Ignore response failures for dropped browser connections.
            }
        } finally {
            $client.Close()
        }
    }
} finally {
    $listener.Stop()
}
