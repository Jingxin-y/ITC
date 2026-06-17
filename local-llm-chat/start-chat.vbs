Option Explicit

Dim shell, fso, folder, server, command
Set shell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
folder = fso.GetParentFolderName(WScript.ScriptFullName)
server = fso.BuildPath(folder, "server.ps1")

command = "powershell.exe -NoProfile -ExecutionPolicy Bypass -File " & Chr(34) & server & Chr(34)
shell.Run command, 0, False
WScript.Sleep 1200
shell.Run "http://127.0.0.1:8787", 1, False
