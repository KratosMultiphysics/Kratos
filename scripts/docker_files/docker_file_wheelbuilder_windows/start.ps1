param([string]$branch = "master",
    [string]$cotire = "OFF",
    [string[]]$pythons = @("38", "37", "36", "35"),
    [int]$cpus = (Get-WmiObject -class Win32_ComputerSystem).numberoflogicalprocessors)

Write-Host "Starting build for branch $( $branch )"
Write-Host "Using $($cpus) cpus"

#Load development env
cmd.exe /c "call `"C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat`" && set > %temp%\vcvars.txt"

Get-Content "$env:temp\vcvars.txt" | Foreach-Object {
    if ($_ -match "^(.*?)=(.*)$") {
        Set-Content "env:\$( $matches[1] )" $matches[2]
    }
}

mkdir c:\kratos
mkdir c:\scripts

cd c:\kratos
git clone --depth 1 --single-branch --branch $branch https://github.com/KratosMultiphysics/Kratos.git kratos

cp c:\kratos\Kratos\scripts\wheels\windows\build.ps1 c:\scripts\build.ps1
& "c:\scripts\build.ps1" -cotire $cotire -pythons $pythons -cpus $cpus