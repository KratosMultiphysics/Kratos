:: this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
:: to use it, copy it to another file and run it

cls

@REM Set variables
if not defined KRATOS_SOURCE set KRATOS_SOURCE=%~dp0..
if not defined KRATOS_BUILD set KRATOS_BUILD=%KRATOS_SOURCE%/build

@REM Set basic configuration
if not defined KRATOS_BUILD_TYPE set KRATOS_BUILD_TYPE=Release
@REM Decomment the following in case of considering MKL
@REM if not defined MKLROOT set MKLROOT=C:\PROGRA~2\Intel\oneAPI\mkl\latest\

@REM rem setting environment variables for using intel MKL
@REM call "%MKLROOT%\env\vars.bat" intel64 lp64

:: you may want to decomment this the first time you compile
@REM Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

sh %KRATOS_BUILD%\configure_MINGW_UCRT64.sh