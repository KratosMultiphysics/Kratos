@echo off
rem Please do not modify this script

rem For any question please contact with us in:
rem  - https://github.com/KratosMultiphysics/Kratos

rem Optional parameters:
rem You can find a list will all the compiation options in INSTALL.md or here:
rem  - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

rem Set compiler
set CC=cl.exe
set CXX=cl.exe

rem Set variables
set KRATOS_SOURCE=~0,-1%/..
rem set KRATOS_SOURCE=..
set KRATOS_BUILD=%KRATOS_SOURCE%/build
set KRATOS_APP_DIR=applications

rem Warning: In windows this option only works if you run through a terminal with admin privileges
set KRATOS_INSTALL_PYTHON_USING_LINKS=ON

rem Set basic configuration
set KRATOS_BUILD_TYPE=Debug
set BOOST_ROOT="D:\Program_Files\boost\boost_1_70_0"
rem set PYTHON_EXECUTABLE="C:\Program Files (x86)\Microsoft Visual Studio\Shared\Python37_64\python.exe"
rem set PYTHON_LIBRARIES="C:\Program Files (x86)\Microsoft Visual Studio\Shared\Python37_64\python37.dll"
set PYTHON_EXECUTABLE=C:\Windows\py.exe


rem Set applications to compile
set KRATOS_APPLICATIONS=
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\GeoMechanicsApplication;

rem Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

rem Configure
@echo on
cmake -G"Visual Studio 16 2019" ^
-H"%KRATOS_SOURCE%" ^
-B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"          ^
-DBoost_INCLUDE_DIR="D:\Program_Files\boost\boost_1_70_0" ^
-DLAPACK_LIBRARIES="D:\Program_Files\lapack\x64\liblapack.lib"    ^
-DBLAS_LIBRARIES="D:\Program_Files\blas\x64\libblas.lib"          ^
-DUSE_EIGEN_MKL=OFF

rem Build
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64
goto:eof

rem Function to add apps
:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
