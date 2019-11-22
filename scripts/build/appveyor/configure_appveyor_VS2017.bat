rem For any question please contact with us in:
rem  - https://github.com/KratosMultiphysics/Kratos

rem Optional parameters:
rem You can find a list will all the compiation options in INSTALL.md or here:
rem   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

rem Set compiler
@echo off
set CC=cl.exe
set CXX=cl.exe

rem Set variables
set KRATOS_BUILD_TYPE=Release
set KRATOS_SOURCE=.
set KRATOS_BUILD=.\build
set KRATOS_APP_DIR=applications

set KRATOS_APPLICATIONS=
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\FluidDynamicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\StructuralMechanicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ContactStructuralMechanicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MeshingApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MeshMovingApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\DEMApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\SwimmingDEMApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\CSharpWrapperApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\SolidMechanicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ConstitutiveModelsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\DelaunayMeshingApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ContactMechanicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\PfemApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\PfemFluidDynamicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\PfemSolidMechanicsApplication;

rem Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

rem Configur
@echo on
 cmake                                              ^
 -G"Visual Studio 15 2017 Win64"                    ^
 -H"%KRATOS_SOURCE%"                                ^
 -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"             ^
 -DBOOST_ROOT="C:\Libraries\boost_1_65_1"           ^
 -DPYTHON_EXECUTABLE="C:\Python36-x64\python.exe"   ^
 -DINCLUDE_FEAST=OFF                                ^
 -DINSTALL_RUNKRATOS=OFF                            ^
 -DUSE_COTIRE=ON

rem Build
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target all_unity