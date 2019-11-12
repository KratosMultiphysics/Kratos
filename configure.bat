rem You can use your interpreter of choice (bash, sh, zsh, ...)

rem For any question please contact with us in:
rem  - https://github.com/KratosMultiphysics/Kratos

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
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\ContactStructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FluidDynamicsApplication;
CALL :add_app %KRATOS_APP_DIR%\MeshingApplication;
CALL :add_app %KRATOS_APP_DIR%\MeshMovingApplication;
CALL :add_app %KRATOS_APP_DIR%\DEMApplication;
CALL :add_app %KRATOS_APP_DIR%\SwimmingDEMApplication;
CALL :add_app %KRATOS_APP_DIR%\CSharpWrapperApplication;
CALL :add_app %KRATOS_APP_DIR%\SolidMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\ConstitutiveModelsApplication;
CALL :add_app %KRATOS_APP_DIR%\DelaunayMeshingApplication;
CALL :add_app %KRATOS_APP_DIR%\ContactMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\PfemApplication;
CALL :add_app %KRATOS_APP_DIR%\PfemFluidDynamicsApplication;
CALL :add_app %KRATOS_APP_DIR%\PfemSolidMechanicsApplication;

rem Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

rem Configure
@echo on
 cmake                                                                                              ^
 -G"Visual Studio 16 2019"                                                                          ^
 -H"%KRATOS_SOURCE%"                                                                                ^
 -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"                                                             ^
 -DBOOST_ROOT="C:\CompiledLibs\boost_1_67_0"                                                        ^
 -DPYTHON_EXECUTABLE="C:\Users\Kratos64\AppData\Local\Programs\Python\Python35\python.exe"          ^
 -DINCLUDE_FEAST=OFF

rem Build
rem cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install

:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof