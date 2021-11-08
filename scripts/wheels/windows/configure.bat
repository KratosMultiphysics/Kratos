@echo off

set CC=cl.exe
set CXX=cl.exe

set KRATOS_SOURCE=%2
set KRATOS_BUILD=%KRATOS_SOURCE%/build
set KRATOS_APP_DIR=applications

set KRATOS_BUILD_TYPE=Release
set BOOST_ROOT=%BOOST%
set PYTHON_EXECUTABLE=%1

set KRATOS_APPLICATIONS=
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication
CALL :add_app %KRATOS_APP_DIR%\FluidDynamicsApplication
CALL :add_app %KRATOS_APP_DIR%\DEMApplication
CALL :add_app %KRATOS_APP_DIR%\ContactStructuralMechanicsApplication
CALL :add_app %KRATOS_APP_DIR%\ParticleMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\ConvectionDiffusionApplication;
CALL :add_app %KRATOS_APP_DIR%\DamApplication;
CALL :add_app %KRATOS_APP_DIR%\PoromechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FSIApplication;
CALL :add_app %KRATOS_APP_DIR%\SwimmingDEMApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\ExternalSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\EigenSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\ConstitutiveLawsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\FemToDemApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\PfemFluidDynamicsApplication;
CALL :add_app %KRATOS_APP_DIR%\DelaunayMeshingApplication;
CALL :add_app %KRATOS_APP_DIR%\MeshingApplication;
CALL :add_app %KRATOS_APP_DIR%\DemStructuresCouplingApplication;
CALL :add_app %KRATOS_APP_DIR%\MeshMovingApplication;
CALL :add_app %KRATOS_APP_DIR%\CSharpWrapperApplication;
CALL :add_app %KRATOS_APP_DIR%\ShapeOptimizationApplication;
CALL :add_app %KRATOS_APP_DIR%\CoSimulationApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\CableNetApplication;
CALL :add_app %KRATOS_APP_DIR%\RANSApplication;
CALL :add_app %KRATOS_APP_DIR%\MappingApplication;
CALL :add_app %KRATOS_APP_DIR%\CompressiblePotentialFlowApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\HDF5Application;
CALL :add_app %KRATOS_APP_DIR%\IgaApplication;
CALL :add_app %KRATOS_APP_DIR%\ChimeraApplication;
CALL :add_app %KRATOS_APP_DIR%\MultilevelMonteCarloApplication;
CALL :add_app %KRATOS_APP_DIR%\StatisticsApplication;
CALL :add_app %KRATOS_APP_DIR%\RomApplication;

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"


echo %KRATOS_SOURCE%
echo %KRATOS_BUILD%\%KRATOS_BUILD_TYPE%

cmake -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"  ^
-DCMAKE_INSTALL_PREFIX=%3                                                                   ^
-DCMAKE_C_FLAGS="/MP%4 /Gm- /Zm10"                                                                  ^
-DCMAKE_CXX_FLAGS="/MP%4 /Gm- /Zm10"                                                                ^
-DBOOST_ROOT=%BOOST_ROOT%                                                                   ^
-DKRATOS_BUILD_TESTING=OFF                                                                  ^
-DINSTALL_RUNKRATOS=OFF

:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
