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
@REM CALL :add_app %KRATOS_APP_DIR%\FluidDynamicsApplication
@REM CALL :add_app %KRATOS_APP_DIR%\DEMApplication
@REM CALL :add_app %KRATOS_APP_DIR%\ContactStructuralMechanicsApplication
@REM CALL :add_app %KRATOS_APP_DIR%\ParticleMechanicsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\ConvectionDiffusionApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\DamApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\PoromechanicsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\FSIApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\SwimmingDEMApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\EigenSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
<<<<<<< HEAD
@REM CALL :add_app %KRATOS_APP_DIR%\ConstitutiveLawsApplication;
@REM @REM CALL :add_app %KRATOS_APP_DIR%\FemToDemApplication;
@REM @REM CALL :add_app %KRATOS_APP_DIR%\PfemFluidDynamicsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\DelaunayMeshingApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\MeshingApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\DemStructuresCouplingApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\MeshMovingApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\CSharpWrapperApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\ShapeOptimizationApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\CoSimulationApplication;
@REM @REM CALL :add_app %KRATOS_APP_DIR%\CableNetApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\RANSApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\MappingApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\CompressiblePotentialFlowApplication;
@REM @REM CALL :add_app %KRATOS_APP_DIR%\HDF5Application;
@REM CALL :add_app %KRATOS_APP_DIR%\IgaApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\ChimeraApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\MultilevelMonteCarloApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\StatisticsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\RomApplication;
=======
CALL :add_app %KRATOS_APP_DIR%\ConstitutiveLawsApplication;
CALL :add_app %KRATOS_APP_DIR%\FemToDemApplication;
CALL :add_app %KRATOS_APP_DIR%\PfemFluidDynamicsApplication;
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
CALL :add_app %KRATOS_APP_DIR%\ShallowWaterApplication;
>>>>>>> master

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"


echo %KRATOS_SOURCE%
echo %KRATOS_BUILD%\%KRATOS_BUILD_TYPE%

cmake -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"  ^
-DCMAKE_INSTALL_PREFIX=%3                                                                   ^
-DUSE_TRIANGLE_NONFREE_TPL=ON                                                               ^
-DCMAKE_C_FLAGS="/MP24 /Gm- /Zm10"                                                          ^
-DCMAKE_CXX_FLAGS="/MP24 /Gm- /Zm10"                                                        ^
-DBOOST_ROOT=%BOOST_ROOT%                                                                   ^
-DKRATOS_BUILD_TESTING=OFF                                                                  ^
-DINSTALL_RUNKRATOS=OFF                                                                     ^
-DKRATOS_GENERATE_PYTHON_STUBS=ON

:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
