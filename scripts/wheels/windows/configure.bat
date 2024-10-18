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
CALL :add_app %KRATOS_APP_DIR%\MPMApplication;
CALL :add_app %KRATOS_APP_DIR%\ConvectionDiffusionApplication;
CALL :add_app %KRATOS_APP_DIR%\DamApplication;
CALL :add_app %KRATOS_APP_DIR%\PoromechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FSIApplication;
CALL :add_app %KRATOS_APP_DIR%\SwimmingDEMApplication;
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
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
CALL :add_app %KRATOS_APP_DIR%\CableNetApplication;
CALL :add_app %KRATOS_APP_DIR%\RANSApplication;
CALL :add_app %KRATOS_APP_DIR%\MappingApplication;
CALL :add_app %KRATOS_APP_DIR%\CompressiblePotentialFlowApplication;
CALL :add_app %KRATOS_APP_DIR%\HDF5Application;
CALL :add_app %KRATOS_APP_DIR%\MedApplication;
CALL :add_app %KRATOS_APP_DIR%\IgaApplication;
CALL :add_app %KRATOS_APP_DIR%\ChimeraApplication;
CALL :add_app %KRATOS_APP_DIR%\StatisticsApplication;
CALL :add_app %KRATOS_APP_DIR%\RomApplication;
CALL :add_app %KRATOS_APP_DIR%\ShallowWaterApplication;
CALL :add_app %KRATOS_APP_DIR%\OptimizationApplication;
CALL :add_app %KRATOS_APP_DIR%\GeoMechanicsApplication;

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
-DKRATOS_BUILD_TESTING=ON                                                                   ^
-DHDF5_ROOT="c:\hdf5\bin"                                                                   ^
-DMED_ROOT="c:\med\bin"                                                                     ^
-DKRATOS_GENERATE_PYTHON_STUBS=ON

:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
