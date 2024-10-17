call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64 || goto :error

set CC=cl.exe
set CXX=cl.exe

set KRATOS_SOURCE=%cd%
set KRATOS_BUILD=%cd%\build
set KRATOS_APP_DIR=applications

set KRATOS_APPLICATIONS=
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\CableNetApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ConvectionDiffusionApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\FluidDynamicsApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\FluidDynamicsBiomedicalApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\FluidDynamicsHydraulicsApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\FSIApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MeshingApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MeshMovingApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\LinearSolversApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\StructuralMechanicsApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\DEMApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ChimeraApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\IgaApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MPMApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MappingApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\CoSimulationApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\StatisticsApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\SwimmingDEMApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ShapeOptimizationApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ConstitutiveLawsApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\RANSApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\CompressiblePotentialFlowApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\RomApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ShallowWaterApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\GeoMechanicsApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\DamApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\PoromechanicsApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\OptimizationApplication;
@REM set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\SystemIdentificationApplication;

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

cmake                                                 ^
  -G"Visual Studio 17 2022"                           ^
  -H"%KRATOS_SOURCE%"                                 ^
  -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"              ^
  -DBOOST_ROOT="%TEMP%\boost"                         ^
  -DCMAKE_CXX_FLAGS="/Od /we4661 /we4804 /WX /wd4996 /GS-" ^
  -DKRATOS_BUILD_TESTING=OFF                          ^
  -DFORCE_LOCAL_ZLIB_COMPILATION=ON                   ^
  -DCMAKE_UNITY_BUILD=ON                                    || goto :error

cmake --build "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%" --target all_build -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64 || goto :error
cmake --build "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64 || goto :error

goto :EOF

:error
echo Failed with error #%errorlevel%.
exit /b %errorlevel%
