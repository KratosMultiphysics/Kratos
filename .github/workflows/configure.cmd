call "%ProgramFiles(x86)%\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64 || goto :error

set CC=cl.exe
set CXX=cl.exe

set KRATOS_SOURCE=%cd%
set KRATOS_BUILD=%cd%\build
set KRATOS_APP_DIR=applications

set KRATOS_APPLICATIONS=
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ConvectionDiffusionApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\FluidDynamicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\FSIApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MeshingApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MeshMovingApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\LinearSolversApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\StructuralMechanicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\DEMApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ChimeraApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\IgaApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ParticleMechanicsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\MappingApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\CoSimulationApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\StatisticsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\SwimmingDEMApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ShapeOptimizationApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\ConstitutiveLawsApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\RANSApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\CompressiblePotentialFlowApplication;
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%KRATOS_APP_DIR%\RomApplication;

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

cmake                                                 ^
  -G"Visual Studio 16 2019"                           ^
  -H"%KRATOS_SOURCE%"                                 ^
  -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"              ^
  -DBOOST_ROOT="%TEMP%\boost"                         ^
  -DINSTALL_RUNKRATOS=OFF                             ^
  -DCMAKE_CXX_FLAGS="/Od /we4661 /we4804 /WX /wd4996" ^
  -DFORCE_LOCAL_ZLIB_COMPILATION=ON                   ^
  -DCMAKE_UNITY_BUILD=ON                                    || goto :error

cmake --build "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%" --target all_build -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64 || goto :error
cmake --build "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64 || goto :error

goto :EOF

:error
echo Failed with error #%errorlevel%.
exit /b %errorlevel%
