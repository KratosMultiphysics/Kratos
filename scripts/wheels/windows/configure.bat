@echo off

set CC=cl.exe
set CXX=cl.exe

set KRATOS_SOURCE=%2
set KRATOS_BUILD=%KRATOS_SOURCE%/build
set KRATOS_APP_DIR=applications

set KRATOS_BUILD_TYPE=Release
set BOOST_ROOT=%BOOST%
set PYTHON_EXECUTABLE=%1

set USE_COTIRE=%3

set KRATOS_APPLICATIONS=
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FluidDynamicsApplication;
CALL :add_app %KRATOS_APP_DIR%\DEMApplication;
CALL :add_app %KRATOS_APP_DIR%\ContactStructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\ParticleMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\ConvectionDiffusionApplication;
CALL :add_app %KRATOS_APP_DIR%\DamApplication;
CALL :add_app %KRATOS_APP_DIR%\PoromechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FSIApplication;
CALL :add_app %KRATOS_APP_DIR%\SwimmingDEMApplication;
CALL :add_app %KRATOS_APP_DIR%\ExternalSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\EigenSolversApplication;

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

cmake -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"          ^
-DINSTALL_EMBEDDED_PYTHON=OFF                                                                       ^
-DLAPACK_LIBRARIES=%LAPACK%                                                                         ^
-DBLAS_LIBRARIES=%BLAS%                                                                             ^
-DINSTALL_RUNKRATOS=OFF                                                                             ^
-DUSE_COTIRE=%USE_COTIRE%                                                                           ^
-DKRATOS_BUILD_TESTING=OFF

:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
