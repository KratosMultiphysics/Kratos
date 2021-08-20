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
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FluidDynamicsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\DEMApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\ContactStructuralMechanicsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\ParticleMechanicsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\ConvectionDiffusionApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\DamApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\PoromechanicsApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\FSIApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\SwimmingDEMApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\ExternalSolversApplication;
@REM CALL :add_app %KRATOS_APP_DIR%\EigenSolversApplication;

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"


echo %KRATOS_SOURCE%
echo %KRATOS_BUILD%\%KRATOS_BUILD_TYPE%

cmake -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"  ^
-DCMAKE_INSTALL_PREFIX=%3                                                                   ^
-DCMAKE_C_FLAGS="/MP24 /Gm-"                                                                  ^
-DCMAKE_CXX_FLAGS="/MP24 /Gm-"                                                                ^
-DBOOST_ROOT=%BOOST_ROOT%                                                                   ^
-DKRATOS_BUILD_TESTING=OFF                                                                  ^
-DINSTALL_RUNKRATOS=OFF

:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
