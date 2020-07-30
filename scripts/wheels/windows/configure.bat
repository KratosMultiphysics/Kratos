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
CALL :add_app %KRATOS_APP_DIR%\DEMApplication;
CALL :add_app %KRATOS_APP_DIR%\ContactStructuralMechanicsApplication;

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

%cmake% -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"        ^
-DINCLUDE_FEAST=OFF                                                                                 ^
-DINSTALL_EMBEDDED_PYTHON=OFF                                                                       ^
-DLAPACK_LIBRARIES=%LAPACK%                                                                         ^
-DBLAS_LIBRARIES=%BLAS%                                                                             ^
-DINSTALL_RUNKRATOS=OFF

:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
