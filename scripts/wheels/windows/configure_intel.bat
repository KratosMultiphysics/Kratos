@echo off

set CC=cl.exe
set CXX=cl.exe

CALL "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

set KRATOS_SOURCE=%2
set KRATOS_BUILD=%KRATOS_SOURCE%/build
set KRATOS_APP_DIR=applications

set KRATOS_BUILD_TYPE=Release
set BOOST_ROOT=%BOOST%
set PYTHON_EXECUTABLE=%1

set KRATOS_APPLICATIONS=
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

echo %KRATOS_SOURCE%
echo %KRATOS_BUILD%\%KRATOS_BUILD_TYPE%

cmake -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"  ^
-DCMAKE_PREFIX_PATH="C:\Program Files (x86)\Intel\oneAPI\tbb\2021.5.2\lib\cmake\tbb"        ^
-DCMAKE_INSTALL_PREFIX=%3                                                                   ^
-DUSE_TRIANGLE_NONFREE_TPL=ON                                                               ^
-DCMAKE_C_FLAGS="/MP24 /Gm- /Zm10"                                                          ^
-DCMAKE_CXX_FLAGS="/MP24 /Gm- /Zm10"                                                        ^
-DBOOST_ROOT=%BOOST_ROOT%                                                                   ^
-DUSE_EIGEN_MKL=ON                                                                          ^
-DUSE_INTEL_TBB=ON                                                                          ^
-DKRATOS_BUILD_TESTING=OFF                                                                  ^
-DINSTALL_RUNKRATOS=OFF

:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
