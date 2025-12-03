@echo off
rem Please do not modify this script

rem For any question please contact with us in:
rem  - https://github.com/KratosMultiphysics/Kratos

rem Optional parameters:
rem You can find a list will all the compiation options in INSTALL.md or here:
rem  - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

rem Set variables
if not defined KRATOS_SOURCE set KRATOS_SOURCE=%~dp0..
if not defined KRATOS_BUILD set KRATOS_BUILD=%KRATOS_SOURCE%/build

rem oneAPI call and set compiler
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
set CC=icx.exe
set CXX=icx.exe
@REM Only Ninja compiler is supported with Intel LLVM
set CMAKE_GENERATOR=Ninja


rem Warning: In windows this option only works if you run through a terminal with admin privileges
rem set KRATOS_INSTALL_PYTHON_USING_LINKS=ON

rem Set basic configuration
if not defined KRATOS_BUILD_TYPE set KRATOS_BUILD_TYPE=Release
if not defined BOOST_ROOT set BOOST_ROOT=C:\CompiledLibs\boost_1_67_0
if not defined PYTHON_EXECUTABLE set PYTHON_EXECUTABLE=C:\Windows\py.exe
if not defined NUMBER_OF_COMPILATION_CORES set NUMBER_OF_COMPILATION_CORES=%NUMBER_OF_PROCESSORS%

rem Set applications to compile
set KRATOS_APP_DIR=applications
set KRATOS_APPLICATIONS=
CALL :add_app %KRATOS_APP_DIR%\LinearSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FluidDynamicsApplication;

rem Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

rem Enable this if your build is slow and you have a multi-core machine
rem set KRATOS_PARALLEL_BUILD_FLAG=/MP4

rem Configure
@echo on
cmake                                              ^
-G "%CMAKE_GENERATOR%"                             ^
-H"%KRATOS_SOURCE%"                                ^
-B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"             ^
-DCMAKE_CXX_FLAGS=" %KRATOS_PARALLEL_BUILD_FLAG% " ^
-DCMAKE_POLICY_VERSION_MINIMUM=3.5                 ^
-DUSE_EIGEN_MKL=OFF                                ^
-DKRATOS_GENERATE_PYTHON_STUBS=ON

rem Build
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install -- -j%NUMBER_OF_COMPILATION_CORES%
goto:eof

rem Function to add apps
:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof