# Contents
* [Basic Configuration](#basic-configuration)
* [Adding Applications](#basic-configuration)
* [Adding Kratos to Path](#adding-kratos-to-path)
* [Examples](#examples)
  * [Linux](#linux)
  * [Windows](#windows)
  * [MacOS](#macos)
* [Advanced Configuration](#advanced-configuration)
  * [Building Environment](#building-environments)
  * [Common Flags](#common-flags)
  * [Compilation Performance](#compilation-performance)
  * [Parallelism](#parallelism)
  * [External Libraries](#external-libraries)
    * [Feast](#feast)
    * [Metis](#metis)
    * [Trilinos](#trilinos)
* [Applications](#applications)


## Basic Configuration

You can find the new kratos configuration file in Kratos `scripts` folder: `standard_configure.sh` for linux, `standard_configure_max.sh` for MacOS, `standard_configure.bat` for win and others.

Out of the box Kratos will try to find all necessary libraries in your system automatically, but we recommend you to copy these scripts and modify it according to your preferences. Please take a look at the following configuration options:

`KRATOS_BUILD_TYPE`

Compilation Type. Options are `Release`,`RelWithDebInfo`,`Debug`,`FullDebug`,`Custom`

**Release**: Full Release with maximum optimization options and no debug Info.

**RelWithDebInfo**: Full Release with optimization and debug info. Adecuate to debug simple problems without losing performance.

**Debug**: Debug build with no optimization flags.

**FullDebug**: Debug build with no optimization falgs, extended debug info and extremly low performance.

**Custom**: No flags are automatically added.

`PYTHON_EXECUTABLE`

Path to the python executable that Kratos will use. We recommend that you manually set this in case you have multiple versions of python in the system.
Ubuntu users need to be extra careful with this as default versions tends to be Python2, while Kratos is preferably compiled with Python3

`BOOST_ROOT`

Path to boost root directory.

## Adding Applications

In order to add an application you can use the provided macro (`add_app [PATH]` for Linux, `CALL :add_app [PATH]` for Win) along with the route folder of the application that you want to compile. Several examples are provided in the configuration files.

Its now also possible to compile applications outside kratos source dir:

Linux:
```shell
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication    
add_app ${KRATOS_APP_DIR}/FluidDynamicApplication
add_app /home/username/development/ExternalApplication  # Example of external Application
```

Windows:
```shell
CALL :add_app %KRATOS_APP_DIR%/ExternalSolversApplication    
CALL :add_app %KRATOS_APP_DIR%/FluidDynamicApplication
CALL :add_app C:/users/username/development/ExternalApplication  # Example of external Application
```

## Adding Kratos to Path

As Kratos is not an executable but a set of modules and libraries, you will need to add them to the path. In order to do that please add the Kratos install folder (If you didn't touch anything should be `$KRATOS_SOURCE/bin/Release`)

```bash
export PYTHONPATH=$PYTHONPATH:$HOME/Kratos/bin/Release
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Kratos/bin/Release/libs
```

If you are in windows instead do:

```cmd
set PYTHONPATH=%PYTHONPATH%;C:/Kratos/bin/Release
set PATH=%PATH%;C:/Kratos/bin/Release/libs
```

## Examples

This examples are also located in 

### Linux

```bash
# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=gcc
export CXX=g++

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE="Release"
export PYTHON_EXECUTABLE="/usr/bin/python3"

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" -DUSE_MPI=OFF

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j4

```

### Windows

```cmd
rem Set compiler
set CC=cl.exe
set CXX=cl.exe

rem Set variables
set KRATOS_SOURCE=~0,-1%/..
set KRATOS_BUILD=%KRATOS_SOURCE%/build
set KRATOS_APP_DIR=applications

rem Set basic configuration
set KRATOS_BUILD_TYPE=Release
set BOOST_ROOT=C:\boost_1_67_0
set PYTHON_EXECUTABLE=C:\Python37\python.exe

rem Set applications to compile
set KRATOS_APPLICATIONS=
CALL :add_app %KRATOS_APP_DIR%\ExternalSolversApplication;
CALL :add_app %KRATOS_APP_DIR%\StructuralMechanicsApplication;
CALL :add_app %KRATOS_APP_DIR%\FluidDynamicsApplication;

rem Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

rem Configure
@echo on
cmake -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"          ^
-DINCLUDE_FEAST=OFF                                                                                 ^
-DLAPACK_LIBRARIES="C:\CompiledLibs\blas_x64\liblapack.lib"                                         ^
-DBLAS_LIBRARIES="C:\CompiledLibs\blas_x64\libblas.lib"

rem Build
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64
goto:eof

rem Function to add apps
:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof

```

### MacOS

```bash
# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE="Release"
export BOOST_ROOT="/path/to/boost"
export PYTHON_EXECUTABLE="/Library/Frameworks/Python.framework/Versions/3.7/bin/python3"

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
/Applications/CMake.app/Contents/bin/cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
 -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 -L/usr/local/opt/llvm/lib" \
 -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 -L/usr/local/opt/llvm/lib" \
 -DLAPACK_LIBRARIES="/usr/lib/liblapack.dylib" \
 -DBLAS_LIBRARIES="/usr/lib/libblas.dylib" \

# Buid
/Applications/CMake.app/Contents/bin/cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j3

```

## Advanced Configuration

### Building Environment

It is possible to configure the build environment for Kratos, that is: where the source is located, which will be the install dir, and how the python files are going to be installed.

`KRATOS_SOURCE=Path`

Path to the source of Kratos. It will target the directory above this script by default.

`KRATOS_BUILD=Path`

Build directory for Kratos. Makefiles, vsprojects and other artifacts will be stored here. It defaults to Kratos/Build

`KRATOS_APP_DIR=Path`

Path where your applications are located. This variable is not necessary but it helps to organize the list of applications to be compiled. It defaults to Kratos/Applications

`KRATOS_INSTALL_PYTHON_USING_LINKS=ON/OFF`

Controls wether the python files are installed by making copies or creating symlinks to the files in the source directory. This options is specially usefull if you are developing python files and don't want to reinstall every time you touch a script.


### Common Flags

It is also possible to use more advanced configuration options. To use any of these options please add them directly to the cmake configuration line just as any other flag


`-DCMAKE_C_COMPILER=String`

Path to the C compiler. Overrides `CC` environment variable

`-DCMAKE_CXX_COMPILER=String`

Path to the C++ compiler. Overrides `CXX` environment variable

`-DCMAKE_INSTALL_PREFIX=String`

Install path for Kratos. If not set the installation will be done in `bin/[build_type]`

`-DCMAKE_C_FLAGS=String`

User defined flags for the C compiler.

`-DCMAKE_CXX_FLAGS=String`

User defined flags for the C++ compiler. 

`-DBOOST_ROOT=String`

Root directory for boost. Overrided by BOOST_ROOT environmental variable if defined.

`-DPYTHON_EXECUTABLE=String`

Python executable to be used. It is recommended to set this option if more than one version of python is present in the system (For example while using ubuntu). Overrided by PYTHON_EXECUTABLE environmental variable if defined.

`-DINSTALL_RUNKRATOS=ON/OFF`

Enables(Default) or Disables the compilation of the embedded python interpreter (aka Runkratos).

`-DKRATOS_BUILD_TESTING=ON/OFF`

Enables(Default) or Disables the compilation of the C++ unitary tests for Kratos and Applications.

### Compilation Performance
`-DUSE_COTIRE=ON/OFF`

Enables or Disables(default) the use of [cotire](https://github.com/sakra/cotire) to speedup compilation by using unitary builds.
Please notice that enabling this options can greatly increase the amount of memory needed to compile some targets, specially if combined with -jx.

In order to install and compile with this switch please use:
```shell
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target all_unity -- -j1
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install/fast -- -j1 
```

Instead of the regular install target.

### Parallelism
`-DUSE_MPI=ON/OFF`

Enables or Disables(default) the modules and code for mpi. This option is needed if you want to compile Trilinos, Metis, etc...

### Logging
`-DKRATOS_COLORED_OUTPUT=ON/OFF`

Enables colored output of the Logger. If switched on, e.g. warning level messages will be printed in yellow to the terminal. Please notice that colored output is not supported by all terminals.

### External libraries
#### Feast
`-DINCLUDE_FEAST`

Enables or Disables(default) the use of FEAST

#### Metis

`-DUSE_METIS_5=ON/OFF`

Specifies if the metis version is 5 or greater (OFF by default).

`-DMETIS_ROOT_DIR=String`

Root directory for Metis library

#### Trilinos

`-DTRILINOS_ROOT=String`

Root directory for Trilinos library.

`-DTRILINOS_INCLUDE_DIR=String`

Not required if `TRILINOS_ROOT` is set. Path to trilinos include dir.

`-DTRILINOS_LIBRARY_DIR=String`

Not required if `TRILINOS_ROOT` is set. Path to trilinos library dir.

`-DTRILINOS_PREFIX=String`
Indicates the prefix of the trilinos libraries in case they have:
```
libepetra.so          -> No prefix
libtrilinos_epetra.so -> -DTRILINOS_PREFIX="trilinos_"
```

## Applications

Specific compilation information about applications can be found in their own directories:

- [EigenSolvers Application](applications/EigenSolversApplication/README.md#build-instructions)
- [HDF5 Application](applications/HDF5Application/README.md#build-instructions)
- [MultilevelMontecarlo Application](applications/MultilevelMonteCarloApplication/README.md#external-libraries)
- [Poromechanics Application](applications/PoromechanicsApplication/README.md#how-to-use-mpi-in-poromechanics-application)
- [Trilinos Application (Aditional notes)](applications/TrilinosApplication/README.md#notes-for-compilation)
