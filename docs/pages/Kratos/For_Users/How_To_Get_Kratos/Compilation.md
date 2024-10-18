---
title: Linux Build
keywords: 
tags: [Linux-Build.md]
sidebar: kratos_for_users
summary: 
---

## Cloning Kratos

In order to obtain the source code of *Kratos*, you will need to clone the repository using git.

You can install git through the following command in *GNU/Linux*:
```consoleell
sudo apt-get install git
```
In *Windows*, you can download it in:

* [Download Git](https://git-scm.com/downloads)

Once git is installed you can fetch the code by using this command in a terminal:

```consoleell
git clone https://github.com/KratosMultiphysics/Kratos Kratos
```

## Kratos Dependencies

### Kratos Core Dependencies
  These are the basic dependecies needed to compile the *Kratos* Core and most of the *Kratos* applications.
  * Python3-dev
  * C++17 compiler
  * CMake
  * Boost (dependencies are header-only, no compilation of boost libraries required)

Additionaly, Visual Studio is required to compile in *Windows*.

- #### Linux/WSL installation

    The command below will install all the packages needed.

    ```consoleell
    sudo apt-get install python3-dev gcc g++ cmake libboost-all-dev
    ```
    Newer versions of boost can be downloaded in:

    http://www.boost.org/users/download/.

- #### Windows installation

    - Visual Studio

        *Visual Studio* is the only compiler officially supported to build *Kratos* under *Windows*. The minimium required version is *Visual Studio 2019* or higher (please ensure that version is at least 16.8 or MSVC version is at least 19.24).

        * [Download Visual Studio](https://visualstudio.microsoft.com/en/thank-you-downloading-visual-studio/?sku=Community&rel=16)

        Since *Visual Studio* is a multi-language IDE, some distributions come without C++ compiler. Please, make sure that you can create a C++ project before continuing, in case C++ packages were missing you will be prompt to download them. You can install the **Desktop development with C++** workload with the Visual Studio Installer to acquire all necessary depencencies to compile C++ projects.

        When compiling *Kratos* in *Windows*, please take into consideration the [Windows Visual Studio compilation configuration](#Windows-Visual-Studio-compilation-configuration).

    - CMake
        * [Download CMake](http://cmake.org/download/)

        Once installing, please **do not forget to mark the option: '''"Add CMake to the system PATH for all users"'''**

        Minimum required version: CMake 3.20

    - Python

        You will need at least *Python* 3.8 (recommended 3.8/3.9/3.10) in your computer in order to compile *Kratos*. You can download python from its official webpage:

        * [Download Python](http://www.python.org/downloads/)

        Please, take special care to download a installer that suits your desired architecture **x86 for 32 bits**  compilations and **x86_64 for 64 bits**  compilations. Otherwise it won't work.

    - Boost

        The next step will consist in obtain Boost. *Kratos Multiphysics* needs *Boost* libraries to support some of its functions. You can use any version from `version 1.67` onward.

        * [Download Boost](http://www.boost.org/users/download/)

        Extract boost, and note the path as it will be needed in the configure stage to set the environmental variable `BOOST_ROOT`.

- #### MinGW
  MingGw compilation details are hidden by default to avoid confusion, please click the button below to show them.
  <details>
    <summary>Show MinGW compilation details</summary>
  
    *MinGW* means minimal GNU for *Windows*. There are different manners of installing, the simplest one using *MSYS2*.
  
    - MSYS2
  
        First, we download *MSYS2* in the following [link](https://www.msys2.org/). This will install *MinGW*, which allows to easiy install packages *a la* Arch-Linux (Pacman package manager). We install it, and with it the first thing we do is to update as follows ([in the *MSYS2* bash](https://www.msys2.org/docs/terminals/)):
        ![](https://www.msys2.org/docs/mintty.png) ![](https://www.msys2.org/docs/launchers.png)
  
        ```consoleell
        pacman -Syu
        ```
  
        It is very relevant to add to the *Windows* `PATH` your `msys64\mingw64\bin` folder in order that the system locates the binaries.
  
    - Git
  
        The first thing you will need is the *Kratos* Multiphysics source code. To download the code you will have to use a git. You can install the default git by using this command:
  
        ```consoleell
        pacman -S git
        ```
  
        Once git is installed you can fetch the code by using these commands:
  
        ```consoleell
        git clone https://github.com/KratosMultiphysics/Kratos Kratos
        ```
    - Dev Packages
  
        You will need a series of packages with some *Kratos* dependencies. These include the compilers (*GCC*,*Clang/LLVM*), *CMake*, *Blas and Lapack* libraries and the *OpenMP* support. The command below will install all the packages needed. The command below will install all the packages needed.
  
        ```consoleell
        pacman -S mingw64/mingw-w64-x86_64-lapack mingw64/mingw-w64-x86_64-openblas mingw64/mingw-w64-x86_64-cmake mingw64/mingw-w64-x86_64-clang mingw64/mingw-w64-x86_64-gcc mingw64/mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-make mingw64/mingw-w64-x86_64-openmp mingw64/mingw-w64-x86_64-dlfcn
        ```
  
    - Python 
        You will need at least *Python* 3.8 (recommended 3.8/3.9/3.10) in your computer in order to compile *Kratos*. You can download python from its official webpage:
  
        * [Download Python](http://www.python.org/downloads/)
  
        Please, take special care to download a installer that suits your desired architecture **x86 for 32 bits**  compilations and **x86_64 for 64 bits**  compilations. Otherwise it won't work.
  
        Unfortunately, we cannot use right now *MSYS2* directly, as the development files are not available (`python3-dev` equivalent to *GNU/Linux*).
  
    - Boost
        The next step will consist in obtain Boost. *Kratos Multiphysics* needs *Boost* libraries to support some of its functions. You can use any version from `version 1.67` onward. For that, we will use `pacman` again:
  
        ```consoleell
        pacman -S mingw64/mingw-w64-x86_64-boost
        ```
  
    ##### Using UCRT64
  
    UCRT (Universal C Runtime) is a newer version which is also used by Microsoft Visual Studio by default, see https://www.msys2.org/docs/environments/. It should work and behave as if the code was compiled with MSVC.
  
    - Better compatibility with MSVC, both at build time and at run time.
    - It only ships by default on Windows 10 and for older versions you have to provide it yourself or depend on the user having it installed.
  
    If using UCRT64 the dependencies will be like:
  
    ```consoleell
    pacman -S ucrt64/mingw-w64-ucrt-x86_64-lapack ucrt64/mingw-w64-ucrt-x86_64-openblas ucrt64/mingw-w64-ucrt-x86_64-cmake ucrt64/mingw-w64-ucrt-x86_64-clang ucrt64/mingw-w64-ucrt-x86_64-gcc ucrt64/mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-make ucrt64/mingw-w64-ucrt-x86_64-openmp ucrt64/mingw-w64-ucrt-x86_64-dlfcn ucrt64/mingw-w64-ucrt-x86_64-boost
    ```
  </details>
  
### Specific Application Dependencies

Some applications have additional dependencies. Please check the `README` files of the applications that are compiled

## Basic Configuration

You can find the new kratos configuration file in *Kratos* `scripts` folder: `standard_configure.sh` for *GNU/Linux*, `standard_configure_mac.sh` for *MacOS*, `standard_configure.bat` for *Windows* and others. In the special case of *Windows* using *MinGW* you will need to copy two scripts (`standard_configure_MINGW.bat` and `standard_configure_MINGW.sh`) both are required, but only the `.bat` file is invoked.

Out of the box *Kratos* will try to find all necessary libraries in your system automatically, but we recommend you to copy these scripts and modify it according to your preferences. Please take a look at the following configuration options:

`KRATOS_BUILD_TYPE`

Compilation Type. Options are `Release`,`RelWithDebInfo`,`Debug`,`FullDebug`,`Custom`

**Release**: Full Release with maximum optimization options and no debug Info.

**RelWithDebInfo**: Full Release with optimization and debug info. Adecuate to debug simple problems without losing performance.

**Debug**: Debug build with no optimization flags.

**FullDebug**: Debug build with no optimization flags, extended debug info and extremly low performance.

**Custom**: No flags are automatically added.

`PYTHON_EXECUTABLE`

Path to the python executable that *Kratos* will use. We recommend that you manually set this in case you have multiple versions of python in the system.
*Ubuntu* users need to be extra careful with this as default versions tends to be Python2, while *Kratos* is compiled with Python3

`BOOST_ROOT`

Don't use this unless you have problems during the compilation. Path to boost root directory, set it if you downloaded but without using `apt-get`.

## Configuration scripts examples

These examples are also located [in the /scripts folder](https://github.com/KratosMultiphysics/Kratos/tree/master/scripts). You can simply create your own copy:

```consoleell
cp /path_to_kratos/scripts/standard_configure.sh /path_to_kratos/scripts/configure.sh
```
Then, these scripts can be launched through the system terminal.

*GNU/Linux*

```consoleell
sh /path_to_kratos/scripts/configure.sh
```

**NOTE**: In case the compiler runs out of memory, try increasing the swap size to at least 16 GB and re-starting the compilation process.

*Windows*

```consoleell
./path_to_kratos/scripts/configure.bat
```

**NOTE**: In case os compiling with *Visual Studio*, after installing *Visual Studio*, in some *Windows* systems the console does not have direct access to the *Visual Studio Compiler*. In order to make sure the compiler is available, try typing `cl`. Use this console to compile *Kratos* if the compiler responds. In case of error, instead of using the standard *Windows* console, open the *Native Tools Command Prompt* console and launch the compilation from there.

The example scripts for every system are shown next.

### GNU/Linux

```console
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

# Set basic configuration
export KRATOS_BUILD_TYPE="Release"
export PYTHON_EXECUTABLE="/usr/bin/python3"

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" -DUSE_MPI=OFF -DUSE_EIGEN_MKL=OFF

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j4
```

### Windows

#### Visual Studio

```cmd
rem Set compiler
set CC=cl.exe
set CXX=cl.exe

rem Set variables
if not defined KRATOS_SOURCE set KRATOS_SOURCE=%~dp0..
if not defined KRATOS_BUILD set KRATOS_BUILD=%KRATOS_SOURCE%/build

rem Warning: In windows this option only works if you run through a terminal with admin privileges
rem set KRATOS_INSTALL_PYTHON_USING_LINKS=ON

rem Set basic configuration
if not defined KRATOS_BUILD_TYPE set KRATOS_BUILD_TYPE=Release
if not defined BOOST_ROOT set BOOST_ROOT=C:\CompiledLibs\boost_1_67_0
if not defined PYTHON_EXECUTABLE set PYTHON_EXECUTABLE=C:\Windows\py.exe

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
cmake -G"Visual Studio 16 2019" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"  ^
-DUSE_EIGEN_MKL=OFF        ^
-DCMAKE_CXX_FLAGS=" %KRATOS_PARALLEL_BUILD_FLAG% "

rem Build
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64
goto:eof

rem Function to add apps
:add_app
set KRATOS_APPLICATIONS=%KRATOS_APPLICATIONS%%1;
goto:eof
```

##### Windows Visual Studio compilation configuration

Some of the parameters detailed in the example script above may vary from system to system, or the *Visual Studio* version.

If you are using *Visual Studio 2017*, the configure command should be:
```
cmake -G"Visual Studio 15 2017" -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"  ^
-DUSE_EIGEN_MKL=OFF
```
You can check the specific *Visual Studio* version that you have installed in your system, by checking the *Visual Studio* Tab 'Help' > 'About Microsoft Visual Studio'.

If you have a 64-bit system, you might need to also specify it in the configure command for some *Visual Studio* versions with the flag ```-A x64```.

```
cmake -G"Visual Studio 15 2017" -A x64 -H"%KRATOS_SOURCE%" -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"  ^
-DUSE_EIGEN_MKL=OFF
```
#### MinGW

In the case of *MinGW* two scripts are required, one is the *Command Prompt* for *Windows*:

```cmd
cls

@REM Set variables
if not defined KRATOS_SOURCE set KRATOS_SOURCE=%~dp0..
if not defined KRATOS_BUILD set KRATOS_BUILD=%KRATOS_SOURCE%/build

@REM Set basic configuration
if not defined KRATOS_BUILD_TYPE set KRATOS_BUILD_TYPE=Release
@REM Decomment the following in case of considering MKL
@REM if not defined MKLROOT set MKLROOT=C:\PROGRA~2\Intel\oneAPI\mkl\latest\

@REM rem setting environment variables for using intel MKL
@REM call "%MKLROOT%\env\vars.bat" intel64 lp64

:: you may want to decomment this the first time you compile
@REM Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

sh %KRATOS_BUILD%\configure_MINGW.sh
```

And the second is the bash script that will be called by the former script (it is similar to the one in *GNU/Linux*):

```console
# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler #NOTE: Currently only GCC is supported, linking error on recent versions of Clang/LLVM, see https://github.com/llvm/llvm-project/issues/53433
export CC=${CC:-gcc}
export CXX=${CXX:-g++}

# Set variables
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON
export KRATOS_SHARED_MEMORY_PARALLELIZATION=${KRATOS_SHARED_MEMORY_PARALLELIZATION:-"OpenMP"}

# Set basic configurati
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"C:/Windows/py.exe"}

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication

# Configure
cmake ..                                                                                            \
-G "MinGW Makefiles"                                                                                \
-DWIN32=TRUE                                                                                        \
-DCMAKE_INSTALL_PREFIX="${KRATOS_SOURCE}/bin/${KRATOS_BUILD_TYPE}"                                  \
-DCMAKE_BUILD_TYPE="${KRATOS_BUILD_TYPE}"                                                           \
-DCMAKE_EXE_LINKER_FLAGS="-s"                                                                       \
-DCMAKE_SHARED_LINKER_FLAGS="-s"                                                                    \
-H"${KRATOS_SOURCE}"                                                                                \
-B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"                                                            \
-DUSE_MPI=OFF                                                                                       \
-DKRATOS_SHARED_MEMORY_PARALLELIZATION="${KRATOS_SHARED_MEMORY_PARALLELIZATION}"                    \
-DKRATOS_GENERATE_PYTHON_STUBS=ON                                                                   \
-DUSE_EIGEN_MKL=OFF

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j$(nproc)
```
### MacOS

```console
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
export PYTHON_EXECUTABLE="/Library/Frameworks/Python.framework/Versions/3.8/bin/python3"

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
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
 -DUSE_EIGEN_MKL=OFF

# Buid
/Applications/CMake.app/Contents/bin/cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j3

```

## Adding Applications

In order to add an application you can use the provided macro (`add_app [PATH]` for *GNU/Linux*, `CALL :add_app [PATH]` for Win) along with the route folder of the application that you want to compile. Several examples are provided in the configuration files.

Its now also possible to compile applications outside *Kratos* source dir:

*GNU/Linux*:
```consoleell
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicApplication
add_app /home/username/development/ExternalApplication  # Example of external Application
```

*Windows*:
```consoleell
CALL :add_app %KRATOS_APP_DIR%/LinearSolversApplication
CALL :add_app %KRATOS_APP_DIR%/FluidDynamicApplication
CALL :add_app C:/users/username/development/ExternalApplication  # Example of external Application
```

For *Windows* with *MinGW* it works in the same way as in *GNU/Linux* as it works as a bash script.

## Post Compilation

As *Kratos* is not an executable but a set of modules and libraries, you will need to add them to the path. In order to do that please add the *Kratos* install folder (If you didn't touch anything should be `$KRATOS_SOURCE/bin/Release`)

### GNU/Linux
```console
export PYTHONPATH=$PYTHONPATH:$HOME/Kratos/bin/Release
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Kratos/bin/Release/libs
```
Or set them permanently by adding these lines in your `~/.bashrc`.

### Windows
In a *Command Prompt:*

```cmd
set PYTHONPATH=%PYTHONPATH%;C:/Kratos/bin/Release
set PATH=%PATH%;C:/Kratos/bin/Release/libs
```
In *Windows Powershell*:
```cmd
$Env:PYTHONPATH+=";C:/Kratos/bin/Release"
$Env:PATH+=";C:/Kratos/bin/Release/libs"
```

Or set them permanently using the  **Edit the system environment variables** option in the Control panel.

You can then test your compilation by executing an example script or trying to import the python module

```python
from KratosMultiphysics import *
```

The result should be:

```
 |  /           |                  
 ' /   __| _` | __|  _ \   __|    
 . \  |   (   | |   (   |\__ \  
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 9.X.Y-4afb88094a-Release-ARM64
           Compiled for GNU/Linux and Python3.8 with GCC-8.5
Compiled with threading and MPI support.
Maximum number of threads: 1.
Running without MPI.
```
## Advanced Configuration

### Compilation of Kratos in parallel

We provide several flavours in order to parallelize *Kratos* compilation. We have divided this option according to the operating system specifics.

#### GNU/Linux

*GNU/Linux* builds should automatically make use of the maximum number of threads in your computer which is passed to the compiler in the `-j$(nproc)` flag on the last line of the configure file:
```
# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j$(nproc)
```

If your *GNU/Linux* flavour does not support the `$(nproc)` shortcut or you simply want to tune this value to some of your liking, you can change it:
```
# Buid (This will make it compile with 2 threads)
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j2
```
**Warning**: Please be carefull while mixing parallel builds with unitay builds. See [below](#unitary-builds)

#### Windows

**NOTE:** The following will only apply in the case of compiling with *Visual Studio*, in case of compiling with *MinGW* the same tips as in *GNU/Linux* can be applied.

*Windows* should detect automatically the number of threads of your computer, but many times this mechanism fails. We included several options in order to force the parallel compilation:

You can force it manually by commenting this lines in the configuration file, and adding a number of processes of your choice:
```ps1
rem Enable this if your build is slow and you have a multi-core machine
rem set KRATOS_PARALLEL_BUILD_FLAG=/MPX
```

This will pass the `/MPX` option directly to `CL.exe`, where `X` is the number of threads you want to use.

If you preffer to interact directly with `MSBuild.exe` you can use either of this options in the cmake build command:
- `/p:CL_MPcount=X`: Enable multiples cpp to be compiled in parallel
- `/m:x`: Enable multiple applications to be compiled in parallel

Example using 4 threads and a single project
```ps1
rem Build
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64 /p:CL_MPcount=4 /m:1
```

Example using 2 threads and 2 project ( total io 4 threads )
```ps1
rem Build
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64 /p:CL_MPcount=2 /m:2
```

Finally you can set parallelism options in the *VisualStudio IDE*.

**Warning**: Please be careful while mixing parallel builds with unitary builds. See [below](#unitary-builds)

#### MacOS

There is no dedicated support for parallel builds in *MacOS*, but *GNU/Linux* options should behave very similarly. If you detect a problem please inform us and we will try to update this section with the specifics.