---
title: Windows Install (MinGW)
keywords: 
tags: [Windows-Install-(MinGW).md]
sidebar: kratos_for_developers
summary: 
---

In this section we are going to go through the process of compiling a basic version of Kratos Multiphysics under Windows environments with MinGW. The process is something in between the compilation in Linux and in Windows using Visual Studio.

_Please notice that this guide is introductory and many of the steps and software used here can be customized. The only hard restriction is to use a C++11 compatible compiler._

## MSYS2

First we download MSYS2 in the following [link](https://www.msys2.org/). This will install MinGW, which allows to easiy install packages *a la* Arch-Linux (Pacman package manager). We install it, and with it the first thing we do is to update as follows (in the MSYS2 bash):

```Shell
pacman -Syu
```
It is very relevant to add to the Windows PATH your `msys64\mingw64\bin` folder in order that the system locates the binaries.

## Git

* Objectives: 
 * Install git
 * Get Kratos Multiphysics source code

The first thing you will need is the *Kratos* Multiphysics source code. To download the code you will have to use a git. You can install the default git by using this command:

```Shell
pacman -S git
```

Once git is installed you can fetch the code by using these commands:

```Shell
git clone https://github.com/KratosMultiphysics/Kratos Kratos
```

## Dev Packages 

* Objectives:
 * Get G++
 * Get Fortran compiler
 * Get LIBBLAS and LIBLAPACK 

You will need a series of packages with some *Kratos* dependencies. The command below will install all the packages needed. It will allow you to compile with python2 or python3 so you will be able to chose later ( we recommend python3 )

```Shell
pacman -S mingw64/mingw-w64-x86_64-lapack mingw64/mingw-w64-x86_64-openblas mingw-w64-x86_64-cotire mingw64/mingw-w64-x86_64-cmake mingw64/mingw-w64-x86_64-clang mingw64/mingw-w64-x86_64-gcc mingw64/mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-make
```

## Python 

- Objectives:
	- Install Python3

You will need any version of python in your computer in order to compile *Kratos*. We strongly recommend *Python* 3, at least 3.3.4 or higher. You can download python from its official webpage:

[Python](http://www.python.org/downloads/)

Please, take special care to download a installer that suits your desired architecture <span style="color:red">x86 for 32 bits</span>  compilations and <span style="color:red">x86_64 for 64 bits</span>  compilations. Otherwise it won't work.

Unfortunately we cannot use right now MSYS2 directly, as the development files are not available (python3-dev equivalent to Linux).

## Boost
* Objectives:
 * Download Boost

The next step will consist in downloading Boost. For that we will use pacman again:

```Shell
pacman -S mingw64/mingw-w64-x86_64-boost
```

## Configure

* Objectives:
 * Configure *Kratos* for the first time compilation

In order to compile *Kratos* for the first time you will need to configure the project.  First, navigate to your `Kratos\scripts` folder and make a copy of the template files:

```Shell
copy standard_configure_MINGW.sh ..\build\configure.sh
copy standard_configure_MINGW.bat ..\build\configure.bat
```

We only execute the configure.bat, but we need to edit both. The configure.bat you just need to specify the location of your msys2 folder so it can find the sh.exe (unfortunately sh.exe cannot be in the PATH, so the two files are required). 

Then, open configure.sh with any text editor and modify the lines that tell cmake where some components are located.
You will need to provide at least '''BOOST_ROOT'' and '''PYTHON_EXECUTABLE'''.

*Kratos* has moved to C++11 recently, Please mind to add the "-std=c++11" to your compiler of choice. If you follow the example below, it is already present ( notice the flag in CMAKE_CXX_FLAGS, highlighted in bold)

For example, it will look something like:

For the `bat`:
```
:: this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
:: to use it, copy it to another file and run it

cls

rem Set variables
set KRATOS_SOURCE=~0,-1%/..
set KRATOS_BUILD=%KRATOS_SOURCE%/build

rem Set basic configuration
if not defined KRATOS_BUILD_TYPE set KRATOS_BUILD_TYPE=Release

:: you may want to decomment this the first time you compile
rem Clean
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

location_msys64\msys64\usr\bin\sh.exe configure.sh
```

For the `sh`:
```
# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=gcc
export CXX=g++

# Set variables
export KRATOS_SOURCE="kratos_folder/Kratos"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-"Release"}
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"python_folder/python.exe"}

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication
add_app ${KRATOS_APP_DIR}/MeshingApplication
add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication

# Configure
cmake ..                                                                                            \
-G "MinGW Makefiles"                                                                                \
-DWIN32=TRUE                                                                                        \
-DCMAKE_EXE_LINKER_FLAGS="-s"                                                                       \
-DCMAKE_SHARED_LINKER_FLAGS="-s"                                                                    \
-DLAPACK_LIBRARIES="msys_folder/msys64/mingw64/lib/liblapack.dll.a"                                 \
-DBLAS_LIBRARIES="msys_folder/msys64/mingw64/lib/libblas.dll.a"                                     \
-H"${KRATOS_SOURCE}"                                                                                \
-B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"                                                            \
-DUSE_MPI=OFF                                                                                       \

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j4
    
// More options ( do not include this line )
```

Note that the ```\``` must be the last character in the line. Even an space after it will cause an error! (and the returned message is completely misleading, so be careful with this!!)

Notice that you can also turn on/off parts of the code according to your necessities: [Compilation options](Compilation-options)

:warning: Cmake requires all definitions in a single line! The line concatenation character ```\``` therefore MUST NOT be followed by any whitespace in the same line as this would prevent the cmake from running the lines below

## Compile 

* Objectives:
 * Compile *Kratos*.

If you followed all the steps correctly, compiling *Kratos* should be as easy as executing the configure script:

```Shell
configure.bat
```
Please, notice that *Kratos* is big and the compilation process can easily take 1 or 2 hours, depending on which applications are being compiled. A typical compilation process with the default configuration takes approximately 45 minutes with a i7 / 8GB Ram computer.

## Setting up your enviroment

* Objectives:
 * Tell Windows how to execute *Kratos*

Once *Kratos* is compiled, you will have to tell the OS where to find the libraries.  You need to add to PYTHONPATH (environment variables) the *Kratos* folder, and to PATH (environment variables) the libs folder (usually the kratos_folder/libs).

If you have enabled the embedded python option -DINSTALL_EMBEDDED_PYTHON=ON, you can also add to PATH (environment variables) the *Kratos* folder in order to have the "runkratos" available as a regular command.

## Test

* Objectives:
 * Tests kratos

To to tests the compilation, you can execute a simple python script containing this line:

```python
from KratosMultiphysics import *
```

We strongly recommend you to run kratos scripts with the "runkratos" binary that will be generated inside your Kratos installation folder. You can also run them by using python (if you have compiled with python version 2.x.x), or python3 (if you have compiled with python version 3.x.x)

* runkratos test.py
* python test.py
* python3 test.py

If everething was ok you will see this message:

```
     |  /           |             
     ' /   __| _` | __|  _ \   __|
     . \  |   (   | |   (   |\__ \ 
    _|\_\_|  \__,_|\__|\___/ ____/
               Multi-Physics 7.0.11016
```