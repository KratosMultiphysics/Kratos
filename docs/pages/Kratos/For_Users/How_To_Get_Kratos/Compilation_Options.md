---
title: Compilation
keywords: 
tags: [Compilation Flags]
sidebar: kratos_for_users
summary: 
---

This section inteds to give you an overview of the compilation options available in Kratos.

## Building Environment

It is possible to configure the build environment for *Kratos*, that is: where the source is located, which will be the install dir, and how the python files are going to be installed.

`BOOST_ROOT=[Path]`
Sets the path to boost root. If left unassigned Kratos will try to detect a usable boost installed in the system.

`PYTHON_EXECUTABLE=[Path]`

Path to the python executable that *Kratos* will use. We recommend that you manually set this in case you have multiple versions of python in the system.
*Ubuntu* users need to be extra careful with this as default versions tends to be Python2, while *Kratos* is compiled with Python3

`KRATOS_BUILD_TYPE=[Release|RelWithDebInfo|Debug|FullDebug]`
Build type of kratos:   
- Release: Maximum optimization (-O2) level with no debug info.
- RelWithDebInfo: Some optimization (-O1) with debug info. Memory may be optimized in when debuging for some specific locations of memory.
- Debug: No optimization (-O0) with debug info.
- FullDebug: No optimization (-O0) with debug info. Additionally `-DNO_DEBUG` is not active and extra assertions from `-DKRATOS_DEBUG` are enabled.

`KRATOS_SOURCE=[Path]`
Path to the source of *Kratos*. It will target the directory above this script by default.

`KRATOS_BUILD=[Path]`
Build directory for *Kratos*. Makefiles, vsprojects and other artifacts will be stored here. It defaults to `Kratos/Build`

`KRATOS_APP_DIR=[Path]`
Path where your applications are located. This variable is not necessary but it helps to organize the list of applications to be compiled. It defaults to `Kratos/Applications`

`KRATOS_INSTALL_PYTHON_USING_LINKS=[ON|OFF]`
Controls wether the python files are installed by making copies or creating symlinks to the files in the source directory. This options is specially usefull if you are developing python files and don't want to reinstall every time you touch a script.

Using this option in *Windows* requires elevated privileges (you must run the script as admin)

## Common Flags

It is also possible to use more advanced configuration options. To use any of these options please add them directly to the cmake configuration line just as any other flag

`-DCMAKE_C_COMPILER=String`
Path to the C compiler. Overrides `CC` environment variable

`-DCMAKE_CXX_COMPILER=String`
Path to the C++ compiler. Overrides `CXX` environment variable

`-DCMAKE_INSTALL_PREFIX=String`
Install path for *Kratos*. If not set the installation will be done in `bin/[build_type]`

`-DCMAKE_C_FLAGS=String`
User defined flags for the C compiler.

`-DCMAKE_CXX_FLAGS=String`
User defined flags for the C++ compiler.

`-DBOOST_ROOT=String`
Root directory for boost. Overrided by `BOOST_ROOT` environmental variable if defined.

`-DPYTHON_EXECUTABLE=String`
Python executable to be used. It is recommended to set this option if more than one version of python is present in the system (For example while using ubuntu). Overrided by `PYTHON_EXECUTABLE` environmental variable if defined.

`-DKRATOS_BUILD_TESTING=ON/OFF`
Enables(Default) or Disables the compilation of the C++ unitary tests for *Kratos* and Applications.

`-DKRATOS_NO_TRY_CATCH=ON/OFF`
Enables or Disables(Default) the prevention of code generation in `KRATOS_TRY` and `KRATOS_CATCH` macros to allow direct debug of the code through gdb without having to break at `__cxa_throw`.

## Unitary Builds
`-DCMAKE_UNITY_BUILD=ON/OFF`
Enables or Disables(default) the use of [cmake unity build](https://cmake.org/cmake/help/latest/prop_tgt/UNITY_BUILD.html) to speedup compilation by using unitary builds.
Please notice that enabling this options can greatly increase the amount of memory needed to compile some targets, specially if combined with `-jx`.

In order to install and compile with this switch please use:

On *GNU/Linux*
```console
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j1
```
On *Windows*
```console
cmake --build "%KRATOS_BUILD%/%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64
```

Instead of the regular install target.

Please, beware that using this flag along with a parallel compilation may cause a **VERY LARGE** use of RAM as we hardcoded *Kratos* compilation so unitary builds try to make as many unitary targets as threads are usable.
We recommed you to disable parallel compilation unless you know what you are doing.

## Parallelism

Building *Kratos* with support for **MPI** requires an advanced configuration of its building script, as well as the building of dependencies and parallel applications.
Here you can find [guidelines](https://github.com/KratosMultiphysics/Kratos/wiki/Compiling-Kratos-with-MPI-support) for the compilation of Kratos and its dependencies.

## Logging
`-DKRATOS_COLORED_OUTPUT=ON/OFF`
Enables colored output of the Logger. If switched on, e.g. warning level messages will be printed in yellow to the terminal. Please notice that colored output is not supported by all terminals.

## TPL-Libraries
*Kratos* can make use of TPL libraries that cannot be included in the main compilation processes for a variety of reasons
If you want to add support for those libraries, we provide switches to enable them.

Please note that **Kratos will NEVER DISTRIBUTE, RELEASE or COMPILE** with these libraries unless explicitly specified, and the use of these libraries may add additional restrictions on top of BSD.

### Tetgen
[Tetgen](http://wias-berlin.de/software/tetgen/) is a library to generate tetrahedral meshes. We provide some utilities that can make use of *Tetgen*. The flags related with *Tetgen* are the following:

`-DUSE_TETGEN_NONFREE_TPL=ON`
Enables/Disables the use of *Tetgen* and its related utilities in the code. If no other options provided *Kratos* will try to find *Tetgen* installed on your system.

`-DUSE_TETGEN_NONFREE_TPL_PATH="${TETGEN_PATH}"`
Tries to use a local version of *Tetgen* from a given `TETGEN_PATH`

`-DUSE_TETGEN_NONFREE_TPL_URL="${TETGEN_URL}"`
Tries to download and use a version of *Tetgen* from a given `TETGEN_URL`

`-DFORCE_TETGEN_NONFREE_TPL_URL`
Forces to re-download and replace an existing version of *Tetgen* obtained through `USE_TETGEN_NONFREE_TPL_URL`

### Triangle
[Triangle](http://www.cs.cmu.edu/~quake/triangle.html) is a library for delaunay triangulation. We provide some utilities in *Kratos* that depend on this library. The flags related with *Triangle* are the following:

`-DUSE_TRIANGLE_NONFREE_TPL`
Enables or disables the use of *Triangle* and its related utilities in the code.