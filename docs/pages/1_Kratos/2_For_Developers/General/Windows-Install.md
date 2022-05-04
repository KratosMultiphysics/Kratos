---
title: Windows Install
keywords: 
tags: [Windows-Install.md]
sidebar: kratos_for_developers
summary: 
---

# How to compile Kratos:

You can find updated instructions to compile kratos in the [Install.md](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md)

<!--
If you already have **VisualStudio** 2010 or 2012 in your system, there are a couple of automatic compiler/installer both for **32** and **64** bits that will set up all the libraries required by *Kratos* and perform the compilation process automatically. The process is fully automated and you will only have to specify the install directory and click "next" You can found the installers here:

- [For 32 Bits](https://web.cimne.upc.edu/users/croig/data/kratos-3.3.dev-win-32.exe)

- [For 64 Bits](https://web.cimne.upc.edu/users/croig/data/kratos-4.0.dev-win-64.exe)

After using these installers you will have the last version available in the repository compiled.

Notice that these installers will overwrite any previous version of the prerequisite libraries on your machine. The complete list of the libraries installed is the following:

- Python 3.3
- SVN 1.8.0.1
- CMake 3.0.2 
- ACML 4.4

# How to compile Kratos: Windows

In this section we are going to go through the process of compiling a basic version of *Kratos Multiphysics* under *Windows* environments. We recommend you to use *Windows 10* but you can compile *Kratos* with *Windows 7* or higher. 

A basic knowledge of *Windows* is assumed ( execute commands in cmd, create directories, etc...)

- **Tested and working configurations**:
	- *Visual Studio 2015* (Update3) Community with Python 3.4 and Boost 1.67.<br>
	- *Visual Studio 2017* Community with python 3.6 and Boost 1.67.

**We strongly recommend you to use a 64 bit system as some files may require large amounts of ram to compile.**

'''It is highly recommended to install Kratos for 64 bit systems. Thus, all dependent components (such as Python) should be installed for 64 bit systems.'''

## Visual Studio

*Visual Studio* is the only compiler officially supported to build *Kratos* under *Windows*.
Since the adoption of **C++11** we support versions 2015 update 3 onwards. We recommend you to use the latest version of visual studio which can be obtained here for free:

* [Download Visual Studio](https://visualstudio.microsoft.com/en/thank-you-downloading-visual-studio/?sku=Community&rel=15)

Since *Visual Studio* is a multi-language IDE, some distributions come without C++ compiler. Please, make sure that you can create a C++ project before continue, in case C++ packages were missing you will be prompt to download them.

## Git

- Objectives: 
	- Install git
	- Get Kratos Multiphysics source code


The first thing you will need is the *Kratos Multiphysics* source code. To download the code you will have to use a git manager. You can install a git manager from the link below. There are may other git clients that you can use. 

[GitKraken](https://www.gitkraken.com/download)

Once a git client is installed you can clone the code from this url:

~~~
https://github.com/KratosMultiphysics/Kratos
~~~

Once this is done, you should have a "Kratos" directory containing Kratos soruces

## CMake 

- Objectives:
	- Install CMake

*CMake* is the tool used to compile *Kratos*. You can obtain it from its official webpage.

[CMake](http://cmake.org/download/)

Once installing, please <span style="color:red"> do not forget to mark the option: '''"Add CMake to the system PATH for all users"'''</span> 

Please notice that if you want to use *python* 3.4 or higher, you will need *CMake* 3.0.2 or higher.

## Python 

- Objectives:
	- Install Python3

You will need any version of python in your computer in order to compile *Kratos*. We strongly recommend *Python* 3, at least 3.3.4 or higher. you can download python from its official webpage:

[Python](http://www.python.org/downloads/)

Please, take special care to download a installer that suits your desired architecture <span style="color:red">x86 for 32 bits</span>  compilations and <span style="color:red">x86_64 for 64 bits</span>  compilations. Otherwise it won't work.

## BLAS and LAPACK

- Objectives:
	- Get LIBBLAS and LIBLAPACK 

*Blas* and *Lapack* are needed for many solvers, specially those present in the *ExternalSolversApplication*, that you will likely need to compile. You can get these libraries from:

- [For 32 bits](http://icl.cs.utk.edu/lapack-for-windows/lapack/)

Under the section: '''"Prebuilt libraries for Microsoft Visual Studio Projects"'''. Please download both dll and lib files for your architecture.

- [For 64 bits]( http://web.cimne.upc.edu/users/maceli/data/libs.7z)

Please notice that temporally we recommend to use an older version of the libs for 64 bits.

Additionally, you will need some extra dependencies for these libs. The easiest way to fulfil them is to install a proper version of *MinGW* in your system (32 or 64). Any distribution should work, you can find one here:

[MinGW]( http://sourceforge.net/projects/mingw-w64/)

**Warning:** After launching the installer, several options must be selected. Choose **version 6.4.0**. Please take special care to select the correct architecture during this installation (32 bits is called `i868`, 64 bits is called `x86_64`).

## Boost 

- Objectives:
	- Download boost libraries

The next step will consist in obtain Boost. *Kratos Multiphysics* needs *Boost* libraries to support some of its functions. You can use any version from `version 1.67` onward.

[Boost](http://www.boost.org/users/download/)

Extract boost, and note the path as it will be needed in the configure stage to set the `-DBOOST_ROOT` variable.

## Compiling Kratos

### Customize configure.bat

- Objectives:
	- Prepare *Kratos/ configuration file

In the Kratos root folder (**C:\kratos\cmake_build**) copy the  `example_configure.bat.do_not_touch`  to `configure.bat`.
This file controls where Kratos is going to search for the libraries, which applications are going to be installed and how the visual studio solution is going to be generated, among other things.

#### Set the Generator

The first thing you need to do is to tell *CMake* that you intend to build a *VisualStudio* project. This is done automatically by *CMake*, but is highly recommended to add it yourself. To do it add `-G` option followed by your target. For example, if you are using *VisualStudio 2015*:

**For 32 bits:**
```bash
  cmake -G "Visual Studio 14 2015" ^
```

**For 64 bits:**
```bash
  cmake -G "Visual Studio 14 2015 Win64" ^
```
**For 64 bits and VS 2017:**
```bash
  cmake -G "Visual Studio 15 2017 Win64" ^
```

You can find more info and a list of available generators here:

[Link](https://cmake.org/cmake/help/v3.5/manual/cmake-generators.7.html)

**Warning:** If you already configured for 32 bits and compiled, **remove** all the files in the folder `cmake_build` except the .bat file used for configuration. Not removing them can lead to strange errors during compilation.

#### Set Libraries

Once the generator is correctly set, you have to make sure that the paths to all libraries are correctly set. Please make sure that you configure.bat
file has the following lines with the correct path:

```bash
  -DBOOST_ROOT="example/boost_1_67_0/" ^
  -DBLAS_LIBRARIES="example/libblas.lib" ^
  -DLAPACK_LIBRARIES="example/liblapack.lib" ^
```

It is possible that if you have multiple python versions in your system CMake detects the wrong one ( typically, the one with the highest version ) to avoid that, please set it manually. For instance, for Python 3.6:

```bash
  -DPYTHON_EXECUTABLE=C:\Python36\python.exe"
```

#### Enable/Disable Setting

Now, you can enable and  disable the applications you may want to compile or not. For instance:

```bash
  -DSTRUCTURAL_APPLICATION=ON/OFF ^
```

We recommend you to enable:

```bash
  -DINSTALL_EMBEDDED_PYTHON=ON ^
```

#### Example

Here we present a full example using these assumptions:

- You want a x64 build
- You use VisualStudio 2015
- You have boost in "C:\boost_1_67_0"
- You have blas and lapack in "C:\external_libraries"
- You have Python36 in "C:\Python36"
- You downloaded Kratos in "C:\Kratos"
 
```bash
  del CMakeCache.txt
  
  cmake -G "Visual Studio 14 2015 Win64"                                              ^
  -DCMAKE_BUILD_TYPE=Release                                                          ^
  -DCMAKE_CXX_FLAGS=" -D_SCL_SECURE_NO_WARNINGS "                                     ^
  -DBOOST_ROOT="C:\boost_1_67_0"                                                      ^
  -DPYTHON_EXECUTABLE="C:\Python36\python.exe"                                        ^
  -DLAPACK_LIBRARIES="C:\external_libraries\liblapack.lib"                            ^
  -DBLAS_LIBRARIES="C:\external_libraries\libblas.lib"                                ^
  -DMESHING_APPLICATION=ON                                                            ^
  -DEXTERNAL_SOLVERS_APPLICATION=ON                                                   ^
  -DPFEM_APPLICATION=ON                                                               ^
  -DSTRUCTURAL_APPLICATION=ON                                                         ^
  -DCONVECTION_DIFFUSION_APPLICATION=ON                                               ^
  -DFLUID_DYNAMICS_APPLICATION=ON                                                     ^
  -DALE_APPLICATION=ON                                                                ^
  -DFSI_APPLICATION=ON                                                                ^
  -DDEM_APPLICATION=OFF                                                               ^
  -DSWIMMING_DEM_APPLICATION=OFF                                                      ^
  -DINSTALL_PYTHON_FILES=ON                                                           ^
  -DINSTALL_EMBEDDED_PYTHON=ON                                                        ^
  -DEXCLUDE_ITSOL=ON                                                                  ^
  -DSOLID_MECHANICS_APPLICATION=ON                                                    ^
  -DCONSTITUTIVE_MODELS_APPLICATION=ON                                                ^
  ..
```

**Warning:** All these options must be written in the same line, the symbol `^` tells the cmake to read the following line as if it was the same. If you add new options, like the ones in red, do not forget to put this symbol at the end of every line and do not write spaces after this symbol.

### Running configure.bat
- Objectives:
	- Configure Kratos

Once the modifications of the file are done, it needs to be executed. You can do that by executing the command:

```bash
  configure
```

In a `cmd`. Please check the output to ensure that all the paths and libraries are the correct ones. If the configuration has been successful you will se a `Configuration Done` near the end.

### Compiling Kratos
- Objectives:
	- Compile and Install Kratos

Once the configuration script has finished without errors a `KratosMultiphysics.sln` will be generated in the `cmake_build` directory. Double click this file and the visual studio project will open.

Please, make sure that the project is set to `Release` and  `Win32` or `x64` and change it if is not. Finally, in order to compile Kratos, right click the `INSTALL` project and select the option `BUILD`

## Post Compilation

- Objectives:
	- Finish last details
	- Copying necessary files

Once the compilation process has finished, you will need to add a couple of directories to your system path:
- The folder with the blas and lapack libraries (For example: C:\external_libraries)
- The bin folder of the MinGW installation (For example C:\Program Files\ming-gw\x86_64-6.4.0\mingw64\bin)

You may also need to copy the python3X.dll library from your python directory
-->

# Test and Usage

- Objectives:
	- Test that Kratos Works
	- Get familiar with the execution process

Once Kratos is compiled and correctly install all its left is to execute some cases. In this section we will descrive how you can make a simple example to test that kratos works, and the general way to run problems: directly from the command line or using GiD.

## Test the compilation

To to test the compilation, you can prepare a simple script (for example `test_kratos.py`) that contains this line:

```bash
  from KratosMultiphysics import *
```

Then you can execute this script as explained in the next section.

## Executing an script from the command line

The most easy way to execute a KratosMultiphysics script from the command line is to prepare a `.bat` file. A `.bat` file is just a series of commands that are executed together and will make the process simpler. This file should contain only two lines. First the path for the kratos executable and libs folders, and the second, the actual command. For instace:

- Your Downloaded Kratos in `C:\Kratos`
- Your script is called `test_kratos.py`

```bash
  set PATH=C:\\Kratos;C:\\Kratos\\libs;%PATH%
  "C:\\Kratos\\runkratos" test_kratos.py
```

We strongly recommend you to run kratos scripts with the `runkratos` binary inside your *Kratos* installation folder, because it gives the correct values to the environment variables.

```bash
  runkratos test.py
```

Or more exactly, you can go to the folder where your case is (input files and main *python* script) and type:

```bash
  path_to_kratos/runkratos test.py
```

You can also run them directly using the python you have installed in your system (provided that the system knows where python is and the environment variables have the correct values assigned, `PYTHONPATH`, `LD_LIBRARY_PATH` and `PATH`).

```bash
  python test.py
  python3 test.py
  python36 test.py
```

If everything is correct you will see this message:

```bash
   |  /           |             
   ' /   __| _` | __|  _ \   __|
   . \  |   (   | |   (   |\__ \ 
  _|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 6.1.XXXXX
```

## Common problems

### Error C1002 / Error C1060: Out of heap space
This problem **can only be solved in 64bit machine**

This error happens because the compiler is not able to index more than **4GB of RAM**. This is due to the fact that, by default visual studio uses a 32 bit toolset regardless of the target platform. *Visual Studio* 2015 and 2017 both have a 64 bit toolset that does not have this problem. In order to activate it:

Add this to the `*.vsproj` file of the project is giving you problems under the configuration you are using ( `Release`, `Debug`, etc...)

```xml
<PreferredToolArchitecture>x64</PreferredToolArchitecture>
```

or by setting:

```bash
set PreferredToolArchitecture=x64
```

before calling `Msbuild.exe` if you are compiling directly from the cmd.

[Source](https://stackoverflow.com/questions/19820718/how-to-make-visual-studio-use-the-native-amd64-toolchain)

### warning C4273: 'round' : inconsistent dll linkage

You have a conflicting declaration of round. This could happen for example if you have multiple versions of Visual Studio in your computer. Please add `-DHAVE_ROUND` to the configure `cxxflags` entry:

```bash
 -DCMAKE_CXX_FLAGS=" -D_SCL_SECURE_NO_WARNINGS -DHAVE_ROUND"
```

