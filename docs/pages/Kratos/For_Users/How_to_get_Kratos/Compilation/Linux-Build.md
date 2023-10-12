---
title: Linux Build
keywords: 
tags: [Linux-Build.md]
sidebar: kratos_for_users
summary: 
---

## How to compile Kratos:

You can find latest Kratos compilation instructions in [Install.md](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md)

## Test and Usage

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
               Multi-Physics 3.3.11016
```

## Common problems

In this section we provide will try to provide solutions to the most common problems and issues that may appear during the compilation process

### Cannot find KratosMultiphysics

Make sure that ```LD_LIBRARY_PATH``` and ```PYTHONPATH``` are pointing to your kratos folder.
To make sure that the variables are set correctly, you can always print their value from the terminal by typing:
```Shell
echo $LD_LIBRARY_PATH
```
```Shell
echo $PYTHONPATH
```
```Shell
echo $PATH
```

Some shell interpreters have issues with the separator token. If environment variables above listed are correctly set, try to remove the ":" if kratos is the only path there.

```Shell
# This may cause problem
echo $PYTHONPATH
:path/to/kratos

# Should be
export PYTHONPATH=path/to/kratos
echo $PYTHONPATH
path/to/kratos
```

### I am getting Python link errors

This errors appear if the version of python used to compile boost is not the same as the one
used by Kratos.

There are several causes that may be causing this. Pleas try the following:

#### Python version mismatch
The most provable reason for this error to happend is a missmatch between the python versions used by Kratos and Boost. Please, double check you have the same version of python in the projet-config.jam (boost) and configure.sh (Kratos) files.

#### Old version of CMake
Please check that CMake version is 3.14 or newer. In general newer versions of Python will need newer versions of CMake to work.

It has been observed that compiling with IDE's ( QTCreator, Netbeans, ...) sometimes causes this error as well.
If you are experiencing this problem, try to modify the configure.sh script and replace '''cmake''' by the absolute path of CMake 3.14 or the any newer version:

```
/path/to/correct/cmake ..                                                     \
-DCMAKE_C_COMPILER=/usr/bin/gcc                                               \
...
```

If for some reason you have to use an older version of CMake you can manually add support for an outdated versiof of python by adding the version to these files:

```
/usr/share/cmake-3.x/Modules/FindPythonLibs.cmake
/usr/share/cmake-3.x/Modules/FindPythonInterp.cmake
```

Please add the version at the begining of the list:

```CMake
SET(_PYTHON1_VERSIONS 1.6 1.5)
SET(_PYTHON2_VERSIONS 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0)
SET(_PYTHON3_VERSIONS 3.4 3.3 3.2 3.1 3.0)
```

### I am getting lots of warnings when I compile Kratos

It is known that in some cases warnings appear while including boost files due to the fact that the flag **"-Wunused-local-typedefs"** is set by default in gcc.

This does not have any impact on the final code, but if you want a cleaner output you can add the flag **"-Wno-unused-local-typedefs"** to the configuration files:

```CMake
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -Wno-unused-local-typedefs"      \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 -Wno-unused-local-typedefs"          \
```

## Compiling with MPI
If you want to compile the MPI-version of Kratos, you need to at least to
* Have a MPI installation (e.g. [OpenMPI](http://www.open-mpi.de/) or [Intel MPI](https://software.intel.com/en-us/intel-mpi-library)). There are no known restrictions on which version to use.
* Set the flag `MPI_NEEDED` to `ON` in the configure-script

This will compile the core of Kratos with MPI (i.e. compile the MPI-Core) and provide basic features of MPI. Check [here](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/mpi/tests/test_mpi_data_communicator_python.py) to see which functionalities are exposed to Python.

If you want to run Kratos in MPI together with e.g. the `FluidDynamicsApplication` or the `StructuralMechanicsApplication`, then you also need [Trilinos](https://trilinos.org/) for solving the distributed system of equations and [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for partitioning of the problem.
In order to use them it is required to compile the `TrilinosApplication` and the `MetisApplication` respectively. These applications require to install the corresponding libraries first.

### MetisApplication
#### Installing METIS using the package-manager
On Ubuntu 18.04, the following command installs the necessary files:
```Shell
sudo apt-get install libmetis-dev
```
#### Compiling METIS yourself
Compiling is straight-forward using the build-instructions that come with the download.

#### Compiling the `MetisApplication`
Afterwards add the following lines to your configure-script:
```Shell
-DMETIS_APPLICATION=ON \
-DUSE_METIS_5=ON \
```

### TrilinosApplication

#### Installing Trilinos using the package-manager
On Ubuntu 18.04, the following command installs the necessary files:
```Shell
sudo apt-get install trilinos-all-dev
```
Add the following to your configure-script:
```Shell
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \
-DTRILINOS_LIBRARY_PREFIX="trilinos_" \
```
#### Compiling Trilinos yourself
Some information on how to compile Trilinos (version 12.10.1) is provided in the following, but this might not be applicable for all systems. Please consult the build-instructions of Trilinos itself if needed or open an issue.

Add the following to your configure-script if you compile Trilinos yourself:
```Shell
-DTRILINOS_ROOT="path/to/trilinos" \
```

Depending on the configuration it might also be necessary to specify `TRILINOS_SOURCE` for CMake to find the correct trilinos-libraries.
***

It is recommended to create a build directory in the root folder of Trilinos. Inside this folder, create a `do-configure.sh` to compile Trilinos.

If you use LAPACK 3.6 or newer, then you have to apply the following fixes: (sources: ​[1](https://github.com/gahansen/Albany/wiki/ALCF-Vesta), [​2](https://www.dealii.org/developer/external-libs/trilinos.html)):

In `trilinos-12.10.1-Source/packages/epetra/src/Epetra_LAPACK_wrappers.h`

Change
```
#define DGGSVD_F77  F77_BLAS_MANGLE(dggsvd,DGGSVD)
#define SGGSVD_F77  F77_BLAS_MANGLE(sggsvd,SGGSVD)
```
to
```
#define DGGSVD_F77  F77_BLAS_MANGLE(dggsvd3,DGGSVD3)
#define SGGSVD_F77  F77_BLAS_MANGLE(sggsvd3,SGGSVD3)
```

This is an example for a `do-configure.sh`. Of course, the paths have to be adjusted
```
TRILINOS_ROOT="${HOME}/software/kratos/trilinos"
EXTRA_LINK_FLAGS=""
EXTRA_ARGS=$@
MPI_ROOT="${HOME}/software/ompi"
METIS_ROOT="{HOME}/software/kratos/parmetis/ParMetis-3.2.0"
LAPACK_ROOT="${HOME}/software/kratos/lapack/lapack-3.7.0"

rm CMakeCache.txt
rm -rf CMakeFiles
rm *.cmake

cmake \
  -D CMAKE_INSTALL_PREFIX:FILEPATH="${TRILINOS_ROOT}" \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D Trilinos_SOURCE_DIR="${TRILINOS_SOURCE}" \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_Fortran_COMPILER=mpif90 \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D Trilinos_ENABLE_CXX11=ON\
  -D MPI_BASE_DIR="${MPI_ROOT}" \
  -D MPI_INCLUDE_DIRS:PATH="${MPI_ROOT}/include" \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D BLAS_LIBRARY_DIRS:FILEPATH="${LAPACK_ROOT}" \
  -D BLAS_LIBRARY_NAMES:STRING="librefblas.a" \
  -D LAPACK_LIBRARY_DIRS:FILEPATH="${LAPACK_ROOT}" \
  -D LAPACK_LIBRARY_NAMES:STRING="liblapack.a" \
  -D TPL_ENABLE_ParMETIS:BOOL=ON \
  -D ParMETIS_INCLUDE_DIRS:PATH="${METIS_ROOT}" \
  -D ParMETIS_LIBRARY_NAMES:STRING="parmetis" \
  -D ParMETIS_LIBRARY_DIRS:PATH="${METIS_ROOT}" \
  -D TPL_ENABLE_METIS:BOOL=ON \
  -D METIS_INCLUDE_DIRS:PATH="${METIS_ROOT}/METISLib" \
  -D METIS_LIBRARY_NAMES:STRING="metis" \
  -D METIS_LIBRARY_DIRS:PATH="${METIS_ROOT}" \
  -D Trilinos_EXTRA_LINK_FLAGS:STRING="$EXTRA_LINK_FLAGS" \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
  -D Amesos_ENABLE_SuperLUDist:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_AztecOO:BOOL=ON \
  -D AztecOO_ENABLE_Teuchos:BOOL=ON \
  -D Trilinos_ENABLE_Didasko:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_Galeri:BOOL=ON \
  -D Trilinos_ENABLE_Ifpack:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
  -D Trilinos_ENABLE_Triutils:BOOL=ON \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  -D CMAKE_Fortran_FLAGS:STRING="-O5 -funroll-all-loops -fPIC" \
  -D CMAKE_C_FLAGS:STRING="-O3 -fPIC -funroll-loops -march=native" \
  -D CMAKE_CXX_FLAGS:STRING="-O3 -fPIC -funroll-loops -ffast-math -march=native -DMPICH_IGNORE_CXX_SEEK" \
 $EXTRA_ARGS \
 ${TRILINOS_SOURCE}
```
Build Trilinos with:
```
sh do-configure.sh
make
make install
```

#### Compiling the `TrilinosApplication`
Afterwards add the following lines to your configure-script:
```Shell
-DTRILINOS_APPLICATION=ON \
```

### Testing MPI compilation
To to tests the compilation, you can execute a simple python script containing these lines:

```
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *
```

### Running Kratos in MPI
To execute Kratos in parallel, you can for example use the `mpiexec` command (running with 4 processes):

`mpiexec -np 4 python3 MainKratos.py --using-mpi`

<!-- ALL THIS IS OUTDATED AND SHOULD NOT BE USED. COMMENTED BECAUSE MAY BE OF SOME USAGE AS FUTURE REFERENCE
In this section we are going to go through the process of compiling a basic version of Kratos Multiphysics under linux environments. Specifically, we explain how to compile in the current Ubuntu LTS (16.04 at this moment), with the latest checked libraries. A basic knowledge of Linux is assumed ( execute commands, create directories, etc...).

_Please notice that this guide is introductory and many of the steps and software used here can be customized. The only hard restriction is to use a C++11 compatible compiler._

## Git

* Objectives:
 * Install git
 * Get Kratos Multiphysics source code

The first thing you will need is the Kratos Multiphysics source code. To download the code you will have to use a git. You can install the default git by using this command:

```Shell
sudo apt-get install git
```

Once git is installed you can fetch the code by using these commands:

```Shell
git clone https://github.com/KratosMultiphysics/Kratos Kratos
```

## Dev Packages

* Objectives:
 * Get Python3-dev
 * Get G++
 * Get Fortran compiler
 * Get LIBBLAS and LIBLAPACK

You will need a series of packages with some kratos dependencies. The command below will install all the packages needed. It will allow you to compile with python2 or python3 so you will be able to chose later ( we recommend python3 )

```Shell
sudo apt-get install python-dev python3-dev gcc g++ gfortran libblas-dev liblapack-dev
```

## Boost
* Objectives:
 * Download Boost

The next step will consist in downloading Boost. Kratos Multiphysics needs Boost libraries to support some of its functions ans you can found them here: http://www.boost.org/users/download/.

We recommend you to use boost 1.67 or newer, earlier versions may cause compiling errors.

It's very important to add the correct path to the boost library in the configure.sh by setting the variable -DBOOST_ROOT. You will see an example in the configure section.

## CMake

* Objectives:
 * Install CMake

CMake is the tool used to compile kratos. To install it, the first option is to execute the following command. Please ensure that you '''use CMake 3.4 or higher''':

```
sudo apt-get install cmake
```

:warning: **If there is no other alternative but CMake 2.8 or earlier please refer to [old verions of CMake](Linux-install#old-version-of-cmake)**

## Configure

* Objectives:
 * Configure Kratos for the first time compilation

In order to compile kratos for the first time you will need to configure the project.  First, navigate to your kratos/cmake_build folder and make a copy of the template file:

```Shell
cp example_configure.sh.do_not_touch configure.sh
```

Then, open configure.sh with any text editor and modify the lines that tell cmake where some components are located.
You will need to provide at least '''BOOST_ROOT'' and '''PYTHON_EXECUTABLE'''.

Kratos has moved to C++11 recently, Please mind to add the "-std=c++11" to your compiler of choice. If you follow the example below, it is already present ( notice the flag in CMAKE_CXX_FLAGS, highlighted in bold)

For example, in ubuntu it will look something like:

```
cmake ..                                                                      \
-DCMAKE_C_COMPILER=/usr/bin/gcc                                               \
-DCMAKE_INSTALL_RPATH="/home/example/kratos/libs"                             \
-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE                                      \
-DCMAKE_CXX_COMPILER=/usr/bin/g++                                             \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3"                                     \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11"                      \
-DBOOST_ROOT="~/compiled_libraries/boost_1_67_0"                              \
-DPYTHON_EXECUTABLE="/usr/bin/python3"                                        \
-DMESHING_APPLICATION=ON                                                      \
-DEXTERNAL_SOLVERS_APPLICATION=ON                                             \

// More options ( do not include this line )
```

Note that the ```\``` must be the last character in the line. Even an space after it will cause an error! (and the returned message is completely misleading, so be careful with this!!)

Notice that you can also turn ON/OFF parts of the code according to your necessities:

```CMake
-DSTRUCTURAL_MECHANICS_APPLICATION=ON/OFF
```

:warning: Cmake requires all definitions in a single line! The line concatenation character ```\``` therefore MUST NOT be followed by any whitespace in the same line as this would prevent the cmake from running the lines below

## Compile

* Objectives:
 * Compile kratos.

If you followed all the steps correctly, compiling kratos should be as easy as executing the configure script:

```Shell
sh configure.sh
```

Please, notice that kratos is big and the compilation process can easily take 1 or 2 hours, depending on which applications are being compiled. A typical compilation process with the default configuration takes approximately 45 minutes with a i7 / 8GB Ram computer.

## Setting up your enviroment

* Objectives:
 * Tell Linux how to execute kratos

Once Kratos is compiled, you will have to tell the OS where to find the libraries. You can do that
by executing these commands. Notice that '''you have to put the same path as in the section "Configure"'''

```Shell
echo "export PYTHONPATH=$PYTHONPATH:/path/to/my/kratos/installation" >> $HOME/.bashrc
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/my/kratos/installation/libs" >> $HOME/.bashrc
```

If you have enabled the embedded python option -DINSTALL_EMBEDDED_PYTHON=ON, you can also add

```Shell
echo "export PATH=$PATH:/path/to/my/kratos/installation" >> $HOME/.bashrc
```

In order to have the "runkratos" available as a regular command.

Now each time you open a terminal these commands will be executed and the paths set automatically.
If you don't want to reset your terminal the first time, just execute:

```Shell
source ~/.bashrc
```

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
               Multi-Physics 3.3.11016
```

# Troubleshooting

In this section we provide will try to provide solutions to the most common problems and issues that may appear during the compilation process

## Cannot find KratosMultiphysics

Make sure that ```LD_LIBRARY_PATH``` and ```PYTHONPATH``` are pointing to your kratos folder.
To make sure that the variables are set correctly, you can always print their value from the terminal by typing:
```Shell
echo $LD_LIBRARY_PATH
```
```Shell
echo $PYTHONPATH
```
```Shell
echo $PATH
```

Some shell interpreters have issues with the separator token. If environment variables above listed are correctly set, try to remove the ":" if kratos is the only path there.

```Shell
# This may cause problem
echo $PYTHONPATH
:path/to/kratos

# Should be
export PYTHONPATH=path/to/kratos
echo $PYTHONPATH
path/to/kratos
```

## I am getting Python link errors

This errors appear if the version of python used to compile boost is not the same as the one
used by Kratos.

There are several causes that may be causing this. Pleas try the following:

### Python version mismatch
The most provable reason for this error to happend is a missmatch between the python versions used by Kratos and Boost. Please, double check you have the same version of python in the projet-config.jam (boost) and configure.sh (Kratos) files.

### Old version of CMake
If the error remains, please check that CMake version is 3.0.2 or newer. If it is not, it will be unable to load python 3.4.
To solve the error please **upgrade to CMake 3.0.2**.

It has been observed that compiling with IDE's ( QTCreator, Netbeans, ...) sometimes causes this error as well.
If you are experiencing this problem, try to modify the configure.sh script and replace '''cmake''' by the absolute path of CMake 3.0.2:

```
/path/to/correct/cmake ..                                                     \
-DCMAKE_C_COMPILER=/usr/bin/gcc                                               \
...
```

If for some reason you have to use an older version of CMake you can manually add support for python 3.4 by adding the version to these files:

```
/usr/share/cmake-2.8/Modules/FindPythonLibs.cmake (line 41)
/usr/share/cmake-2.8/Modules/FindPythonInterp.cmake (line 36)
```

Please add the version at the begining of the list:

```CMake
SET(_PYTHON1_VERSIONS 1.6 1.5)
SET(_PYTHON2_VERSIONS 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0)
SET(_PYTHON3_VERSIONS 3.4 3.3 3.2 3.1 3.0)
```

## I am getting lots of warnings when I compile Kratos

It is known that in some cases warnings appear while including boost files due to the fact that the flag **"-Wunused-local-typedefs"** is set by default in gcc.

This does not have any impact on the final code, but if you want a cleaner output you can add the flag **"-Wno-unused-local-typedefs"** to the configuration files:

```CMake
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -Wno-unused-local-typedefs"      \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 -Wno-unused-local-typedefs"          \
```

# Compiling with MPI
If you want to compile the MPI-version of Kratos, you need to at least to
* Have a MPI installation (e.g. [OpenMPI](http://www.open-mpi.de/) or [Intel MPI](https://software.intel.com/en-us/intel-mpi-library)). There are no known restrictions on which version to use.
* Set the flag `MPI_NEEDED` to `ON` in the configure-script

This will compile the core of Kratos with MPI (i.e. compile the MPI-Core) and provide basic features of MPI. Check [here](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/mpi/tests/test_mpi_data_communicator_python.py) to see which functionalities are exposed to Python.

If you want to run Kratos in MPI together with e.g. the `FluidDynamicsApplication` or the `StructuralMechanicsApplication`, then you also need [Trilinos](https://trilinos.org/) for solving the distributed system of equations and [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for partitioning of the problem.
In order to use them it is required to compile the `TrilinosApplication` and the `MetisApplication` respectively. These applications require to install the corresponding libraries first.

### MetisApplication
#### Installing METIS using the package-manager
On Ubuntu 18.04, the following command installs the necessary files:
```Shell
sudo apt-get install libmetis-dev
```
#### Compiling METIS yourself
Compiling is straight-forward using the build-instructions that come with the download.

#### Compiling the `MetisApplication`
Afterwards add the following lines to your configure-script:
```Shell
-DMETIS_APPLICATION=ON \
-DUSE_METIS_5=ON \
```

### TrilinosApplication

#### Installing Trilinos using the package-manager
On Ubuntu 18.04, the following command installs the necessary files:
```Shell
sudo apt-get install trilinos-all-dev
```
Add the following to your configure-script:
```Shell
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \
-DTRILINOS_LIBRARY_PREFIX="trilinos_" \
```
#### Compiling Trilinos yourself
Some information on how to compile Trilinos (version 12.10.1) is provided in the following, but this might not be applicable for all systems. Please consult the build-instructions of Trilinos itself if needed or open an issue.

Add the following to your configure-script if you compile Trilinos yourself:
```Shell
-DTRILINOS_ROOT="path/to/trilinos" \
```

Depending on the configuration it might also be necessary to specify `TRILINOS_SOURCE` for CMake to find the correct trilinos-libraries.
***

It is recommended to create a build directory in the root folder of Trilinos. Inside this folder, create a `do-configure.sh` to compile Trilinos.

If you use LAPACK 3.6 or newer, then you have to apply the following fixes: (sources: ​[1](https://github.com/gahansen/Albany/wiki/ALCF-Vesta), [​2](https://www.dealii.org/developer/external-libs/trilinos.html)):

In `trilinos-12.10.1-Source/packages/epetra/src/Epetra_LAPACK_wrappers.h`

Change
```
#define DGGSVD_F77  F77_BLAS_MANGLE(dggsvd,DGGSVD)
#define SGGSVD_F77  F77_BLAS_MANGLE(sggsvd,SGGSVD)
```
to
```
#define DGGSVD_F77  F77_BLAS_MANGLE(dggsvd3,DGGSVD3)
#define SGGSVD_F77  F77_BLAS_MANGLE(sggsvd3,SGGSVD3)
```

This is an example for a `do-configure.sh`. Of course, the paths have to be adjusted
```
TRILINOS_ROOT="${HOME}/software/kratos/trilinos"
EXTRA_LINK_FLAGS=""
EXTRA_ARGS=$@
MPI_ROOT="${HOME}/software/ompi"
METIS_ROOT="{HOME}/software/kratos/parmetis/ParMetis-3.2.0"
LAPACK_ROOT="${HOME}/software/kratos/lapack/lapack-3.7.0"

rm CMakeCache.txt
rm -rf CMakeFiles
rm *.cmake

cmake \
  -D CMAKE_INSTALL_PREFIX:FILEPATH="${TRILINOS_ROOT}" \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D Trilinos_SOURCE_DIR="${TRILINOS_SOURCE}" \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_Fortran_COMPILER=mpif90 \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D Trilinos_ENABLE_CXX11=ON\
  -D MPI_BASE_DIR="${MPI_ROOT}" \
  -D MPI_INCLUDE_DIRS:PATH="${MPI_ROOT}/include" \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D BLAS_LIBRARY_DIRS:FILEPATH="${LAPACK_ROOT}" \
  -D BLAS_LIBRARY_NAMES:STRING="librefblas.a" \
  -D LAPACK_LIBRARY_DIRS:FILEPATH="${LAPACK_ROOT}" \
  -D LAPACK_LIBRARY_NAMES:STRING="liblapack.a" \
  -D TPL_ENABLE_ParMETIS:BOOL=ON \
  -D ParMETIS_INCLUDE_DIRS:PATH="${METIS_ROOT}" \
  -D ParMETIS_LIBRARY_NAMES:STRING="parmetis" \
  -D ParMETIS_LIBRARY_DIRS:PATH="${METIS_ROOT}" \
  -D TPL_ENABLE_METIS:BOOL=ON \
  -D METIS_INCLUDE_DIRS:PATH="${METIS_ROOT}/METISLib" \
  -D METIS_LIBRARY_NAMES:STRING="metis" \
  -D METIS_LIBRARY_DIRS:PATH="${METIS_ROOT}" \
  -D Trilinos_EXTRA_LINK_FLAGS:STRING="$EXTRA_LINK_FLAGS" \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
  -D Amesos_ENABLE_SuperLUDist:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_AztecOO:BOOL=ON \
  -D AztecOO_ENABLE_Teuchos:BOOL=ON \
  -D Trilinos_ENABLE_Didasko:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_Galeri:BOOL=ON \
  -D Trilinos_ENABLE_Ifpack:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
  -D Trilinos_ENABLE_Triutils:BOOL=ON \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  -D CMAKE_Fortran_FLAGS:STRING="-O5 -funroll-all-loops -fPIC" \
  -D CMAKE_C_FLAGS:STRING="-O3 -fPIC -funroll-loops -march=native" \
  -D CMAKE_CXX_FLAGS:STRING="-O3 -fPIC -funroll-loops -ffast-math -march=native -DMPICH_IGNORE_CXX_SEEK" \
 $EXTRA_ARGS \
 ${TRILINOS_SOURCE}
```
Build Trilinos with:
```
sh do-configure.sh
make
make install
```

#### Compiling the `TrilinosApplication`
Afterwards add the following lines to your configure-script:
```Shell
-DTRILINOS_APPLICATION=ON \
```

### Testing MPI compilation
To to tests the compilation, you can execute a simple python script containing these lines:

```
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *
```

### Running Kratos in MPI
To execute Kratos in parallel, you can for example use the `mpiexec` command (running with 4 processes):

`mpiexec -np 4 python3 MainKratos.py --using-mpi`

Additional information regarding the compilation can be found [here](https://github.com/KratosMultiphysics/Kratos/blob/master/cmake_build/README.txt)
-->