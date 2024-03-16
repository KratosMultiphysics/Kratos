---
title: the PaStiX solver in Kratos
keywords: 
tags: [How-to-use-the-PaStiX-solver-in-Kratos.md]
sidebar: kratos_for_developers
summary: 
---

# Contents
1. [How to use the PaStiX solver in Kratos:][pastix1]
    1. [Compile the Scotch library (shared memory parallelism)][pastix1a]
    2. [Compile the PaStiX library (shared memory parallelism)][pastix1b]
    3. [Compile the Kratos with PaStiX (shared memory parallelism)][pastix1c]
    4. [Additional hints][pastix1d]

[pastix1]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-use-the-PaStiX-solver-in-Kratos#how-to-use-the-pastix-solver-in-kratos
[pastix1a]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-use-the-PaStiX-solver-in-Kratos#compile-the-scotch-library-shared-memory-parallelism
[pastix1b]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-use-the-PaStiX-solver-in-Kratos#compile-the-pastix-library-shared-memory-parallelism
[pastix1c]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-use-the-PaStiX-solver-in-Kratos#compile-the-kratos-with-pastix-shared-memory-parallelism
[pastix1d]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-use-the-PaStiX-solver-in-Kratos#additional-hints

# How to use the PaStiX solver in Kratos 

The **PaStiX (Parallel Sparse matriX package)** is a open-source scientific library that provides a high performance parallel solver for very large sparse linear systems based on direct methods. Numerical algorithms are implemented in single or double precision (real or complex) using *LLt*, *LDLt* and *LU* with static pivoting (for non symmetric matrices having a symmetric pattern). This solver provides also an adaptive blockwise *iLU(k)* factorization that can be used as a parallel preconditioner using approximated supernodes to build a coarser block structure of the incomplete factors ([link](http://pastix.gforge.inria.fr/files/README-txt.html)).

It can be linked to the Kratos by the **ExternalSolversApplication**. The following steps have to be made in order to use be able to use it in *Ubuntu* version 14.04 LTS 64 bit or 16.04 LTS 64 bit. 

## Compile the Scotch library (shared memory parallelism)

*PaStiX* needs the *Scotch* library as a prerequisite to work.

You can download version 6.0.4 here: [http://gforge.inria.fr/projects/scotch/](http://gforge.inria.fr/projects/scotch/) 

Copy it to a folder of your choice. In the following we are using `~/compiled_libraries` as a reference.

Extract the files by the command: 

```console
tar -xf scotch_6.0.4.tar.gz
```

Navigate into the source directory by: 

```console
cd scotch_6.0.4/src
```

Copy the correct **Makefile.inc** for the version of your system to the source directory. Ín the case of the *Ubuntu* versions mentioned above it is the **Makefile.inc.x86-64_pc_linux2**:  (**NOTE**: if you want to use shared libraries **.so** you will copy the **Makefile.inc.x86-64_pc_linux2.shlib**)

```console
cp ./Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc
```

Add -fPIC to CFLAGS in **Makefile.inc** for shared libraries. Optional libraries may be installed (see scotch_6.0.4/INSTALL.txt): 

```console
sudo apt-get install flex bison zlib1g-dev
```

Now *Scotch* can be compiled by (we are still in the `/src` directory): 

```console
make scotch
```

And can be installed by: 

```console
make prefix=~/path-to-local-scotch-install-dir install
```

## Compile the PaStiX library (shared memory parallelism)

Download the *PaStiX* version 5.2.3 from:  [https://gforge.inria.fr/frs/?group_id=186](https://gforge.inria.fr/frs/?group_id=186) 

Copy it to a folder of your choice. In the following we are using `~/compiled_libraries` as a reference.

Extract the files by the command: 

```console 
tar -xf pastix_5.2.3.tar.bz2
```

Navigate into the source directory by: 

```console
cd pastix_5.2.3/src
```

The different config files can be found in the config folder. You have to copy the one fitting to your system to the `/src` folder. In the case of the systems mentioned above it ist the **LINUX-GNU.in** file. Copy it to the source directory as **config.in** by: 

```console
cp ./config/LINUX-GNU.in config.in
```

You have to edit the following blocks of the **config.in** file: 

```console
###################################################################                               
#                  SETTING INSTALL DIRECTORIES                    #                               
###################################################################                               
ROOT          = ${HOME}/path-to-local-pastix-install-dir
INCLUDEDIR    = ${ROOT}/include
LIBDIR        = ${ROOT}/lib
BINDIR        = ${ROOT}/bin
```

```console
###################################################################                             
#                  SHARED LIBRARY GENERATION                      #                             
###################################################################                             
SHARED=1                                                                                       
SOEXT=.so                                                                                      
SHARED_FLAGS =  -shared -Wl,-soname,__SO_NAME__                                                
CCFDEB       := ${CCFDEB} -fPIC
CCFOPT       := ${CCFOPT} -fPIC
CFPROG       := ${CFPROG} -fPIC
```

```console
###################################################################
#                          MPI/THREADS                            #
###################################################################
# Uncomment the following lines for sequential (NOMPI) version
VERSIONMPI  = _nompi
CCTYPES    := $(CCTYPES) -DFORCE_NOMPI
MPCCPROG    = $(CCPROG)
MCFPROG     = $(CFPROG)
MPCXXPROG   = $(CXXPROG)
...some code...
#CCPASTIX   := $(CCPASTIX) -DCUDA_SM_VERSION=20
#NVCCOPT    := $(NVCCOPT) -arch sm_20
```

```console
###################################################################                             
#                      GRAPH PARTITIONING                         #                             
###################################################################
...some code...
# Uncomment on of this blocks
#scotch
CCPASTIX   := $(CCPASTIX) -I$(SCOTCH_INC) -DWITH_SCOTCH
EXTRALIB   := $(EXTRALIB) -L$(SCOTCH_LIB) -lscotch -lscotcherrexit
#ptscotch
#CCPASTIX   := $(CCPASTIX) -I$(SCOTCH_INC) -DDISTRIBUTED -DWITH_SCOTCH
#if scotch >= 6.0
#EXTRALIB   := $(EXTRALIB) -L$(SCOTCH_LIB) -lptscotch -lscotch -lptscotcherrexit
#else
#EXTRALIB   := $(EXTRALIB) -L$(SCOTCH_LIB) -lptscotch -lptscotcherrexit
```

```console
###################################################################
#                Portable Hardware Locality                       #
###################################################################
# By default PaStiX uses hwloc to bind threads,
# comment this lines if you don't want it (not recommended)
#HWLOC_HOME = ~/compiled_libraries/hwloc-1.11.5 #/opt/hwloc/
#HWLOC_INC  = $(HWLOC_HOME)/include
#HWLOC_LIB  = $(HWLOC_HOME)/lib
#CCPASTIX   := $(CCPASTIX) -I$(HWLOC_INC) -DWITH_HWLOC
#EXTRALIB   := $(EXTRALIB) -L$(HWLOC_LIB) -lhwloc
```

**NOTE**: In order to compile the shared libraries **.so** instead of the static ones with `-fPIC` as above, you will need to modify:

Now you can compile the *PaStiX* library (we are still in the /src directory): 

```console
make SCOTCH_HOME=~/path-to-local-scotch-install-dir
```

You can install it to the directory defined in the config.in file by: 

```console
make SCOTCH_HOME=~/path-to-local-scotch-install-dir install
```

## Compile the Kratos with PaStiX (shared memory parallelism)

Finally we can compile the Kratos with the *PaStiX* library. In order to do this, you have to add the following lines to you ' configure.sh` file: 

```console
-DINCLUDE_PASTIX=ON                                       \
-DPASTIX_INSTALL_DIR="/path/to/your/pastix/installation" \ (e.g. "~/software/pastix_5.2.3-install/install") 
-DSCOTCH_INSTALL_DIR="~/path-to-local-scotch-install-dir/lib" \ (e.g. "~/software/scotch_6.0.4-install/lib")
```

**NOTE**: Remember to add the shared libraries to your runtime path.

**NOTE**: Be careful, if you already installed the **Scotch** libraries using the packages available in *Ubuntu* maybe you can have some problems to compile *Kratos*, in that case we recommend to install **Scotch** directly in the system by running `sudo make install`.

You can now use the direct and iterative *PaStiX* solver by adding the following block to your **ProjectParameters.json** file: 

```json
"solver_type": "PastixSolver",
"solution_method": "Direct",
"tolerance": 0.000001,
"max_iteration": 100,
"gmres_krylov_space_dimension": 100,
"ilu_level_of_fill": 1,
"is_symmetric": false,
"verbosity": 0,
"scaling": false,
"block_size": 1,
"use_block_matrices_if_possible": true
```

## Additional hints

Sometimes you can find additional issues with the PasTiX compilation. These hints define some possible solution to these problems.

### Compile OpenBLAS

Some solver can became slower once you compile PaStiX, this is due to the incompatibilty of the [BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms) BLAS spaces and the OpenMP. One alternative is to compile OpenBLAS in order to tackle this problem.

First we go the OpenBLAS [webpage](http://www.openblas.net/) webpage, and download the last version. 

Extract the files by the command: 

```console
 tar -xf v0.2.19.tar.gz
```

We can download it too using git:

```console
git clone https://github.com/xianyi/OpenBLAS.git
```

Now we go the folder:

```console
cd OpenBLAS
```

And we compile the library using the following command (`make install` might require `sudo`)

```console
make -j8  USE_THREAD=0 OPENBLAS_NUM_THREADS=1
make install
```

We just need to modify the `config.in` and we recompile PasTiX again: 

```console
###################################################################
#                              BLAS                               #
###################################################################

# Choose Blas library (Only 1)
# Do not forget to set BLAS_HOME if it is not in your environnement
BLAS_HOME=/your_directory/OpenBLAS/
#----  Blas    ----
BLASLIB  = -L${BLAS_HOME} -lopenblas
#---- Gotoblas ----
#BLASLIB  = -L${BLAS_HOME} -lgoto
#----  MKL     ----
#Uncomment the correct line
#BLASLIB  = -L$(BLAS_HOME) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#BLASLIB  = -L$(BLAS_HOME) -lmkl_intel -lmkl_sequential -lmkl_core
#----  Acml    ----
#BLASLIB  = -L$(BLAS_HOME) -lacml
```