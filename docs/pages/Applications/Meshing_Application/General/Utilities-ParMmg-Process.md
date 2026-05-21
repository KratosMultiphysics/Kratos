---
title: Utilities ParMmg Process
keywords: 
tags: [Utilities-ParMmg-Process.md]
sidebar: mmg_application
summary: 
---

# Introduction

ParMMG is the MPI parallel version of the remeshing library MMG. In order to use ParMmg, please checkout the guide on how to install and use MMG first:

[MMG remeshing process install guide](%5BUtilities%5D-MMG-Process)

As of now, ParMMG is separate library that depends on MMG and it therefore requires another installation. ParMmg does not have all functionalities from MMG, but it currently features:
* 3D volume remeshing of tetrahedras

# How to install

In order to install ParMmg in our system, we need to clone the repository. 

    git clone https://github.com/MmgTools/ParMmg

Once cloned, we can change directory to ParMmg, and create a build folder where the libraries will be written. 

    cd ParMmg;
    mkdir build;
    cd build;
      
Now, we will create and edit a configure.sh file:

    touch configure.sh

In the configure file we will write the following:

    cmake ..                                              \
    -DCMAKE_BUILD_TYPE=Release                            \
    -USE_POINTMAP=ON                                      \
    -DCMAKE_CXX_FLAGS="-O3 -msse3 -fPIC -fopenmp"         \
    -DCMAKE_C_FLAGS="-O3 -msse3 -fPIC -fopenmp"           \
    -DLIBPARMMG_SHARED=ON                                 \
    -DLIBPARMMG_STATIC=OFF                                \
    -DDOWNLOAD_MMG=OFF                                    \
    -DMMG_DIR="/path/to/mmg"                              \
    -DMMG_BUILDDIR="/path/to/mmg/build"                   \
    -DDOWNLOAD_METIS=OFF                                  \
    -DMETIS_DIR="/usr/include"                            \

    make -j4


Where we have set the path to the MMG installation. Please set it with the correct path in your machine. ParMmg will be linked to a certain version of MMG. For instance, v1.3.0 of ParMmg needs version v5.5.2 of MMG. This can be checked in the file:

    ParMmg/CMakeLists.txt

Specified as `GIT_TAG` in the MMG section. We can now compile:

    sh configure.sh


After compiling, we can now compile Kratos with the following extra variables. Note that Kratos will have to be compiled in MPI with the flag `-DUSE_MPI=ON`

    -DUSE_MPI=ON                      \
    -DINCLUDE_PMMG=ON                 \
    -DPMMG_ROOT=/path/to/parmmg/build \

Where we have set the path to ParMmg build folder with the flag `-DPMMG_ROOT`. Note that this will assume that libraries are located in `ParMmg/build/lib` and the include folder in `ParMmg/build/include`. If this is not the case, you can specify the custom library and include dir with the following flags (remove PMMG_ROOT if this is the case):

    -DUSE_MPI=ON                                                  \
    -DINCLUDE_PMMG=ON                                             \
    -DPMMG_INCLUDE_DIR=/path/to/custom/parmmg/build/include/      \
    -DPMMG_LIBRARY=/path/to/custom/parmmg/build/lib/libparmmg.so  \

Lastly, remember to add the `/lib` folder to  your LD_LIBRARY_PATH:

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/parmmg/build/lib

You can add this command to your `.bashrc` to have it always on login. 