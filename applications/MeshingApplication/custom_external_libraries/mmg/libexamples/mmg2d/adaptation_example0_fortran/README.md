# Basic use of the mmg2d library for a Fortran call

## I/ Implementation
To call the **mmg3d** library, you must:  
  1. build mesh and sol at MMG5 format;
  2. call the MMG5 library;
  3. get the final mesh and sol.

### example0_a  
  We read mesh and solution files using the **MMG2D_loadMesh** and **MMG2D_loadSol** functions.
  Results are saved using **MMG2D_saveMesh** and **MMG2D_saveSol** functions.

### example0_b
  The mesh and solution are hard coded.    
  They are build in MMG5 format using API functions and are recovered by the same way.  
  We show how to recover the mesh/sol.

## II/ Compilation
  1. Build and install the **mmg2d** shared and static library. We suppose in the following that you have installed the **mmg2d** library in the **_$CMAKE_INSTALL_PREFIX_** directory (see the [installation](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) section of the setup guide);
  2. compile the main.F90 file specifying:
    * the **mmg2d** include directory with the **-I** option;
    * the **mmg2d** library location with the **-L** option;
    * the **mmg2d** library name with the **-l** option;
    * with the shared library, you must add the ***_$CMAKE_INSTALL_PREFIX_** directory to your **LD_LIBRARY_PATH**.

> Example 1  
>  Command line to link the application with the **mmg2d** static library
> ```Shell
> gfortran -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg2d main.F90 -L$CMAKE_INSTALL_PREFIX/lib -lmmg2d -lm
> ```

> Example 2  
>  Command line to link the application with the **mmg2d** shared library:  
> ```Shell
> gfortran -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg2d main.F90 -L$CMAKE_INSTALL_PREFIX/lib -lmmg2d
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```
