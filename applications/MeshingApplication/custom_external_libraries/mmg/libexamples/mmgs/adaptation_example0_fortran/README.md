# Basic use of the mmgs library for a Fortran call

## I/ Implementation
To call the **mmgs** library, you must:  
  1. build mesh and sol at MMG5 format;
  2. call the MMG5 library;
  3. get the final mesh and sol.

### example0_a  
  We read mesh and solution files using the **MMGS_loadMesh** and **MMGS_loadSol** functions.
  Results are saved using **MMGS_saveMesh** and **MMGS_saveSol** functions.

### example0_b
  The mesh and solution are hard coded.    
  They are build in MMG5 format using API functions and are recovered by the same way.  
  We show how to recover the mesh/sol.

## II/ Compilation
  1. Build and install the **mmgs** shared and static library. We suppose in the following that you have installed the **mmgs** library in the **_$CMAKE_INSTALL_PREFIX_** directory (see the [installation](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) section of the setup guide);
  2. compile the main.c file specifying:
    * the **mmgs** include directory with the **-I** option;
    * the **mmgs** library location with the **-L** option;
    * the **mmgs** library name with the **-l** option;
    * for the static library you must also link the executable with, if used for the **mmgs** library compilation, the scotch and scotcherr libraries and with the math library;
    * with the shared library, you must add the ***_$CMAKE_INSTALL_PREFIX_** directory to your **LD_LIBRARY_PATH**.

> Example 1  
>  Command line to link the application with the **mmgs** static library (we supposed here that the scotch library is installed in the **_$SCOTCH_PATH_** directory):  
> ```Shell
> gfortran -I$CMAKE_INSTALL_PREFIX/include main.c -L$CMAKE_INSTALL_PREFIX/lib -L$SCOTCH_PATH -lmmgs -lscotch -lscotcherr -lm
> ```

> Example 2  
>  Command line to link the application with the **mmgs** shared library:  
> ```Shell
> gfortran -I$CMAKE_INSTALL_PREFIX/include main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmgs
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```
