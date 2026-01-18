# Basic use of the mmg library for an adaptation testcase in C++

## I/ Implementation
To call the **mmg** library, you must:  
  1. build mesh and sol at MMG2D, MMG3D or MMGS format;
  2. call the MMG2D, MMGS or MMG3D library;
  3. get the final mesh and sol.

  We read mesh and solution files using the **MMG<2D/S/3D>_loadMesh** and **MMG<2D/S/3D>_loadSol** functions.
  Results are saved using **MMG<2D/S/3D>_saveMesh** and **MMG<2D/S/3D>_saveSol** functions.

## II/ Compilation
  1. Build and install the **mmg** shared and static library. We suppose in the following that you have installed the **mmg** library in the **_$CMAKE_INSTALL_PREFIX_** directory (see the [installation](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) section of the setup guide);
  2. compile the main.cpp file specifying:
    * the **mmg** include directory with the **-I** option;
    * the **mmg** library location with the **-L** option;
    * the **mmg** library name with the **-l** option;
    * for the static library you must also link the executable with, if used for the **mmg** library compilation, the scotch and scotcherr libraries and with the math library;
    * with the shared library, you must add the ***_$CMAKE_INSTALL_PREFIX_** directory to your **LD_LIBRARY_PATH**.

> Example 1  
>  Command line to link the application with the **mmg** static library (we supposed here that the scotch library is installed in the **_$SCOTCH_PATH_** directory):  
> ```Shell
> g++ -I$CMAKE_INSTALL_PREFIX/include/mmg/ main.c -L$CMAKE_INSTALL_PREFIX/lib -L$SCOTCH_PATH -lmmg -lscotch -lscotcherr -lm
> ```

> Example 2  
>  Command line to link the application with the **mmg** shared library:  
> ```Shell
> g++ -I$CMAKE_INSTALL_PREFIX/include/mmg/ main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmg
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```
