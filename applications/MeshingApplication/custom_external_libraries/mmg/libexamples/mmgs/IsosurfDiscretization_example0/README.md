# Example of basic use of library **libmmgs** for a level-set discretization test case

## I/ Implementation
  We read the mesh and level-set values at the mesh nodes and then, we discretize the implicit function and optimize the mesh.

## II/ Compilation
  1. Build and install the **mmgs** shared and/or static library. We suppose in the following that you have installed the **mmgs** library in the **_$CMAKE_INSTALL_PREFIX_** directory (see the [installation](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) section of the setup guide);
  2. compile the main.c file specifying:
    * the **mmgs** include directory with the **-I** option;
    * the **mmgs** library location with the **-L** option;
    * the **mmgs** library name with the **-l** option;
    * for the static library you must also link the executable with the math library and, if used for the **mmgs** library compilation, the **scotch** and **scotcherr** libraries;
    * with the shared library, you must add the ***_$CMAKE_INSTALL_PREFIX_** directory to your **LD_LIBRARY_PATH**.

> Example 1  
>  Command line to link the application with the **mmgs** static library (we supposed here that the scotch library is installed in the **_$SCOTCH_PATH_** directory):  
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/mmg/mmgs main.c -L$CMAKE_INSTALL_PREFIX/lib -L$SCOTCH_PATH -lmmgs -lscotch -lscotcherr -lm
> ```

> Example 2  
>  Command line to link the application with the **mmgs** shared library:  
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/mmg/mmgs main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmgs
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```
