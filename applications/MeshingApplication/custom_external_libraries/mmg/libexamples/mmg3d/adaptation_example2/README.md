# Example of advanced use of library libmmg3d

## I/ Implementation
  We read the mesh and solution in files.

  * First we remesh in debug mode:
    * we ask for a minimal size of 0.001, a maximal size of 40, a gradation of 2 and a global hausdorff value (applied on all the boundaries) of 0.1;
    * we save the mesh and solution.

  * Second, we remesh in normal mode, with specified memory and lower verbosity:
    * in addition to previous parameters, we ask that all boundary triangles of ref 36 respect a hasdorff number of 0.01 and all boundary triangles of ref 38 a hasdorff number of 1. For both reference, we set the minimal (resp. maximal) edge size to 0.00001 (resp. 40).  
      The local hausdorff number on ref 38 has no effects because it is higher than the previous value (without local value, we apply global hausdorff (0.1)) and this value is now contained in the metric;
    * we don't save results but we reset the computed metric and reapply the initial constant metric of size 10;
    * we perform the last wave of refinment. Now we can see the effect of the local hausdorff number on ref 38;
    * we save the mesh and solution

## II/ Compilation
  1. Build and install the **mmg3d** shared and static library. We suppose in the following that you have installed the **mmg3d** library in the **_$CMAKE_INSTALL_PREFIX_** directory (see the [installation](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) section of the setup guide);
  2. compile the main.c file specifying:
    * the **mmg3d** include directory with the **-I** option;
    * the **mmg3d** library location with the **-L** option;
    * the **mmg3d** library name with the **-l** option;
    * for the static library you must also link the executable with, if used for the **mmg3d** library compilation, the scotch and scotcherr libraries and with the math library;
    * with the shared library, you must add the ***_$CMAKE_INSTALL_PREFIX_** directory to your **LD_LIBRARY_PATH**.

> Example 1  
>  Command line to link the application with the **mmg3d** static library (we supposed here that the scotch library is installed in the **_$SCOTCH_PATH_** directory):  
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg3d main.c -L$CMAKE_INSTALL_PREFIX/lib -L$SCOTCH_PATH -lmmg3d -lscotch -lscotcherr -lm
> ```

> Example 2  
>  Command line to link the application with the **mmg3d** shared library:  
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg3d main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmg3d
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```
