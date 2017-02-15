# Example of use of libmmg3d5

## I/ Object
  Here, we transform the library in an executable similar to the **mmg3d** application.  
  We provide a mesh and a sol file (SphereIso0.5.mesh/sol) to test the program.

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

## III/ Execution
Because it contains hard coded paths to the mesh and solution files, the tests must be run from a subdirectory of the root of your **mmg** project.

> Example  
> Assuming that your **mmg** project is cloned into the **_mmg_** directory (default case), you can run the test from the **_mmg/build/_** or **_mmg/libexamples_** directories but not from the **_mmg/libexamples/mmg3d/_** folder.