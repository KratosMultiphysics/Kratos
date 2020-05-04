# Input and output for multiple solutions (we do not call any mmg library here)

## I/ Implementation
  1. Read a mesh and a solution file with multiple solutions at MMG5 format
     (using the **MMG3D_loadMesh** and **MMG3D_loadAllSols** functions).
  2. Use the Mmg setters to set this solutions in another solution structure.
  3. Use the Mmg getters to get the solutions in another new solution structure.
  4. Save this last solution structure and the multiple solutions stored using
     the **MMG3D_saveAllSols** function.

## II/ Compilation
  1. Build and install the **mmg3d** shared and static library. We suppose in the following that you have installed the **mmg3d** library in the **_$CMAKE_INSTALL_PREFIX_** directory (see the [installation](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) section of the setup guide);
  2. compile the main.c file specifying:
    * the **mmg3d** include directory with the **-I** option;
    * the **mmg3d** library location with the **-L** option;
    * the **mmg3d** library name with the **-l** option;
    * with the shared library, you must add the ***_$CMAKE_INSTALL_PREFIX_** directory to your **LD_LIBRARY_PATH**.

> Example 1  
>  Command line to link the application with the **mmg3d** static library
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg3d main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmg3d -lm
> ```

> Example 2  
>  Command line to link the application with the **mmg3d** shared library:  
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg3d main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmg3d
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```
