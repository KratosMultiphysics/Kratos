# Basic use of the **mmg3d** library

## I/ Implementation
To call the **mmg3d** library, you must:  
  1. build mesh and sol at MMG5 format;
  2. call the MMG5 library;
  3. get the final mesh and sol.

  Results are saved in the **_mesh.o.mesh_** and **_mesh.o.sol_** files.  

### example0_a  
  We read mesh and solution files (**_cube.mesh_** and **_cube.sol_**) using the **MMG5_loadMesh** and **MMG5_loadMet** functions.
  Results are saved using **MMG5_saveMesh** and **MMG5_saveMet** functions.

### example0_b
  The mesh and solution are hard coded.    
  They are build in MMG5 format using API functions and are recovered by the same way.  
  We show how to recover the mesh/sol by writting it in mesh.o.mesh/sol file.

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