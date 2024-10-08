# Basic use of the **mmg2d** library

## Adaptation example 3
No mesh or sol file necessary, only compilation of main.c.

## Compilation
  1. Build and install the **mmg2d** shared and static library. We suppose in the following that you have installed the **mmg2d** library in the **_$CMAKE_INSTALL_PREFIX_** directory (see the [installation](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) section of the setup guide);
  2. compile the main.c file specifying:
    * the **mmg2d** include directory with the **-I** option;
    * the **mmg2d** library location with the **-L** option;
    * the **mmg2d** library name with the **-l** option;
    * with the shared library, you must add the ***_$CMAKE_INSTALL_PREFIX_** directory to your **LD_LIBRARY_PATH**.

>  Command line to link the application with the **mmg2d** static library
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg2d main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmg2d -lm
> ```

>  Command line to link the application with the **mmg2d** shared library:  
> ```Shell
> gcc -I$CMAKE_INSTALL_PREFIX/include/mmg/mmg2d main.c -L$CMAKE_INSTALL_PREFIX/lib -lmmg2d -lm
> export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
> ```
