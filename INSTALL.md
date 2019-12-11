# Contents
* [Basic Configuration](#basic-configuration)
* [Adding Applications](#basic-configuration)
* [Advanced Configuration](#advanced-configuration)
  * [Common Flags](#common-flags)
  * [Compilation Performance](#compilation-performance)
  * [Parallelism](#parallelism)
  * [External Libraries](#external-libraries)
    * [Feast](#feast)
    * [Metis](#metis)
    * [Trilinos](#trilinos)
* [Applications](#applications)


## Basic Configuration

You can find the new kratos configuration file in Kratos root folder (`configure.sh` for linux and `configure.bat` for win).
Out of the box Kratos will try to find all necessary libraries in your system automatically, but we recommend to take a look at the following configuration options:

`KRATOS_BUILD_TYPE`

Compilation Type. Options are `Release`,`RelWithDebInfo`,`Debug`,`FullDebug`,`Custom`

**Release**: Full Release with maximum optimization options and no debug Info.

**RelWithDebInfo**: Full Release with optimization and debug info. Adecuate to debug simple problems without losing performance.

**Debug**: Debug build with no optimization flags.

**FullDebug**: Debug build with no optimization falgs, extended debug info and extremly low performance.

**Custom**: No flags are automatically added.

`PYTHON_EXECUTABLE`

Path to the python executable that Kratos will use. We recommend that you manually set this in case you have multiple versions of python in the system.
Ubuntu users need to be extra careful with this as default versions tends to be Python2, while Kratos is preferably compiled with Python3

`BOOST_ROOT`

Path to boost root directory.

## Adding Applications

In order to add an application you can use the provided macro (`add_app [PATH]` for Linux, `CALL :add_app [PATH]` for Win) along with the route folder of the application that you want to compile. Several examples are provided in the configuration files.

Its now also possible to compile applications outside kratos source dir:

Linux:
```shell
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication    
add_app ${KRATOS_APP_DIR}/FluidDynamicApplication
add_app /home/username/development/ExternalApplication  # Example of external Application
```

Windows:
```shell
CALL :add_app %KRATOS_APP_DIR%/ExternalSolversApplication    
CALL :add_app %KRATOS_APP_DIR%/FluidDynamicApplication
CALL :add_app C:/users/username/development/ExternalApplication  # Example of external Application
```

## Advanced Configuration

It is also possible to use more advanced configuration options. To use any of these options please add them directly to the cmake configuration line just as any other flag

### Common Flags

`-DCMAKE_C_COMPILER=String`

Path to the C compiler. Overrides `CC` environment variable

`-DCMAKE_CXX_COMPILER=String`

Path to the C++ compiler. Overrides `CXX` environment variable

`-DCMAKE_INSTALL_PREFIX=String`

Install path for Kratos. If not set the installation will be done in `bin/[build_type]`

`-DCMAKE_C_FLAGS=String`

User defined flags for the C compiler.

`-DCMAKE_CXX_FLAGS=String`

User defined flags for the C++ compiler. 

`-DBOOST_ROOT=String`

Root directory for boost. Overrided by BOOST_ROOT environmental variable if defined.

`-DPYTHON_EXECUTABLE=String`

Python executable to be used. It is recommended to set this option if more than one version of python is present in the system (For example while using ubuntu). Overrided by PYTHON_EXECUTABLE environmental variable if defined.

`-DINSTALL_PYTHON_USING_LINKS=ON/OFF`

Enables or Disables(default) the creation of links to the source python files instead of copying them during the installation. Useful while developing code. Does not work on windows.

`-DINSTALL_RUNKRATOS=ON/OFF`

Enables(Default) or Disables the compilation of the embedded python interpreter (aka Runkratos).

`-DKRATOS_BUILD_TESTING=ON/OFF`

Enables(Default) or Disables the compilation of the C++ unitary tests for Kratos and Applications.

### Compilation Performance
`-DUSE_COTIRE=ON/OFF`

Enables or Disables(default) the use of [cotire](https://github.com/sakra/cotire) to speedup compilation by using unitary builds.
Please notice that enabling this options can greatly increase the amount of memory needed to compile some targets, specially if combined with -jx.

In order to install and compile with this switch please use:
```shell
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target all_unity -- -j1
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install/fast -- -j1 
```

Instead of the regular install target.

### Parallelism
`-DUSE_MPI=ON/OFF`

Enables or Disables(default) the modules and code for mpi. This option is needed if you want to compile Trilinos, Metis, etc...

### External libraries
#### Feast
`-DINCLUDE_FEAST`

Enables or Disables(default) the use of FEAST

#### Metis

`-DUSE_METIS_5=ON/OFF`

Specifies if the metis version is 5 or greater (OFF by default).

`-DMETIS_ROOT_DIR=String`

Root directory for Metis library

#### Trilinos

`-DTRILINOS_ROOT=String`

Root directory for Trilinos library.

`-DTRILINOS_INCLUDE_DIR=String`

Not required if `TRILINOS_ROOT` is set. Path to trilinos include dir.

`-DTRILINOS_LIBRARY_DIR=String`

Not required if `TRILINOS_ROOT` is set. Path to trilinos library dir.

`-DTRILINOS_PREFIX=String`
Indicates the prefix of the trilinos libraries in case they have:
```
libepetra.so          -> No prefix
libtrilinos_epetra.so -> -DTRILINOS_PREFIX="trilinos_"
```

## Applications

Specific compilation information about applications can be found in their own directories:

- [EigenSolvers Application](applications/EigenSolversApplication/README.md#build-instructions)
- [HDF5 Application](applications/HDF5Application/README.md#build-instructions)
- [MultilevelMontecarlo Application](applications/MultilevelMonteCarloApplication/README.md#external-libraries)
- [Poromechanics Application](applications/PoromechanicsApplication/README.md#how-to-use-mpi-in-poromechanics-application)
- [Trilinos Application (Aditional notes)](applications/TrilinosApplication/README.md#notes-for-compilation)