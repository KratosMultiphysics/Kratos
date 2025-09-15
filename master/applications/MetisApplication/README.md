# Metis Application

The Metis application provides an interface to the
[METIS](https://github.com/KarypisLab/METIS) graph partitioning library,
used within Kratos Multiphysics to generate mesh partitions for MPI runs.

The applicaiton is organized in a series of processes, each one of them providing an interface to a different partitioning algorithm.
The algorithm most commonly used for finite elements is `metis_divide_heterogeneous_input_process.h`,
which generates a partition of the mesh based on the nodal graph and can be used for meshes combining multiple element types.

## Build instructions

### Direct GNU/Linux install

Building the `MetisApplication` requires *METIS* libraries and their dependencies already installed on the system. The easiest way to get them is by the package manager of the *GNU/Linux* distribution. For example, in *Ubuntu GNU/Linux*:

```console
sudo apt install libmetis-dev
```

However, there may be situations where downloading the packages may not be possible. In this case, other (potentially trickier) option is to download the source code and build the libraries. For more detailed and updated instructions for compiling *METIS*, *GKlib*, and other necessary packages, refer to [Compiling Kratos with MPI support](https://github.com/KratosMultiphysics/Kratos/wiki/Compiling-Kratos-with-MPI-support) in the wiki.

Assuming that the dependencies are installed, the following steps are:

1. Add the `MetisApplication` to the list of applications to be compiled in the building script for Kratos,
as described in the [install instructions](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md#adding-applications).

```console
export KRATOS_APPLICATIONS=
...
add_app ${KRATOS_APP_DIR}/MetisApplication
```

2. Tell CMake where are located the libraries and headers of *METIS*, if CMake does not find them automatically. For example:

```console
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}" \
...
-DMETIS_ROOT_DIR="${HOME}/Projects/METIS/build"
```

### Spack

*Spack* is a multi-platform package manager that builds and installs multiple versions and configurations of software. It works on *GNU/Linux*, *macOS*, and many supercomputers. *Spack* is non-destructive: installing a new version of a package does not break existing installations, so many configurations of the same package can coexist.

To install *Spack* you just need to run the following command:

```console
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
```

To use it you will need to add the corresponding environment variables (or call `spack_location/spack/bin/spack`):

```console
. spack_location/spack/share/spack/setup-env.sh
```

Then in order to install *METIS*:

```console
spack install parmetis
```

If you want that the libraries are added automatically to `LD_LIBRARY_PATH` to run the following commands before loading the modules:

```console
spack config add modules:prefix_inspections:lib64:[LD_LIBRARY_PATH]
spack config add modules:prefix_inspections:lib:[LD_LIBRARY_PATH]
```

Now you just need to load *METIS*:

```console
spack load parmetis
```

Once to compile `MetisApplication` just remember to add the application to the `cofigure` bash script:

```console
export KRATOS_APPLICATIONS=
...
add_app ${KRATOS_APP_DIR}/MetisApplication
```
