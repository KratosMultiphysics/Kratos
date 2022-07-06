# Metis Application

The Metis application provides an interface to the
[METIS](https://github.com/KarypisLab/METIS) graph partitioning library,
used within Kratos Multiphysics to generate mesh partitions for MPI runs.

The applicaiton is organized in a series of processes, each one of them providing an interface to a different partitioning algorithm.
The algorithm most commonly used for finite elements is `metis_divide_heterogeneous_input_process.h`,
which generates a partition of the mesh based on the nodal graph and can be used for meshes combining multiple element types.

## Build instructions

Building the `MetisApplication` requires METIS libraries and their dependencies already installed on the system.
The easiest way to get them is by the package manager of the Linux distribution.
For example, in Ubuntu:

```bash
sudo apt install libmetis-dev
```

However, there may be situations where downloading the packages may not be possible.
In this case, other (potentially trickier) option is to download the source code and build the libraries.
For more detailed and updated instructions for compiling METIS, GKlib, and other necessary packages,
refer to [Compiling Kratos with MPI support](https://github.com/KratosMultiphysics/Kratos/wiki/Compiling-Kratos-with-MPI-support) in the wiki.

Assuming that the dependencies are installed, the following steps are:

1.Add the `MetisApplication` to the list of applications to be compiled in the building script for Kratos,
as described in the [install instructions](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md#adding-applications).

```bash
export KRATOS_APPLICATIONS=
...
add_app ${KRATOS_APP_DIR}/MetisApplication
```

2.Tell CMake where are located the libraries and headers of METIS, if CMake does not find them automatically. For example:

```bash
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}" \
...
-DMETIS_ROOT_DIR="${HOME}/Projects/METIS/build"
```
