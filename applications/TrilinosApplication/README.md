# Trilinos Application

## Distributed matrices and linear solvers

Kratos Multiphysics uses several components of the [Trilinos project](https://trilinos.org/) for its MPI capabilities. The most relevant of these are the __distributed-memory matrix and vector classes__ from the Epetra module and the __linear solvers__ provided by the *AztecOO*, *Amesos* and *ML* packages. In addition, it also contains the **MPI** version of the **AMGCL** solver.

### Main solver classes

- `TrilinosLinearSolver`
- `AztecSolver`
- `AmesosSolver`
- `Amesos2Solver`
- `MultiLevelSolver`
- `AMGCL Mpi-based Solver`

## Kratos interface extension

The Trilinos application also provides __MPI versions of most of the core classes of Kratos__, adapted to work with Epetra distributed matrices where necessary. Hence it provide its own version of the following Kratos classes:

### Builder and solvers

- `TrilinosResidualBasedBuilderAndSolver`
- `TrilinosEliminationBuilderAndSolver`
- `TrilinosBlockBuilderAndSolver`
- `TrilinosBlockBuilderAndSolverPeriodic`

### Convergence Criterias

- `TrilinosDisplacementCriteria`
- `TrilinosResidualCriteria`
- `TrilinosAndCriteria`
- `TrilinosOrCriteria`
- `TrilinosMixedGenericCriteria`

### Solving Strategies

- `TrilinosSolvingStrategy`
- `TrilinosLinearStrategy`
- `TrilinosNewtonRaphsonStrategy`

For more information about these please refer to their serial version (without _Trilinos_ prefix) in the main Kratos documentation.

## Components reference

- [__Epetra__](https://trilinos.github.io/epetra.html)
- [__Aztec__](https://trilinos.github.io/aztecoo.html)
- [__Amesos__](https://trilinos.github.io/amesos.html)
- [__Amesos2__](https://trilinos.github.io/amesos2.html)
- [__Teuchos__](https://trilinos.github.io/teuchos.html)

## Build instructions

Building the `TrilinosApplication` requires the *Trilinos* libraries and their dependencies already installed on the system.

### Direct GNU/Linux install

The easiest way to get them is by the package manager of the *GNU/Linux* distribution. For example, in *Ubuntu GNU/Linux*:

```Shell
sudo apt install trilinos-all-dev
```

However, there may be situations where downloading the packages may not be possible.

In this case, other (potentially trickier) option is to download the source code and build the libraries. For more detailed and updated instructions for compiling *Trilinos* and other necessary packages, refer to [Compiling Kratos with MPI support](https://github.com/KratosMultiphysics/Kratos/wiki/Compiling-Kratos-with-MPI-support), in the wiki.

Assuming that the dependencies are installed, the following steps are:

1. Add the `TrilinosApplication` to the list of applications to be compiled in the building script for Kratos,
as described in the [install instructions](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md#adding-applications).

```console
export KRATOS_APPLICATIONS=
...
add_app ${KRATOS_APP_DIR}/TrilinosApplication
```

2. Tell cmake where are located the libraries and headers of *Trilinos*.
If Trilinos was compiled (instead of downloaded with a package manager), it is usually enough to point the `TRILINOS_ROOT` variable to the build directory. For example:

```console
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}" \
...
-DTRILINOS_ROOT="${HOME}/Projects/Trilinos/build"
```

Or, if *Trilinos* is installed with a package manager, then libraries and headers may be in different locations.
Moreover, the name of the libraries may not be standard.
In this case, instead of setting `TRILINOS_ROOT`, set `-DTRILINOS_INCLUDE_DIR=String` with the path to the include dir, `-DTRILINOS_LIBRARY_DIR=String` with the path to the library dir, and set `-DTRILINOS_LIBRARY_PREFIX=String` with the prefix to use when looking for the *Trilinos* libraries, i.e.,

```console
libepetra.so  # No need to set TRILINOS_PREFIX
libtrilinos_epetra.so  # -DTRILINOS_PREFIX="trilinos_"
```

For example, in the case of Ubuntu with *Trilinos* installed by `sudo apt install trilinos-all-dev`:

```console
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos"
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu"
-DTRILINOS_LIBRARY_PREFIX="trilinos_"
```

_Notes_:

-*Trilinos* is a large project and not all of its packages are being used in *Kratos*.
Check the [docker of the CI](https://github.com/KratosMultiphysics/Kratos/blob/master/scripts/docker_files/docker_file_ci_ubuntu_20_04/DockerFile) to see which packages are necessary in order to compile the TrilinosApplication.

At the moment, the list of required packages is:

```console
sudo apt install \
        libtrilinos-amesos-dev \
        libtrilinos-amesos2-dev \
        libtrilinos-aztecoo-dev \
        libtrilinos-epetra-dev \
        libtrilinos-epetraext-dev \
        libtrilinos-ifpack-dev \
        libtrilinos-ml-dev \
        libtrilinos-teuchos-dev \
        libtrilinos-tpetra-dev \
        libtrilinos-kokkos-dev \
        libtrilinos-kokkos-kernels-dev \
        libtrilinos-shylu-dev
```

-It is possible to do a minimal installation of the TrilinosApplication only using the Epetra package.
Other packages can be disabled with the following flags:

```console
-DTRILINOS_EXCLUDE_ML_SOLVER=ON  # exclude the interface to the Trilinos ML solver package
-DTRILINOS_EXCLUDE_AZTEC_SOLVER=ON  # exclude solvers from the Trilinos AztecOO package
-DTRILINOS_EXCLUDE_AMESOS_SOLVER=ON  # exclude solvers using features of the Trilinos Amesos package
-DTRILINOS_EXCLUDE_AMESOS2_SOLVER=ON  # exclude solvers using features of the Trilinos Amesos2 package
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

Then in order to install *Trilinos*:

```console
spack install trilinos
```

In case we want to use [*MUMPS*](https://graal.ens-lyon.fr/MUMPS/index.php) and [SuperLUDist](https://portal.nersc.gov/project/sparse/superlu/) with *Trilinos*, you can install them all together and properly liked with:

```console
spack install trilinos+mumps+superlu-dist
```

If you want that the libraries are added automatically to `LD_LIBRARY_PATH` to run the following commands before loading the modules:

```console
spack config add modules:prefix_inspections:lib64:[LD_LIBRARY_PATH]
spack config add modules:prefix_inspections:lib:[LD_LIBRARY_PATH]
```

Now you just need to load *Trilinos*:

```console
spack load trilinos
```

or (if installed):

```console
spack load trilinos+mumps+superlu-dist
```

Once to compile `TrilinosApplication` just remember to add the application to the `configure` bash script:

```console
export KRATOS_APPLICATIONS=
...
add_app ${KRATOS_APP_DIR}/TrilinosApplication
```