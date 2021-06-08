## Trilinos Application

### Distributed matrices and linear solvers

Kratos Multiphysics uses several components of the [Trilinos project](https://trilinos.org/) for its MPI capabilities. The most relevant of these are the __distributed-memory matrix and vector classes__ from the Epetra module and the __linear solvers__ provided by the AztecOO, Amesos and ML packages. In addition, it also contains the MPI version of the AMGCL solver.

#### Main solver classes
- TrilinosLinearSolver
- AztecSolver
- AmesosSolver
- Amesos2Solver
- MultiLevelSolver
- AMGCL Mpi-based Solver

### Kratos interface extension

The Trilinos application also provides __MPI versions of most of the core classes of Kratos__, adapted to work with Epetra distributed matrices where necessary. Hence it provide its own version of the following Kratos classes:

#### Builder and solvers
- TrilinosResidualBasedBuilderAndSolver
- TrilinosEliminationBuilderAndSolver
- TrilinosBlockBuilderAndSolver
- TrilinosBlockBuilderAndSolverPeriodic

#### Convergence Criterias
- TrilinosDisplacementCriteria
- TrilinosResidualCriteria
- TrilinosAndCriteria
- TrilinosOrCriteria
- TrilinosMixedGenericCriteria

#### Solving Strategies
- TrilinosSolvingStrategy
- TrilinosLinearStrategy
- TrilinosNewtonRaphsonStrategy

For more information about these please refer to their serial version (without _Trilinos_ prefix) in the main Kratos documentation.

### Components reference
* [__Epetra__](https://trilinos.github.io/epetra.html)
* [__Aztec__](https://trilinos.github.io/aztecoo.html)
* [__Amesos__](https://trilinos.github.io/amesos.html)
* [__Amesos2__](https://trilinos.github.io/amesos2.html)
* [__Teuchos__](https://trilinos.github.io/teuchos.html)


### Build instructions
First add the `TrilinosApplication` to the compiled apps, as described in the [install instructions](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md#adding-applications).

Use the following settings to add Trilinos:

`-DTRILINOS_ROOT=String`

Root directory for Trilinos library.

`-DTRILINOS_INCLUDE_DIR=String`

Not required if `TRILINOS_ROOT` is set. Path to trilinos include dir.

`-DTRILINOS_LIBRARY_DIR=String`

Not required if `TRILINOS_ROOT` is set. Path to trilinos library dir.

`-DTRILINOS_LIBRARY_PREFIX=String`
Indicates the prefix of the trilinos libraries in case they have:
```
libepetra.so          -> No prefix
libtrilinos_epetra.so -> -DTRILINOS_PREFIX="trilinos_"
```
If trilinos was installed using the package manager usually the following lines have to be used:
```
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \
-DTRILINOS_LIBRARY_PREFIX="trilinos_" \
```

### Notes for compilation
Trilinos is a large project and not all of its packages are being used in Kratos. Check the [docker of the CI](https://github.com/KratosMultiphysics/Kratos/blob/master/scripts/docker_files/docker_file_ci_ubuntu_20_04/DockerFile) to see which packages are necessary in order to compile the TrilinosApplication.
Furthermore it is possible to do a minimal installation of the TrilinosApplication with only using the Epetra package. Using the other packages is enabled by default, but can be disabled with the following flags:
- *TRILINOS_EXCLUDE_ML_SOLVER*: Setting this flag to `ON` in the configure file will exclude the interface to the Trilinos ML solver package
- *TRILINOS_EXCLUDE_AZTEC_SOLVER*: Setting this flag to `ON` in the configure file will exclude solvers from the Trilinos AztecOO package
- *TRILINOS_EXCLUDE_AMESOS_SOLVER*: Setting this flag to `ON` in the configure file will exclude solvers using features of the Trilinos Amesos package
- *TRILINOS_EXCLUDE_AMESOS2_SOLVER*: Setting this flag to `ON` in the configure file will exclude solvers using features of the Trilinos Amesos2 package


From Ubuntu 18.04 onwards, Trilinos can be installed with the following command:

```Shell
sudo apt-get install trilinos-all-dev
```


