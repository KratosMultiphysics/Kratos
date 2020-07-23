## Trilinos Application

### Distributed matrices and linear solvers

Kratos Multiphysics uses several components of the [Trilinos project](https://trilinos.org/) for its MPI capabilities. The most relevant of these are the __distributed-memory matrix and vector classes__ from the Epetra module and the __linear solvers__ provided by the AztecOO, Amesos and ML packages. In addition, it also contains the MPI version of the AMGCL solver.

#### Main solver classes
- TrilinosLinearSolver
- AztecSolver
- AmesosSolver
- MultiLevelSolver
- AMGCL Mpi-based Solver

### Kratos interface extension

The Trilinos application also provides __MPI versions of most of the core classes of Kratos__, adapted to work with Epetra distributed matrices where necessary. Hence it provide its own version of the following Kratos (and application) interface elements:

#### Builder and solvers
- TrilinosResidualBasedBuilderAndSolver
- TrilinosEliminationBuilderAndSolver
- TrilinosBlockBuilderAndSolver
- TrilinosBlockBuilderAndSolverPeriodic

#### Convergence Criterias
- TrilinosDisplacementCriteria
- TrilinosUPCriteria
- TrilinosResidualCriteria
- TrilinosAndCriteria
- TrilinosOrCriteria

#### Solving Strategies
- TrilinosSolvingStrategy
- TrilinosLinearStrategy
- TrilinosNewtonRaphsonStrategy
- TrilinosLaplacianMeshMovingStrategy
- TrilinosStructuralMeshMovingStrategy

For more information about these please refer to their serial version (without trilinos prefix) in the main Kratos documentation.

### Components reference
* [__Epetra__](https://trilinos.org/packages/epetra).
* [__Aztec__](https://trilinos.org/packages/aztec).
* [__Amesos__](https://trilinos.org/packages/amesos).
* [__Teuchos__](https://trilinos.org/packages/teuchos).

### Notes for compilation
Trilinos is a large project and not all of its packages are being used in Kratos. Check the [configuration file of the CI (Travis)](https://github.com/KratosMultiphysics/Kratos/blob/master/.travis.yml) to see which packages are necessary in order to compile the TrilinosApplication.
Furthermore it is possible to do a minimal installation of the TrilinosApplication with only using the Epetra package. Using the other packages is enabled by default, but can be disabled with the following flags:
- *TRILINOS_EXCLUDE_ML_SOLVER*: Setting this flag to `ON` in the configure file will exclude the interface to the Trilinos ML solver package
- *TRILINOS_EXCLUDE_AZTEC_SOLVER*: Setting this flag to `ON` in the configure file will exclude solvers from the Trilinos AztecOO package
- *TRILINOS_EXCLUDE_AMESOS_SOLVER*: Setting this flag to `ON` in the configure file will exclude solvers using features of the Trilinos Amesos package
