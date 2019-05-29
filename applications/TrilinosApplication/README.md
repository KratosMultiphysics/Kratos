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
- TrilinosConvectionDiffusionStrategy

#### Convergence Criterias
- TrilinosDisplacementCriteria
- TrilinosUPCriteria
- TrilinosResidualCriteria
- TrilinosAndCriteria
- TrilinosOrCriteria

#### Solving Strategies
- TrilinosSolvingStrategy
- TrilinosFSStrategy
- TrilinosLinearStrategy
- TrilinosNewtonRaphsonStrategy
- TrilinosLaplacianMeshMovingStrategy
- TrilinosStructuralMeshMovingStrategy

For more information about these please refer to their serial version (without trilinos prefix) in the main Kratos documentation.

### Components reference
* [__Zoltan__](https://trilinos.org/packages/zoltan).
* [__Epetra__](https://trilinos.org/packages/epetra).
* [__Aztec__](https://trilinos.org/packages/aztec).
* [__Amesos__](https://trilinos.org/packages/amesos).
* [__Teuchos__](https://trilinos.org/packages/teuchos).
