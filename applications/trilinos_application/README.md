## Trilinos Application

This application contains the wrappers for different components of [Trilinos project](https://trilinos.org/). Kratos use of trilinos is oriented
parallel solvers and mainly uses two of trilinos components.

Trilinos application provide wrappers for the following functionalities 

#### Solvers
- TrilinosLinearSolver
- AztecSolver
- AmesosSolver
- MultiLevelSolver
- AMGCL Mpi-based Solver

### Kratos Interface Extension

Trilinos in Kratos extends the normal kratos interface to be used with trilinos backend. Hance it provide its own
version of common Kratos ( and application ) interface elements:

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
