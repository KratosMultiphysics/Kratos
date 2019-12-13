## CoSimulation Application

The CoSimulation Application contains the core developments in coupling black-box solvers within Kratos Multiphysics.

### Features:

- Coupling (weak / strong)

- Various features available for CoSimulation
    - Convergence Accelerators
    - Convergence Criteria
    - Predictors

- Support for MPI parallelization (if invoved solvers support it)

- Coupling of Kratos <=> Kratos without overhead

### IO with external solvers
Different ways of IO with external solvers are possible
- directly in python (see SDof solvers)
- using the CoSimIO. This is the preferred way for C++, C or Fortran based solvers
- using the legacy EMPIRE-API
