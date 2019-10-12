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
- includes are header-only (in simple cases even only one file)
- IO is done Asynchronous, hence less chances for deadlocking

- Currently available:
    - IO through the legacy-Interface of EMPIRE
        - Carat++
        - OpenFOAM++
        - FASTEST

- Future:
    - PythonIO
    - VTKFileIO
    - VTUFileIO
    - SocketIO
    - MPIIO