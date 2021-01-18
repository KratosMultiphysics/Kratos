## CoSimulation Application

The CoSimulation Application contains the core developments in coupling black-box solvers and other software-tools within Kratos Multiphysics.

## Main Features:

- Various features available for CoSimulation:
    - [Coupling Algorithms](python_scripts/coupled_solvers)
    - [Wrappers for various solvers and other software-tools](python_scripts/solver_wrappers)
    - [Convergence Accelerators](python_scripts/convergence_accelerators)
    - [Convergence Criteria](python_scripts/convergence_criteria)
    - [Predictors](python_scripts/predictors)

- Besides the above mentioned functionalities the modular design of the application makes it straight forward to add a new or customized version of e.g. a ConvergenceAccelerator. It is not necessary to have those custom python scripts inside the CoSimulation, it is sufficient that they are in a directory that is included in the _PYTHONPATH_.

- Support for MPI parallelization. This is independent of whether or not the ued solvers support/run in MPI.

- Coupling of Kratos <=> Kratos without overhead since the same database is used and data duplication is avoided.

- The [MappingApplication](../MappingApplication) is used for mapping between nonmatching grids.

## Examples
This section is currently under construction.
Please refer to the [tests](tests) for examples of how the coupling can be configured.
Especially the [Mok-FSI](tests/fsi_mok) and the [Wall-FSI](tests/fsi_wall) tests are very suitable for getting a basic understanding.

## How to couple a new (external) solver / software-tool?
The CoSimulation Application is very modular and designed to be extended to coupling of more solvers / software-tools.This requires basically two components:

The interface between the CoSimulation and a solver is the done with the [**SolverWrapper**](python_scripts/base_classes/co_simulation_solver_wrapper.py). This wrapper is specific to every solver and calls the solver-custom methods based on the input of CoSimulation.

The second component necessary is an [**IO**](python_scripts/base_classes/co_simulation_io.py). This component is used by the SolverWrapper and is responsible for the exchange of data (e.g. mesh, field-quantities, geomety etc) between the solver and CoSimulation.

In principle three different options are possible for exchanging data with CoSimulation:

- For very simple solvers IO can directly be done in python inside the SolverWrapper, which makes a separate IO superfluous (see e.g. a [python-only single degree of freedom solver](python_scripts/solver_wrappers/sdof))
- Using the [_CoSimIO_](https://github.com/KratosMultiphysics/CoSimIO). This which is the preferred way of exchanging data with CoSimulation. It is currently available for _C++_, _C_, _Fortran_ and _Python_.
- Using a custom solution based on capabilities that are offered by the solver that is to be coupled.
