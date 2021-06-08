# CoSimulation Application

The CoSimulation Application contains the core developments in coupling black-box solvers and other software-tools within Kratos Multiphysics.

## Overview

- [List of features](#list-of-features)
- [Dependencies](#dependencies)
- [Examples](#examples)
- [User Guide](#user-guide)
  - [Running a coupled simulation](#running-a-coupled-simulation)
  - [Setting up a coupled simulation](#setting-up-a-coupled-simulation)
  - [The json configuration file](#the-json-configuration-file)
- [Developer Guide](#developer-guide)
  - [Structure of the Application](#structure-of-the-application)
  - [How to couple a new solver / software-tool?](#how-to-couple-a-new-solver--software-tool)
- [References](#references)

## List of features

- Various features available for CoSimulation:
  - [Coupling Algorithms](python_scripts/coupled_solvers)
  - [Wrappers for various solvers and other software-tools](python_scripts/solver_wrappers)
  - [Data Transfer Operators (including Mapping)](python_scripts/data_transfer_operators)
  - [Convergence Accelerators](python_scripts/convergence_accelerators)
  - [Convergence Criteria](python_scripts/convergence_criteria)
  - [Predictors](python_scripts/predictors)

- Besides the above mentioned functionalities the modular design of the application makes it straight forward to add a new or customized version of e.g. a ConvergenceAccelerator. It is not necessary to have those custom python scripts inside the _CoSimulationApplication_, it is sufficient that they are in a directory that is included in the _PYTHONPATH_ (e.g. the working directory).

- Support for MPI parallelization. This is independent of whether or not the ued solvers support/run in MPI.

- Coupling of Kratos <=> Kratos without overhead since the same database is used and data duplication is avoided.

- The [MappingApplication](../MappingApplication) is used for mapping between nonmatching grids.

## Dependencies

The CoSimulation Application itself doesn't have any dependencies (except the `KratosCore` / `KratosMPICore` for serial/MPI-compilation).

For running coupled simulations the solvers to be used have to be available. Those dependencies are python-only.

The [MappingApplication](../MappingApplication) is required when mapping is used.


## Examples

This section is currently under construction.
Please refer to the [tests](tests) for examples of how the coupling can be configured.
Especially the [Mok-FSI](tests/fsi_mok) and the [Wall-FSI](tests/fsi_wall) tests are very suitable for getting a basic understanding.


## User Guide

### Running a coupled simulation
### Setting up a coupled simulation
### The json configuration file


## Developer Guide



### Structure of the Application

The _CoSimulationApplication_ consists of the following main components (taken from [1]):
- **SolverWrapper**: Baseclass and CoSimulationApplication-interface for all solvers/codes participating in the coupled simulation, each solver/code has its own specific version.
- **CoupledSolver**: Implements coupling schemes such as weak/strong coupling withGauss-Seidel/Jacobi pattern.  It derives from SolverWrapper such that it can beused in nested coupled simulations.
- **IO**: Responsible for communicating and data exchange with external solvers/codes
- **DataTransferOperator**: Transfers data from one discretization to another, e.g.by use of mapping techniques
- **CouplingOperation**: Tool for customizing coupled simulations
- **ConvergenceAccelerator**: Accelerating the solution in strongly coupled simula-tions by use of relaxation techniques
- **ConvergenceCriteria**: Checks if convergence is achieved in a strongly coupledsimulation.
- **Predictor**: Improves the convergence by using a prediction as initial guess for thecoupled solution

The following UML diagram shows the relation between these components:

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Examples/blob/master/co_simulation/CoSimulation_uml.png" style="width: 400px;"/>
</p>


### How to couple a new solver / software-tool?

The CoSimulation Application is very modular and designed to be extended to coupling of more solvers / software-tools.This requires basically two components:

The interface between the CoSimulation and a solver is the done with the [**SolverWrapper**](python_scripts/base_classes/co_simulation_solver_wrapper.py). This wrapper is specific to every solver and calls the solver-custom methods based on the input of CoSimulation.

The second component necessary is an [**IO**](python_scripts/base_classes/co_simulation_io.py). This component is used by the SolverWrapper and is responsible for the exchange of data (e.g. mesh, field-quantities, geomety etc) between the solver and CoSimulation.

In principle three different options are possible for exchanging data with CoSimulation:

- For very simple solvers IO can directly be done in python inside the SolverWrapper, which makes a separate IO superfluous (see e.g. a [python-only single degree of freedom solver](python_scripts/solver_wrappers/sdof))
- Using the [_CoSimIO_](https://github.com/KratosMultiphysics/CoSimIO). This which is the preferred way of exchanging data with CoSimulation. It is currently available for _C++_, _C_, and _Python_. The _CoSimIO_ is included as the [KratosCoSimIO](python_scripts/solver_wrappers/kratos_co_sim_io.py) and can be used directly.
- Using a custom solution based on capabilities that are offered by the solver that is to be coupled.

#### Interface of SolverWrapper

The [**SolverWrapper**](python_scripts/base_classes/co_simulation_solver_wrapper.py) is the interface in the _CoSimulationApplication_ to all involved codes / solvers.




provides the following interface (adapted from [1]):

- **Initialize**:
- **Finalize**:
- **AdvanceInTime**:
- **InitializeSolutionStep**:
- **Predict**:
- **SolveSolutionStep**:
- **FinalizeSolutionStep**:
- **OutputSolutionStep**:


#### Remote controlled CoSimulation
A unique feature of Kratos CoSimulation (in combination with the _CoSimIO_) is the remotely controlled CoSimulation.


## References

- [1] Bucher et al., _Realizing CoSimulation in and with a multiphysics framework_, conference proceedings, IX International Conference on Computational Methods for Coupled Problems in Science and Engineering, 2021, under review