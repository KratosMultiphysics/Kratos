# CoSimulation Application

The CoSimulation Application contains the core developments in coupling black-box solvers and other software-tools within Kratos Multiphysics.

## Overview

- [List of features](#list-of-features)
- [Dependencies](#dependencies)
- [Examples](#examples)
- [User Guide](#user-guide)
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

This section guides users of the _CoSimulationApplication_ to setting up and performing coupled simulations. The overall workflow is the same as what is used for most Kratos applications. It consists of the following files:

- **MainKratosCoSim.py** This file is to be executed with python to run the coupled simulation
- **ProjectParametersCoSim.json** This file contains the configuration for the coupled simulation

### Setting up a coupled simulation
For running a coulpled simulation at least the two files above are required. In addition, the input for the solvers / codes participating in the coupled simulation are necessary.

The **MainKratosCoSim.py** file looks like this (see also [here](python_scripts/MainKratosCoSim.py):

~~~py
import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

"""
For user-scripting it is intended that a new class is derived
from CoSimulationAnalysis to do modifications
Check also "kratos/python_scripts/analysis-stage.py" for available methods that can be overridden
"""

parameter_file_name = "ProjectParametersCoSim.json"
with open(parameter_file_name,'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

simulation = CoSimulationAnalysis(parameters)
simulation.Run()
~~~

It can be executed with python:

~~~
python MainKratosCoSim.py
~~~

If the coupled simulation runs in a distributed environment (MPI) then MPI is required to launch the script

~~~
mpiexec -np 4 python MainKratosCoSim.py --using-mpi
~~~~

Not the passing of the `--using-mpi` flag which tells Kratos that it runs in MPI.



### The json configuration file

The configuration of the coupled simulation is written in `json` format, same as for the rest of Kratos.

The following sections are important:
- _problem_data_: this setting contains global settings of the coupled problem.
  ~~~js
  "start_time" : 0.0,
  "end_time" : 15.0,
  "echo_level" : 0, // verbosity, higher values mean more output
  "print_colors" : true, // use colors in the prints
  "parallel_type" : "OpenMP" // or "MPI"
  ~~~

- _solver_settings_: the settings of the coupled solver.
  ~~~js
  "type" : "coupled_solvers.gauss_seidel_weak", // type of the coupled solver, see python_scripts/coupled_solvers
  "predictors" : [], // list of predictors
  "num_coupling_iterations" : 10, // max number of coupling iterations, only available for strongly coupled solvers
  "convergence_accelerators" : [] // list of convergence accelerators, only available for strongly coupled solvers
  "convergence_criteria" : [] // list of convergence criteria, only available for strongly coupled solvers
  "data_transfer_operators" : {} // map of data transfer operators (e.g. mapping)
  "coupling_sequence" : [] // list specifying in which order the solvers are called
  "solvers" : {} // map of solvers participating in the coupled simulation, specifying their input and interfaces
  ~~~~

See the next section for a basic example with more explanations.
### Basic FSI example

This example is the Wall FSI benchmark which coupled Kratos solvers:

~~~js
{
    "problem_data" :
    {
        "start_time" : 0.0,
        "end_time" : 3.0,
        "echo_level" : 0, // printing no additional output
        "print_colors" : true, // using colors for prints
        "parallel_type" : "OpenMP"
    },
    "solver_settings" :
    {
        "type" : "coupled_solvers.gauss_seidel_weak", // weakly coupled simulation, no interface convergence is checked
        "echo_level" : 0, // no additional output from the coupled solver
        "predictors" : [ // using a predictor to improve the stability of the simulation
            {
                "type" : "average_value_based",
                "solver"         : "fluid",
                "data_name"      : "load"
            }
        ],
        "data_transfer_operators" : {
            "mapper" : {
                "type" : "kratos_mapping",
                "mapper_settings" : {
                    "mapper_type" : "nearest_neighbor" // using a simple mapper, see the README in the MappingApplications
                }
            }
        },
        "coupling_sequence":
        [
        {
            "name": "structure", // the structural solver comes first
            "input_data_list": [ // before solving, the following data is imported in the structural solver
                {
                    "data"              : "load",
                    "from_solver"       : "fluid",
                    "from_solver_data"  : "load", // the fluid loads are mapped onto the structure
                    "data_transfer_operator" : "mapper", // using the mapper defined above (nearest neighbor)
                    "data_transfer_operator_options" : ["swap_sign"] // in Kratos, the loads have the opposite sign, hence it has to be swapped
                }
            ],
            "output_data_list": [ // after solving, the displacements are mapped to the fluid solver
                {
                    "data"           : "disp",
                    "to_solver"      : "fluid",
                    "to_solver_data" : "disp",
                    "data_transfer_operator" : "mapper"
                }
            ]
        },
        {
            "name": "fluid", // the fluid solver solves after the structure
            "output_data_list": [],
            "input_data_list": []
        }
        ],
        "solvers" : // here we specify the solvers, their input and interfaces for CoSimulation
        {
            "fluid":
            {
                "type" : "solver_wrappers.kratos.fluid_dynamics_wrapper", // using the Kratos FluidDynamicsApplication for the fluid
                "solver_wrapper_settings" : {
                    "input_file"  : "fsi_wall/ProjectParametersCFD" // input file for the fluid solver
                },
                "data" : { // definition of interfaces used in the simulation
                    "disp" : {
                        "model_part_name" : "FluidModelPart.NoSlip2D_FSI_Interface",
                        "variable_name" : "MESH_DISPLACEMENT",
                        "dimension" : 2
                    },
                    "load" : {
                        "model_part_name" : "FluidModelPart.NoSlip2D_FSI_Interface",
                        "variable_name" : "REACTION",
                        "dimension" : 2
                    }
                }
            },
            "structure" :
            {
                "type" : "solver_wrappers.kratos.structural_mechanics_wrapper", // using the Kratos StructuralMechanicsApplication for the structure
                "solver_wrapper_settings" : {
                    "input_file"  : "fsi_wall/ProjectParametersCSM" // input file for the structural solver
                },
                "data" : { // definition of interfaces used in the simulation
                    "disp" : {
                        "model_part_name" : "Structure.GENERIC_FSI_Interface",
                        "variable_name" : "DISPLACEMENT",
                        "dimension" : 2
                    },
                    "load" : {
                        "model_part_name" : "Structure.GENERIC_FSI_Interface",
                        "variable_name" : "POINT_LOAD",
                        "dimension" : 2
                    }
                }
            }
        }
    }
}
~~~


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