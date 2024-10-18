# CoSimulation Application

The CoSimulation Application contains the core developments in coupling black-box solvers and other software-tools within Kratos Multiphysics.

<a name="overview"></a>


## Overview
- [CoSimulation Application](#cosimulation-application)
  - [Overview](#overview)
  - [List of features](#list-of-features)
  - [Dependencies](#dependencies)
  - [Examples](#examples)
  - [User Guide](#user-guide)
    - [Setting up a coupled simulation](#setting-up-a-coupled-simulation)
    - [The JSON configuration file](#the-json-configuration-file)
    - [Basic FSI example](#basic-fsi-example)
  - [Developer Guide](#developer-guide)
    - [Structure of the Application](#structure-of-the-application)
    - [How to couple a new solver / software-tool?](#how-to-couple-a-new-solver--software-tool)
      - [Interface of SolverWrapper](#interface-of-solverwrapper)
      - [Remote controlled CoSimulation](#remote-controlled-cosimulation)
    - [Using a solver in MPI](#using-a-solver-in-mpi)
  - [References](#references)

<a name="list-of-features"></a>


## List of features

- Various features available for CoSimulation:
  - [Coupling Algorithms](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/coupled_solvers)
  - [Wrappers for various solvers and other software-tools](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/solver_wrappers)
  - [Data Transfer Operators (including Mapping)](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/data_transfer_operators)
  - [Convergence Accelerators](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/convergence_accelerators)
  - [Convergence Criteria](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/convergence_criteria)
  - [Predictors](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/predictors)

- Support for MPI parallelization. This is independent of whether or not the ued solvers support/run in MPI.

- Coupling of Kratos <=> Kratos without overhead since the same database is used and data duplication is avoided.

- The [MappingApplication](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MappingApplication) is used for mapping between nonmatching grids.


<a name="dependencies"></a>


## Dependencies

The CoSimulation Application itself doesn't have any dependencies (except the `KratosCore` / `KratosMPICore` for serial/MPI-compilation).

For running coupled simulations the solvers to be used have to be available. Those dependencies are python-only.

The [MappingApplication](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MappingApplication) is required when mapping is used.


<a name="examples"></a>


## Examples

The examples can be found in the [examples repository](https://github.com/KratosMultiphysics/Examples/tree/master/co_simulation).
Please also refer to the [tests](tests) for examples of how the coupling can be configured.
Especially the [Mok-FSI](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/tests/fsi_mok) and the [Wall-FSI](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/tests/fsi_wall) tests are very suitable for getting a basic understanding.


<a name="user-guide"></a>


## User Guide

This section guides users of the _CoSimulationApplication_ to setting up and performing coupled simulations. The overall workflow is the same as what is used for most Kratos applications. It consists of the following files:

- **MainKratosCoSim.py** This file is to be executed with python to run the coupled simulation
- **ProjectParametersCoSim.json** This file contains the configuration for the coupled simulation


<a name="user-guide-setting-up-a-coupled-simulation"></a>


### Setting up a coupled simulation
For running a coulpled simulation at least the two files above are required. In addition, the input for the solvers / codes participating in the coupled simulation are necessary.

The **MainKratosCoSim.py** file looks like this (see also [here](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/MainKratosCoSim.py):

```py
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
```

It can be executed with python:

```
python MainKratosCoSim.py
```

If the coupled simulation runs in a distributed environment (MPI) then MPI is required to launch the script

```
mpiexec -np 4 python MainKratosCoSim.py --using-mpi
```

Not the passing of the `--using-mpi` flag which tells Kratos that it runs in MPI.


<a name="user-guide-the-json-configuration-file"></a>


### The JSON configuration file

The configuration of the coupled simulation is written in `json` format, same as for the rest of Kratos.

It contains two settings:
- _problem_data_: this setting contains global settings of the coupled problem.
  ```json
  "start_time" : 0.0,
  "end_time" : 15.0,
  "echo_level" : 0, // verbosity, higher values mean more output
  "print_colors" : true, // use colors in the prints
  "parallel_type" : "OpenMP" // or "MPI"
  ```

- _solver_settings_: the settings of the coupled solver.
  ```json
  "type" : "coupled_solvers.gauss_seidel_weak", // type of the coupled solver, see python_scripts/coupled_solvers
  "predictors" : [], // list of predictors
  "num_coupling_iterations" : 10, // max number of coupling iterations, only available for strongly coupled solvers
  "convergence_accelerators" : [] // list of convergence accelerators, only available for strongly coupled solvers
  "convergence_criteria" : [] // list of convergence criteria, only available for strongly coupled solvers
  "data_transfer_operators" : {} // map of data transfer operators (e.g. mapping)
  "coupling_sequence" : [] // list specifying in which order the solvers are called
  "solvers" : {} // map of solvers participating in the coupled simulation, specifying their input and interfaces
  ```

See the next section for a basic example with more explanations.


<a name="user-guide-basic-fsi-example"></a>


### Basic FSI example

This example is the Wall FSI benchmark, see [1], chapter 7.5.3. The Kratos solvers are used to solve this problem. The input files for this example can be found [here](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/tests/fsi_wall)

```json
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
```


<a name="developer-guide"></a>


## Developer Guide


<a name="developer-guide_structure-of-the-application"></a>


### Structure of the Application

The _CoSimulationApplication_ consists of the following main components (taken from [2]):
- **SolverWrapper**: Baseclass and CoSimulationApplication-interface for all solvers/codes participating in the coupled simulation, each solver/code has its own specific version.
- **CoupledSolver**: Implements coupling schemes such as weak/strong coupling with *Gauss-Seidel/Jacobi* pattern. It derives from SolverWrapper such that it can beused in nested coupled simulations.
- **IO**: Responsible for communicating and data exchange with external solvers/codes
- **DataTransferOperator**: Transfers data from one discretization to another, e.g. by use of mapping techniques
- **CouplingOperation**: Tool for customizing coupled simulations
- **ConvergenceAccelerator**: Accelerating the solution in strongly coupled simulations by use of relaxation or Quasi-Newton techniques
- **ConvergenceCriteria**: Checks if convergence is achieved in a strongly coupled simulation.
- **Predictor**: Improves the convergence by using a prediction as initial guess for the coupled solution

The following UML diagram shows the relation between these components:

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Documentation/blob/master/Readme_files/CoSimulationApplication/CoSimulation_uml.png?raw=true" style="width: 300px;"/>
</p>

Besides the functionalities [listed above](#list-of-features), the modular design of the application makes it straightforward to add a new or customized version of e.g. a _ConvergenceAccelerator_. It is not necessary to have those custom python scripts inside the _CoSimulationApplication_, it is sufficient that they are in a directory that is included in the _PYTHONPATH_ (e.g. the working directory).


<a name="developer-guide_how-to-couple-a-new-solver--software-tool"></a>


### How to couple a new solver / software-tool?

The _CoSimulationApplication_ is very modular and designed to be extended to coupling of more solvers / software-tools. This requires basically two components on the Kratos side:

The interface between the _CoSimulationApplication_ and a solver is done with the [**SolverWrapper**](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/base_classes/co_simulation_solver_wrapper.py). This wrapper is specific to every solver and calls the solver-custom methods based on the input of CoSimulation.

The second component necessary is an [**IO**](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/base_classes/co_simulation_io.py). This component is used by the SolverWrapper and is responsible for the exchange of data (e.g. mesh, field-quantities, geometry etc) between the solver and the _CoSimulationApplication_.

In principle three different options are possible for exchanging data with CoSimulation:

- For very simple solvers IO can directly be done in python inside the SolverWrapper, which makes a separate IO superfluous (see e.g. a [python-only single degree of freedom solver](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/solver_wrappers/sdof))
- Using the [_CoSimIO_](https://github.com/KratosMultiphysics/CoSimIO). This which is the preferred way of exchanging data with the _CoSimulationApplication_. It is currently available for _C++_, _C_, and _Python_. The _CoSimIO_ is included as the [KratosCoSimIO](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/solver_wrappers/kratos_co_sim_io.py) and can be used directly. Its modular and Kratos-independent design as _detached interface_ allows for easy integration into other codes.
- Using a custom solution based on capabilities that are offered by the solver that is to be coupled.

The following picture shows the interaction of these components with the _CoSimulationApplication_ and the external solver:

<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Readme_files/CoSimulationApplication/detached_interface.png" style="width: 300px;"/>
</p>


<a name="developer-guide_interface-of-solverwrapper"></a>


#### Interface of SolverWrapper

The [**SolverWrapper**](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/base_classes/co_simulation_solver_wrapper.py) is the interface in the _CoSimulationApplication_ to all involved codes / solvers. It provides the following interface (adapted from [2]), which is also called in this order:

- **Initialize**: This function is called once at the beginning of the simulation, it e.g .reads the input files and prepares the internal data structures.
- The solution loop is split into the following six functions:
  - **AdvanceInTime**: Advancing in time and preparing the data structure for the next time step.
  - **InitializeSolutionStep**: Applying boundary conditions
  - **Predict**: Predicting the solution of this time step to accelerate the solution.\
  iterate until convergence in a strongly coupled solution:
    - **SolveSolutionStep**: Solving the problem for this time step. This is the only function that can be called multiple times in an iterative (strongly coupled) solution procedure.
  - **FinalizeSolutionStep**: Updating internals after solving this time step.
  - **OutputSolutionStep**: Writing output at the end of a time step
- **Finalize**: Finalizing and cleaning up after the simulation

Each of these functions can implement functionalities to communicate with the external solver, telling it what to do. However, this is often skipped if the data exchange is used for the synchronization of the solvers. This is often done in "classical" coupling tools. I.e. the code to couple internally duplicates the coupling sequence and synchronizes with the coupling tool through the data exchange.

An example of a _SolverWrapper_ coupled to an external solver using this approach can be found [here](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/solver_wrappers/external/external_solver_wrapper.py). Only the mesh exchange is done explicitly at the beginning of the simulation, the data exchange is done inside _SolveSolutionStep_.

The coupled solver has to duplicate the coupling sequence, it would look e.g. like this (using _CoSimIO_) for a weak coupling:
```py
# solver initializes ...

CoSimIO::ExportMesh(...) # send meshes to the CoSimulationApplication

# start solution loop
while time < end_time:
    CoSimIO::ImportData(...) # get interface data

    # solve the time step

    CoSimIO::ExportData(...) # send new data to the CoSimulationApplication
```

An example for an FSI problem where the structural solver of Kratos is used as an external solver can be found [here](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/tests/structural_mechanics_analysis_with_co_sim_io.py). The _CoSimIO_ is used for communicating between the _CoSimulationApplication_ and the structural solver

While this approach is commonly used, it has the significant drawback that the coupling sequence has to be duplicated, which not only has the potential for bugs and deadlocks but also severely limits the useability when it comes to trying different coupling algorithms. Then not only the input for the _CoSimulationApplication_ has to be changed but also the source code in the external solver!

Hence a better solution is proposed in the next section:


<a name="developer-guide_remote-controller-cosimulation"></a>


#### Remote controlled CoSimulation
A unique feature of Kratos CoSimulation (in combination with the _CoSimIO_) is the remotely controlled CoSimulation. The main difference to the "classical" approach which duplicates the coupling sequence in the external solver is to give the full control to CoSimulation. This is the most flexible approach from the point of CoSimulation, as then neither the coupling sequence nor any other coupling logic has to be duplicated in the external solver.

In this approach the external solver registers the functions necessary to perform coupled simulations through the _CoSimIO_. These are then called remotely through the _CoSimulationApplication_. This way any coupling algorithm can be used without changing anything in the external solver.

```py
# defining functions to be registered
def SolveSolution()
{
    # external solver solves timestep
}

def ExportData()
{
    # external solver exports data to the CoSimulationApplication
}

# after defining the functions they can be registered in the CoSimIO:

CoSimIO::Register(SolveSolution)
CoSimIO::Register(ExportData)
# ...

# After all the functions are registered and the solver is fully initialized for CoSimulation, the Run method is called
CoSimIO::Run() # this function runs the coupled simulation. It returns only after finishing
```

A [simple example of this can be found in the _CoSimIO_](https://github.com/KratosMultiphysics/CoSimIO/blob/master/tests/integration_tutorials/cpp/run.cpp).

The _SolverWrapper_ for this approach sends a small control signal in each of its functions to the external solver to tell it what to do. This could be implemented as the following:

```py
class RemoteControlSolverWrapper(CoSimulationSolverWrapper):
    # ...
    # implement other methods as necessary
    # ...

    def InitializeSolutionStep(self):
        data_config = {
            "type"           : "control_signal",
            "control_signal" : "InitializeSolutionStep"
        }
        self.ExportData(data_config)

    def SolveSolutionStep(self):
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            # first tell the controlled solver to import data
            data_config = {
                "type"            : "control_signal",
                "control_signal"  : "ImportData",
                "data_identifier" : data_name
            }
            self.ExportData(data_config)

            # then export the data from Kratos
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

        # now the external solver solves
        super().SolveSolutionStep()

        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            # first tell the controlled solver to export data
            data_config = {
                "type"            : "control_signal",
                "control_signal"  : "ExportData",
                "data_identifier" : data_name
            }
            self.ExportData(data_config)

            # then import the data to Kratos
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)
```

A full example for this can be found [here](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/tests/structural_mechanics_analysis_remote_controlled.py).

If it is possible for an external solver to implement this approach, it is recommended to use it as it is the most robust and flexible.

Nevertheless both approaches are possible with the _CoSimulationApplication_.


<a name="developer-guide_using-a-solver-in-mpi"></a>


### Using a solver in MPI
By default, each _SolverWrapper_ makes use of all ranks in MPI. This can be changed if e.g. the solver that is wrapped by the _SolverWrapper_ does not support MPI or to specify to use less rank.

The base _SolverWrapper_ provides the `_GetDataCommunicator` function for this purpose. In the baseclass, the default _DataCommunicator_ (which contains all ranks in MPI) is returned. The _SolverWrapper_ will be instantiated on all the ranks on which this _DataCommunicator_ is defined (i.e. on the ranks where `data_communicator.IsDefinedOnThisRank() == True`).

If a solver does not support MPI-parallelism then it can only run on one rank. In such cases it should return a _DataCommunicator_ which contains only one rank. For this purpose the function `KratosMultiphysics.CoSimulationApplication.utilities.data_communicator_utilities.GetRankZeroDataCommunicator` can be used. Other custom solutions are also possible, see for example the [structural_solver_wrapper](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/python_scripts/solver_wrappers/kratos/structural_mechanics_wrapper.py).

<a name="references"></a>

## References

- [1] Wall, Wolfgang A., _Fluid structure interaction with stabilized finite elements_, PhD Thesis, University of Stuttgart, 1999, http://dx.doi.org/10.18419/opus-127
- [2] Bucher et al., _Realizing CoSimulation in and with a multiphysics framework_, conference proceedings, IX International Conference on Computational Methods for Coupled Problems in Science and Engineering, 2021, https://www.scipedia.com/public/Bucher_et_al_2021a
