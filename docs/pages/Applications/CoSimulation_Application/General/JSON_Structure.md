---
title: Set-up Example
keywords: json
tags: [JSON_Structure.md]
sidebar: cosimulation_application
summary:
---
# Overview
This page will guide the users how to configure the json for the CoSimulation.

## Problem Data
The json file starts with the problem data. All the inputs in the problem data has to be identical as the other solvers used for the CoSimulation. For an example, the problem data of the [FSI Mok benchmark](../Examples/remote_CoSim_Mok_FSI.html)  between the json files of the Fluid Dynamics Application (ProjectParametersCFD.json), the Structural Mechanics Application (ProjectParametersCSM.json) and the CoSimulation Application (ProjectParametersCoSim.json) should be identical.
```json
{
    "problem_data" :
    {
        "start_time" : 0.0,
        "end_time" : 25.0,
        "echo_level" : 0,
        "print_colors" : true,
        "parallel_type" : "OpenMP"
    },
```
## Solver Settings

### Coupling strategy
Users can choose the coupling strategy at the "type" field. Certain strategy gives additional settings, such as the maximum number of iteration for each time step. These addtional settings can be found at the respective page of the coupling strategy.
```json
    "solver_settings" :
    {
        "type" : "<coupling_strategy>",
        "echo_level" : 1,
        "num_coupling_iterations" : 12,

```
The available names for ```<coupling_strategy>``` are listed below:
- [coupled_solvers.gauss_seidel_strong](../Coupling_Strategy/Gauss_Seidel_Strong.html)
- [coupled_solvers.gauss_seidel_weak](../Coupling_Strategy/Gauss_Seidel_Weak.html)
- [coupled_solvers.jacobi_strong](../Coupling_Strategy/Jacobi_Strong.html)
- [coupled_solvers.jacobi_weak](../Coupling_Strategy/Jacobi_Weak.html)
- [coupled_solvers.feti_dynamic_coupled_solver](../Coupling_Strategy/Feti_Dynamic_Coupled_Solver.html)

### Coupling Sequence
This part of the settings dictate the communication between the solvers used in the CoSimulation. The order of the solver executed will follow the order of the items inside the ```"coupling_sequence"``` field. Each item correspond to one solver. It has the ```name``` of the solver, the ```input_data_list``` and the ```output_data_list```.
The ```name``` field correspond on a solver naming that the user created in the ```"solvers"``` part of the solver settings. More details about how the solver is defined is explained [here](#solvers).
The data to be transfered are managed in the ```input_data_list``` and ```output_data_list```. The ```input_data_list``` contains the data to be used for the current solver. On the other hand, the ```output_data_list``` contains the data from this solver to be sent to another solver.

Each item inside ```"input_data_list"``` contain:
- ```"data"```                      : variable of the input data.
- ```"from_solver_data"```          : source for the input data.
- ```"from_solver" ```              : the solver of which source comes from.
- ```"data_transfer_operator"```    : the transfer operator for the data.

Similarly, ```"output_data_list"``` contain:
- ```"data"```                    : variable of the output data.
- ```"to_solver_data"```          : destination of the output data.
- ```"to_solver" ```              : the solver of the output destination.
- ```"data_transfer_operator"```  : the transfer operator for the data.

Optionally, the ```"input_data_list"``` and ```"output_data_list"``` may also include:
- ```"data_transfer_operator_options"```    : additional option that the user can use from the ```"data_transfer_operator"```
- ```"before_data_transfer_operations"```   : operation to be executed before the data transfer.
- ```"after_data_transfer_operations"```    : operation to be executed after the data transfer.


#### Coupling Sequence Example
Here's an example. The solvers used in this simulation are ```"fluid_solver"``` and ```"structure_solver"```. The Cosimulation will first execute the fluid solver, and then the structural solver. In this [first example](#coupling_sequence-setup) the details of the data to be transfered are managed solely inside ```"structure_solver"```. With this set-up, the structural solver will first gather the data listed in the ```"input_data_list"```. With the data naming convention made in the [solvers description](#solvers-example), the sequence are:
1. Run ```"fluid_solver"``` simulation.
2. Execute ```"Coupling_Operation_A"```.
3. Find a data called ```"fluid_solver_load"``` from ```"fluid_solver"```.
4. Execute ```"Data_Transfer_Operator_A"``` with its optional settings, namely ```"swap_sign"``` and ```"use_transpose"``` . Transfer the data into ```"structure_solver_load"```.
6. Execute ```"Coupling_Operation_B"```.
7. Run ```"structure_solver"``` simulation.
8. Copy ```"structure_solver_displacement"``` from the ```"structure_solver"```.
9. Execute ```"Data_Transfer_Operator_A"```. Transfer the data into ```"fluid_solver_displacement"``` of ```"fluid_solver"```.


**[here](#coupling-operations) are the details for the naming convention for the solvers used in this setup.*

**[here](#data-transfer-operator-example) are the details of the data transfer operator used in this setup.*

**[here](#coupling-operations) are the details of the coupling operation used in this setup.*

Alternatively, the coupling sequence could be set-up like [this](#coupling_sequence-alternative-1) or [this](#coupling_sequence-alternative-2). Both of them has the same meaning as the [first example](#coupling_sequence-example), but defined differently.

##### "coupling_sequence" setup
```json
        "coupling_sequence":
        [
            {
                "name": "fluid_solver",
                "input_data_list"  : [],
                "output_data_list" : []
            },
            {
                "name": "structure_solver",
                "input_data_list": [
                    {
                        "data"             : "structure_solver_load",
                        "from_solver"      : "fluid_solver",
                        "from_solver_data" : "fluid_solver_load",
                        "data_transfer_operator" : "Data_Transfer_Operator_A",
                        "data_transfer_operator_options" : ["swap_sign", "use_transpose"],
                        "before_data_transfer_operations" : ["Coupling_Operation_A"],
                        "after_data_transfer_operations" : ["Coupling_Operation_B"]
                    }
                ],
                "output_data_list": [
                    {
                        "data"           : "structure_solver_displacement",
                        "to_solver"      : "fluid_solver",
                        "to_solver_data" : "fluid_solver_displacement",
                        "data_transfer_operator" : "Data_Transfer_Operator_A"
                    }
                ]
            }
        ],
```
##### "coupling_sequence" alternative 1
```json
        "coupling_sequence":
        [
            {
                "name": "fluid_solver",
                "input_data_list"  : [],
                "output_data_list" : [
                    {
                        "data"           : "fluid_solver_load",
                        "to_solver"      : "structure_solver",
                        "to_solver_data" : "structure_solver_load",
                        "data_transfer_operator" : "Data_Transfer_Operator_A",
                        "data_transfer_operator_options" : ["swap_sign", "use_transpose"],
                        "before_data_transfer_operations" : ["Coupling_Operation_A"],
                        "after_data_transfer_operations" : ["Coupling_Operation_B"]
                    }
                ]
            },
            {
                "name": "structure_solver",
                "input_data_list": [],
                "output_data_list": [
                    {
                        "data"           : "structure_solver_displacement",
                        "to_solver"      : "fluid_solver",
                        "to_solver_data" : "fluid_solver_displacement",
                        "data_transfer_operator" : "Data_Transfer_Operator_A"
                    }
                ]
            }
        ],
```
##### "coupling_sequence" alternative 2
```json
        "coupling_sequence":
        [
            {
                "name": "fluid_solver",
                "input_data_list"  : [
                    {
                        "data"           : "fluid_solver_displacement",
                        "from_solver"      : "structure_solver",
                        "from_solver_data" : "structure_solver_displacement",
                        "data_transfer_operator" : "Data_Transfer_Operator_A"
                    }
                ],
                "output_data_list" : [
                    {
                        "data"           : "fluid_solver_load",
                        "to_solver"      : "structure_solver",
                        "to_solver_data" : "structure_solver_load",
                        "data_transfer_operator" : "Data_Transfer_Operator_A",
                        "data_transfer_operator_options" : ["swap_sign", "use_transpose"],
                        "before_data_transfer_operations" : ["Coupling_Operation_A"],
                        "after_data_transfer_operations" : ["Coupling_Operation_B"]
                    }
                ]
            },
            {
                "name": "structure_solver",
                "input_data_list": [],
                "output_data_list": []
            }
        ],
```
### Solvers
In this part, the user define the solvers to be used for the CoSimulation. For convenience, the user is free to assign a name for the solvers. For an example [here](#solvers-example), the fluid dynamics solver is named as ```"fluid_solver"```, and the structural mechanics solver is named as ```"structure_solver"```. The assigned names will only be used inside the json file of CoSimulation. After assigning a name, the user must declare the details of the solver. Here are the required fields:
- ```"type"```  : is the solver wrapper for the application or the solver. Here are the available solver wrappers:
    - [solver_wrappers.kratos.convection_diffusion_wrapper](../Coupling_Operations/Convection-Diffusion.html)
    - [solver_wrappers.kratos.dem_wrapper](../Coupling_Operations/DEM.html)
    - [solver_wrappers.kratos.fluid_dynamics_wrapper](../Coupling_Operations/Fluid_Dynamics.html)
    - [solver_wrappers.kratos.mpm_dirichlet_wrapper](../Coupling_Operations/MPM_Dirichlet.html)
    - [solver_wrappers.kratos.mpm_neumann_wrapper](../Coupling_Operations/MPM_Neumann.html)
    - [solver_wrappers.kratos.pfem_fluid_dynamics_wrapper](../Coupling_Operations/PFEM_Fluid_Dynamics.html)
    - [solver_wrappers.kratos.potential_flow_wrapper](../Coupling_Operations/Potential_Flow.html)
    - [solver_wrappers.kratos.structural_mechanics_wrapper](../Coupling_Operations/Structural_Mechanics.html)
    - [solver_wrappers.rigid_body.rigid_body_solver_wrapper](../Coupling_Operations/Rigid_Body_Solver.html)
    - [solver_wrappers.sdof.sdof_solver_wrapper](../Coupling_Operations/SDOF_Solver.html)
- ```"solver_wrapper_settings"```: contain the settings for the solver wrapper, such as the name of the ```"input_file"``` for the solver's project parameters.
- ```"data"```: contain the variables of the solver. To declare a variable, the user assign a name for it and pinpoint the variable source.
    - ```"model_part_name"```: the model part of which the variable is located.
    - ```"dimension"```: the dimention of the variable.
    - ```"variable_name"```: Kratos variable name of the solver. (not the assigned name by the user)

#### "solvers" example
Here is an example used for the coupling sequence described [here](#coupling-sequence-example)
```json
        "solvers" :
        {
            "fluid_solver":
            {
                "type" : "solver_wrappers.kratos.fluid_dynamics_wrapper",
                "solver_wrapper_settings" : {
                    "input_file"  : "ProjectParametersCFD"
                },
                "data" : {
                    "fluid_solver_displacement" : {
                        "model_part_name"   : "FluidModelPart.interface",
                        "dimension" : 2,
                        "variable_name" : "MESH_DISPLACEMENT"
                    },
                    "fluid_solver_load" : {
                        "model_part_name"   : "FluidModelPart.interface",
                        "dimension" : 2,
                        "variable_name" : "REACTION"
                    },
                    "fluid_solver_velocity" : {
                        "model_part_name"   : "FluidModelPart.interface",
                        "dimension" : 2,
                        "variable_name" : "VELOCITY"
                    }
                }
            },
            "structure_solver" :
            {
                "type" : "solver_wrappers.kratos.structural_mechanics_wrapper",
                "solver_wrapper_settings" : {
                    "input_file"  : "ProjectParametersCSM"
                },
                "data" : {
                    "structure_solver_displacement" : {
                        "model_part_name"   : "Structure.interface",
                        "dimension" : 2,
                        "variable_name" : "DISPLACEMENT"
                    },
                    "structure_solver_load" : {
                        "model_part_name"   : "Structure.interface",
                        "dimension" : 2,
                        "variable_name" : "POINT_LOAD"
                    }
                }
            }
        }
```
### Convergence Criteria
The user choose the criteria for convergence here. After the type is chosen, the variable for the convergence criteria is set by selecting the solver in the "solver" field and the variable from the solver in the "data_name" field. In the example below, the displacement of the fluid solver (```"fluid_solver_displacement"```) is used for determining the convergence.
```json
        "convergence_criteria" : [
            {
                "type"          : "convergence_type",
                "solver"        : "fluid_solver",
                "data_name"     : "fluid_solver_displacement",
                "abs_tolerance" : 1e-6,
                "rel_tolerance" : 1e-6
            }
        ],
```
The available ```"type"``` are listed below:
- [absolute_norm_energy_conjugate](../Convergence_Criteria/absolute_norm_energy_conjugate.html)
- [relative_norm_initial_residual](../Convergence_Criteria/relative_norm_initial_residual.html)
- [relative_norm_previous_residual](../Convergence_Criteria/relative_norm_previous_residual.html)

### Convergence Accelerators
The convergence of the coupling strategy such as the strong coupling can be improved by using convergence accelerator.
```json
        "convergence_accelerators" : [
            {
                "type"      : "mvqn",
                "solver"    : "fluid",
                "data_name" : "disp"
            }
        ],
```
The available convergence accelerators are listed below:
- [aitken](../Convergence_Accelerators/Aitken.html)
- [anderson](../Convergence_Accelerators/Anderson.html)
- [constant_relaxation](../Convergence_Accelerators/Constant_Relaxation.html)
- [iqnils](../Convergence_Accelerators/IQNILS.html)
- [mvqn](../Convergence_Accelerators/MVQN.html)

### Predictors
Predictors helps the solver by predicting the solution for the next timestep . It can be used both in a weak coupling and a strong coupling.
The available predictors are listed below:
- [average_value_based](../Predictors/Average_Value_Based.html)
- [linear](../Predictors/Linear.html)
- [linear_derivative_based](../Predictors/Linear_Derivative_Based.html)

```json
    "predictors" : [
        {
            "type" : "average_value_based",
            "solver"         : "fluid",
            "data_name"      : "reaction"
        }
    ],
```


### Data Transfer Operators
Transfering data between two simulation with the same solver may be simple, but it is a different case for data transfer between two different solver. Take an FSI simulation for an example. The fluid solver would need a fine mesh, especially at the interface between the fluid domain and the structural domain. Meanwhile, structural solver don't need a fine mesh to solve the structure. This leads to a complication at the interface, where both solver has their own mesh. To solve this, the user can transform the data before sending it to another solver.

The available data transfer operators are listed below:
- [copy](../Data_Transfer_Operators/copy.html)
- [copy_single_to_distributed_vectorial](../Data_Transfer_Operators/copy_single_to_distributed_vectorial.html)
- [copy_with_empty_ranks](../Data_Transfer_Operators/copy_with_empty_ranks.html)
- [kratos_mapping](../Data_Transfer_Operators/kratos_mapping.html)
- [sum_distributed_to_single](../Data_Transfer_Operators/sum_distributed_to_single.html)
- [sum_distributed_to_single_vectorial](../Data_Transfer_Operators/sum_distributed_to_single_vectorial.html)
- [sum_many_to_many](../Data_Transfer_Operators/sum_many_to_many.html)
- [transfer_one_to_many](../Data_Transfer_Operators/transfer_one_to_many.html)

#### Data transfer operator example
Here is an example used for the coupling sequence described [here](#coupling-sequence-example)
```json
        "data_transfer_operators" : {
            "Data_Transfer_Operator_A" : {
                "type" : "kratos_mapping",
                "mapper_settings" : {
                    "mapper_type" : "nearest_neighbor"
                }
            }
        },
```

### Coupling Operations
The operators that has been established here can be executed before or after transfering the data. Take a look at the [coupling sequence](#coupling-sequence) for more details.
```json
    "coupling_operations" : {
        "Coupling_Operation_A" : {
            "type"      : "compute_resultants",
            "solver"    : "fluid",
            "data_name" : "force"
        },
        "Coupling_Operation_B" : {
            "type"      : "impose_mesh_displacement",
            "solver"    : "fluid",
            "data_name" : "disp"
        }
    },
```
The available coupling operations are listed below:
- [compute_boundary_force](../Coupling_Operations/Compute_Boundary_Force.html)
- [compute_normals](../Coupling_Operations/Compute_Normals.html)
- [compute_resultants](../Coupling_Operations/Compute_Resultants.html)
- [convert_distributed_values_to_point](../Coupling_Operations/Convert_Distributed_Values_to_Point.html)
- [coupling_output](../Coupling_Operations/Coupling_Output.html)
- [create_point_load_model_part](../Coupling_Operations/Create_Point_Load_Model_Part.html)
- [distribute_point_values](../Coupling_Operations/Distribute_Point_Values.html)
- [elemental_data_to_nodal_data](../Coupling_Operations/Elemental_Data_to_Nodal_Data.html)
- [impose_mesh_displacement](../Coupling_Operations/Impose_Mesh_Displacement.html)
- [print_iteration_number](../Coupling_Operations/Print_Iteration_Number.html)
- [reset_pfem_kinematics](../Coupling_Operations/Reset_PFEM_Kinematics.html)
- [scaling](../Coupling_Operations/Scaling.html)

