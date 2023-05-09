---
title: Rigid Body Solver
keywords: 
tags: [Rigid_Body_Solver.md]
sidebar: cosimulation_application
summary: 
---
# Overview
The Rigid Body solver can be seen as a combination of 6 single-degree-of-freedom. (SDOF) solvers, solving independently the three displacements and three rotations of a rigid body. It is meant to be used with CoSimulationApplication (e.g. rigid object submerged in a certain flow). For standalone usage, it can also be called directly, since the class acts at the same time as a solver and as an analysis.

**Note that the Rigid Body solver is not suitable for solving a rigid body movement with more than one active rotation.** It solves the displacements and rotations indentpendently, therefore, any coupling effect are not considered.

     ---------------        Sub model part: RigidBody                  \ 
    |               |       Node ID: 1                                 |
    |    Node 1     |       Specific variables:                        |
    |       X       |           FORCE, MOMENT,                         |
    |    (0,0,0)    |           PRESCRIBED_FORCE, PRESCRIBED_MOMENT    |
    |               |           EFFECTIVE_FORCE, EFFECTIVE_MOMENT      |    Model part: Main
     ---------------            BODY_FORCE, BODY_MOMENT                |    Nodes IDs: 1, 2
        |       |                                                       >   General variables:
       <_      _|_                                                     |    DISPLACEMENT, ROTATION,
        _>    |___|                                                    |    VELOCITY, ANGULAR_VELOCITY
       <_     | | |                                                    |    ACCELERATION, ANGULAR_ACCELERATION
         >      |           Sub model part: RootPoint                  |
        |       |           Node ID: 2                                 |
    /////// X ////////      Specific variables:                        |
         Node 2                 REACTION, REACTION_MOMENT              |
         (0,0,0)                PRESCRIBED_DISPLACEMENT,               /
                                PRESCRIBED_ROTATION
 

The solver consist of two sub model parts, with both of them residing in one node (0,0,0). The sub model parts and its variables are:
- "RigidBody" (representing the body itself)
    Sub model part specific variables:
    - FORCE
    - MOMENT
    - PRESCRIBED_FORCE
    - PRESCRIBED_MOMENT
    - EFFECTIVE_FORCE
    - EFFECTIVE_MOMENT
    - BODY_FORCE
    - BODY_MOMENT             
- "RootPoint" (representing the attachment point). 
    Sub model part specific variables:
    - REACTION
    - REACTION_MOMENT  
    - PRESCRIBED_DISPLACEMENT
    - PRESCRIBED_ROTATION
An ilustration of the model and its sub model parts can be seen on the sketch above.



It is possible to apply forces to the rigid body as well as displacements to the reference point (e.g. to be used as a TMD). The "PRESCRIBED_*" variables have the same behaviour as their original version but are necessary to avoid overwriting data in some cases (e.g. when a force comes from another solver with CoSimulation but an extra force must be prescribed directly from the project parameters).

## Problem Data
In this part, the user can determine the simulation time. It is the same as the problem data from the other solvers. 

```json
{
    "problem_data": {
        "problem_name": "problem_name",
        "start_time": 0.0,
        "end_time": 1.0,
        "echo_level" : 0

    },
```
## Solver Settings
```json
    "solver_settings": {
        "domain_size": 3,
        "echo_level": 1,
        "buffer_size": 3,
```
### Model Import settings
The Rigid Body solver do not need any additional input file. The model are described in this json file, with the active degree of freedoms configurated [here](#setting-up-the-degree-of-freedoms). For normal simulation, the ```"model_import_settings"``` need only the ```"input_type"``` to be set as ```"none"```. 

For restarting purposes, set the ```"input_type"``` to be ```"rest"``` and the needed information regarding the location of the restart file.
```json
        "model_import_settings": {
            "input_type": "none",
            "input_filename": "Main",
            "restart_load_file_label": "9.0",
            "input_output_path": "restart/"
        },
```
### Time integration parameters
This part sets up the time integration parameter. The Rigid Body solver uses the **Generalized Alpha Time Integration Scheme**.
```json
        "time_integration_parameters": {
            "rho_inf": 1,
            "time_step": 0.005
        },
```
### Setting up the degree of freedoms
This part sets up the active degree of freedom and its structural parameters. As mentioned before, there are three displacements and three rotations and it is not suitable to have more than one active rotation at a time. The user can adjust the degree of fredom to be active or not by setting the ```"constrained"``` parameter to be ```false```. If the degree of fredom is not listed in ```"active_dofs"```, it is set to be inactive by default.
```json

        "active_dofs": [
            {
                "dof": "displacement_y",
                "constrained": false,
                "system_parameters": {
                    "mass": 3.0,
                    "stiffness": 200,
                    "damping": 0.3
                }
            }
        ]
    },
```
## Applying Boundary Conditions
The boundary conditions are set in the ```"processes"``` field. Here the user can manage the body forces, apply forces and displacement.
### Body Force
The body force is managed in the ```"gravity"``` parameter. The body force is applied to the Rigidbody submodel part. It is applied to the ```"BODY_FORCE"``` variables, with the direction as the index. An example of a body force applied in the Y direction is shown below.
```json
    "processes": {
        "gravity": [
            {
                "python_module": "process_factory",
                "kratos_module": "KratosMultiphysics",
                "process_name": "ApplyConstantScalarValueProcess",
                "Parameters": {
                    "model_part_name": "Main.RigidBody",
                    "variable_name": "BODY_FORCE_Y",
                    "is_fixed": true,
                    "value": -981
                }
            }
        ],
```
### Initial condition
The user can apply initial conditions such as initial displacement and velocity in the ```"initial_conditions_process_list"``` parameter. An example of initial displacement in the X direction is shown below.
```json
        "initial_conditions_process_list": [
            {
                "python_module": "process_factory",
                "kratos_module": "KratosMultiphysics",
                "process_name": "ApplyConstantScalarValueProcess",
                "Parameters": {
                    "model_part_name": "Main.RigidBody",
                    "variable_name": "DISPLACEMENT_X",
                    "value": 1
                }
            }
        ],
```
### Prescribed Forces and Displacements
The ```"boundary_conditions_process_list"``` list all the prescribed forces and displacements. The PRESCRIBED_FORCE, PRESCRIBED_MOMENT, PRESCRIBED_DISPLACEMENT and PRESCRIBED_ROTATION are to be used here for the user prescribed variables. The solver also has FORCE and MOMENT in the Rigidbody submodel part. These variable are to be used in a CoSimulation, where the forces from another solver are transfered there. If the user has prescribed forces, it will be added to the forces obtained from another solver. Here are the examples for the prescribed displacement and force. The user can set a constant value or a function for the prescribed values.
```json
        "boundary_conditions_process_list": [
            {
                "python_module": "assign_vector_variable_process",
                "kratos_module": "KratosMultiphysics",
                "process_name": "AssignVectorVariableProcess",
                "Parameters": {
                    "model_part_name": "Main.RootPoint",
                    "variable_name": "PRESCRIBED_DISPLACEMENT",
                    "interval": [
                        3.76991,
                        "End"
                    ],
                    "constrained": [
                        false,
                        false,
                        false
                    ],
                    "value": [
                        "0.05*sin(20*t)",
                        0.0,
                        0.0
                    ]
                }
            },
            {
                "python_module": "assign_vector_variable_process",
                "kratos_module": "KratosMultiphysics",
                "process_name": "AssignVectorVariableProcess",
                "Parameters": {
                    "model_part_name": "Main.RigidBody",
                    "variable_name": "PRESCRIBED_FORCE",
                    "interval": [
                        3.14159,
                        "End"
                    ],
                    "constrained": [
                        false,
                        false,
                        false
                    ],
                    "value": [
                        0.0,
                        0.0,
                        "100*sin(10*t)"
                    ]
                }
            }
        ],
```
### auxiliar_process_list
```json
        "auxiliar_process_list": []
    }
}
```
## Output Processes
The output of the solver are managed here.
### Variable output
```json
    "output_processes": [
        {
            "python_module": "point_output_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "PointOutputProcess",
            "Parameters": {
                "model_part_name": "Main.RigidBody",
                "entity_type": "node",
                "interval": [
                    0.0,
                    "End"
                ],
                "position": [
                    0,
                    0,
                    0
                ],
                "output_variables": [
                    "DISPLACEMENT",
                    "ROTATION"
                ],
                "output_file_settings": {
                    "file_name": "Displacement",
                    "output_path": "results/FSI_RBS"
                }
            }
        },
```
### Creating a Restart File

```json
        {
            "python_module": "save_restart_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "SaveRestartProcess",
            "Parameters": {
                "model_part_name": "Main",
                "echo_level": 0,
                "serializer_trace": "no_trace",
                "restart_save_frequency": 1.0,
                "restart_control_type": "time",
                "save_restart_files_in_folder": true,
                "output_path": "restart/",
                "max_files_to_keep": 20
            }
        }
    ]
}
```
moment of intertia of the other rotation should be effected 

around the ther axis changes and it is not considered in the RBS