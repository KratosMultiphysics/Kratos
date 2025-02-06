---
title: Overview
keywords: process mpm output
tags: [output mpm process]
sidebar: mpm_application
summary: 
---

In `KratosMultiphysics`, boundary and initial conditions are typically imposed by means of Python classes extending the base [`Process`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/processes/process.h) class, which is implemented in the `KratosCore`.
More details about the `Process` class can be found the [dedicated documentation webpage](../../../Kratos/Processes/process) and in the [Kratos Crash Course](../../../Kratos/For_Users/Crash_Course/5_Simulation_Loop). The parameters required by python processes must be defined in the `"processes": {}` section of the `ProjectParameters.json` file (more details about this input file can be found [here](../Input_Files/json#projectparametersjson)).


In the `MPMApplication`, the following processes are implemented.

* [AssignGravityToMaterialPointProcess](./assign_gravity) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/assign_gravity_to_material_point_process.py)) - Process assigning gravity to the material point elements of an MPM submodel part
* [AssignInitialVelocityToMaterialPointProcess](./assign_initial_velocity) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/assign_initial_velocity_to_material_point_process.py)) - Process assigning an initial velocity to the material point elements of an MPM submodel part
* **Grid-based boundary conditions** (or conforming/fitted boundary conditions) - Processes imposing boundary conditions on the nodes of a background grid submodel part
    - [AssignVectorVariableProcess](./Grid-based_Boundary_Conditions/fixed_displacement_boundary_condition) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/assign_vector_variable_process.py)) - Process implemented in the Kratos Core for imposing a value to the nodes of a model part and here used for imposing a zero displacement boundary condition
    - [ApplyMPMSlipBoundaryProcess](./Grid-based_Boundary_Conditions/slip_boundary_condition) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/apply_mpm_slip_boundary_process.py)) - Process imposing zero displacement in the direction orthogonal to the conditions of a background grid submodel part (no constraints in the tangent direction)
* **Material Point-based boundary conditions** (or non-conforming/unfitted boundary conditions) - Processes imposing boundary conditions by means of (moving) boundary material points
    - [ApplyMPMParticleDirichletConditionProcess](./Material_Point-based_Boundary_Conditions/penalty) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/apply_mpm_particle_dirichlet_condition_process.py)) - Process imposing Dirichlet boundary conditions by means of the Penalty method
