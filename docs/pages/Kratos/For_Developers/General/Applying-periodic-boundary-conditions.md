---
title: Applying Periodic Boundary Conditions
keywords: 
tags: [Applying Periodic Boundary Conditions]
sidebar: kratos_for_developers
summary: 
---

In KratosMultiphysics, application of periodic boundary conditions can be done using the `apply_periodic_condition_process`. This process uses linear master-slave constraints to apply a periodic condition between the specified boundary (sub)modelparts. 

An example entry for such a process in the the processes list in the `ProjectParameters.json` is given below. 

```json
{
    "python_module": "apply_periodic_condition_process",
    "kratos_module": "KratosMultiphysics",
    "process_name": "ApplyPeriodicBoundaryConditionProcess",
    "Parameters": {
        "computing_model_part_name": "computing_domain",
        "model_part_name": "Structure",
        "first_model_part_name": "GENERIC_horiSurface",
        "second_model_part_name": "GENERIC_vertSurface",
        "interval": [0.0, 1e30],
        "variable_names": ["DISPLACEMENT"],
        "transformation_settings": {
            "rotation_settings": {
                "center": [0.0, 0.0, 0.0],
                "axis_of_rotation": [0.0, 0.0, 1.0],
                "angle_degree": 90.0
            }
        },
        "search_settings": {
            "max_results": 100,
            "tolerance": 1E-2
        }
    }
}
```

From the above json parameters :  `python_module` , `kratos_module`, `process_name` are always the same.
`Parameters` dictionary changes depending on the case. Explanation of the entries of `Parameters` dictionary is as follows : 

`computing_model_part_name` : Name of the computing modelpart from the solver.

`model_part_name` : This is the name of the modelpart which is specified in the `solver_settings` dictionary

Periodic boundary conditions is always applied to a pair of modelparts.

`first_model_part_name` : Name of the first (sub)modelpart in the pair. (internally the master modelpart for the master-slave relation)

`second_model_part_name` : Name of the second (sub)modelpart in the pair. (internally the slave modelpart for the master-slave relation)

`interval` : The time interval during which the periodic boundary condition is active.

`variable_names` : The list of variable for which are to be periodic on first and second (sub)modelparts. **IMPORTANT:** These variables should be Degree of Freedom (DOFs) for the given problem for the periodic condition to be effective.

`transformation_settings` : This defines how the modelparts with names `first_model_part_name` and `second_model_part_name` are geometrically related to each other. This can be of two types : rotational or translational. Rotational relation requires  `axis_of_rotation`, `center` and `angle` where as the translational relation requires : `dir_of_translation` and `magnitude` to be specified. 

The relationship between the nodes on the two specified (sub)modelparts is found using nearest element relation ship. That is nodes on the second (here slave) (sub)modelpart are projected on to the first (sbu)modelpart and are located inside a triangle/rectangle/line to find its corresponding master nodes and the master-slave relation is formed accordingly. For this search `search_settings` is used. 

`search_settings` : recommended only for advanced users, this contains the setting for the nearest neighbor search for formulating the master-slave relation. 

## Examples cases 

Two example cases are provided in the Examples repository under KratosMultiphyiscs. 
1.  A structural example can be found [here](https://github.com/KratosMultiphysics/Examples/tree/master/structural_mechanics/use_cases/periodic_bc_example). The example is of a disk with unit thickness fixed on its axis and is subjected to uniform and steady centrifugal force. Here instead of the whole disk only a quarter portion of the disk is modeled by specifying periodic boundary conditions.

_Problem Setup_
<img src="https://github.com/KratosMultiphysics/Examples/blob/master/structural_mechanics/use_cases/periodic_bc_example/data/centrifugal_force_vectors.jpg" width="300">


<img src="https://github.com/KratosMultiphysics/Examples/blob/master/structural_mechanics/use_cases/periodic_bc_example/data/result_periodic_viz.jpg" width="800">

_Reconstructed Solution of Half Disk_


2.  A fluid example can be found [here](https://github.com/KratosMultiphysics/Examples/tree/master/fluid_dynamics/use_cases/kelvin_helmholtz_instability). The example is a simulation of [Kelvin-Helmholtz instability](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Helmholtz_instability). This requires a infinitely long domain whose upper and lower halves move at different (usually opposite) velocities. In this example only a portion of the domain (unit square) is simulated with periodic boundary condition applied on the two sides.

<img src="https://github.com/KratosMultiphysics/Examples/blob/master/fluid_dynamics/use_cases/kelvin_helmholtz_instability/data/result_viz_with_more_eddies.gif" width="500">

_Instability Simulation Result(Vorticity magnitude)_

Here the instabilities can be seen forming at the middle of the domain and the fluid mixing.
