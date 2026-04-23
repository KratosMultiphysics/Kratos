---
title: json files
keywords: mpm json
tags: [mpm json]
sidebar: mpm_application
summary: 
---

In this section we discuss the role played by the `ProjectParameters.json` and `ParticleMaterials.json` files, which define, respectively, the settings and the material properties of the problem to be solved.

The [`JSON` format](https://en.wikipedia.org/wiki/JSON) is an open standard file format and data interchange format that uses human-readable text to store and transmit data objects consisting of name–value pairs and arrays. Kratos uses a thin wrapper arround this syntax, the `Parameters` class.

A detailed description of the `Parameters` class can be found [here](../../../Kratos/For_Users/Crash_Course/Input_Output_and_Visualization/JSON_Configuration_File), along with [a tutorial on how to use it to read `JSON` files](../../../Kratos/For_Users/Crash_Course/Input_Output_and_Visualization/Project_Parameters).


## `ProjectParameters.json`

The input file `ProjectParameters.json` has the following structure

```json
{
    "analysis_stage"   : "KratosMultiphysics.MPMApplication.mpm_analysis",
    "problem_data"     : {
        "problem_name"  : "name",
        "parallel_type" : "OpenMP",
        "echo_level"    : 0,
        "start_time"    : 0.0,
        "end_time"      : 1
    },
    "solver_settings"  : {
        ...
    },
    "processes"        : {
        "constraints_process_list" : [],
        "loads_process_list"       : [],
        "list_other_processes"     : [],
        "gravity"                  : []
    },
    "output_processes" : {
        "gid_output_processes"   : [],
        "vtk_output_processes"   : [],
        "other_output_processes" : []
    }
}
```

The `ProjectParameters.json` typically contains five main blocks:

* `analysis_stage`, a string that specifies the Python module containing the implementation of the `MpmAnalysis` class. This class extends the `AnalysisStage` class from the Kratos Core and defines the overall workflow of the simulation.

* `problem_data`, a set of parameters defining the general settings for the Kratos simulation:

    - `problem_name`: a label that identifies the problem to be solved.

    - `parallel_type`: a label that specifies the parallelization method used for solving the problem. Currently, the `MPMApplication` only supports shared memory parallelization, meaning that the only valid value is `"OpenMP"`.

    - `echo_level`: an integer that controls the verbosity level of the simulation output. Higher values result in more messages printed to the standard output, providing greater detail on the simulation's status.
        * Minimum: `0` - silent
        * Maximum: `4` - detailed output

    -  `start_time`: a number defining the starting time of the simulation.

    - `end_time`: a number defining the ending time of the simulation.

* `solver_settings`, a set of attribute-value pairs defining the settings for the solvers, like analysis type, linear solver, etc. These parameters are discussed in detail [here](../MPM_Solver/mpm_solver).
* `processes`, a set of lists containing the processes to apply, e.g., boundary and initial conditions. Within the `MPMApplication`, processes are usually divided into four lists.
    1. `constraints_process_list`: this list contains all the processes imposing Dirichlet boundary conditions;
    2. `loads_process_list`: this list contains all the processes defining the external loads acting on the body;
    3. `list_other_processes`: this list contains other processes not related to boundary conditions;
    4. `gravity`: this list contains the [`AssignGravityToMaterialPoints`](../Processes/assign_gravity) process.

   Other custom lists can also be defined. A detailed description of the processes implemented within the `MPMApplication` can be found [here](../Processes/Overview.md).
* `output_processes`, a set of lists containing the output processes for the post-processing of the results. Within the `MPMApplication`, output processes are usually divided into three lists.
   1. `gid_output_processes`: this list contains the [`GiDOutputProcess`](../../../Kratos/Processes/Output_Process/GiD_Output_Process) and the [`MPMGiDOutputProcess`](../Output_Processes/mpm_gid_output_process) for exporting data that can be visualized with [GiD](https://www.gidsimulation.com/);
   2. `vtk_output_processes`: this lists contains the [`VtkOutputProcess`](../../../Kratos/Processes/Output_Process/VTK_Output_Process) and the [`MPMVtkOutputProcess`](../Output_Processes/mpm_vtk_output_process) for exporting data in the `vtk` format that can be visualized using post-processing tools such as [Paraview](https://www.paraview.org/) and [VisIt](https://visit-dav.github.io/visit-website/index.html);
   3. `other_output_processes`: this lists contains other output processes.

   A detailed description of all the output processes that can be used with the `MPMApplication` can be found [here](../Output_Processes/Overview.md).

## `ParticleMaterials.json`

The input file `ParticleMaterials.json` has the following structure

```json
{
    "properties" : [{
        "model_part_name" : "Initial_MPM_Material.Material_Domain_1",
        "properties_id"   : 1,
        "Material"        : {
            "constitutive_law" : {
                "name" : "LinearElasticIsotropicPlaneStrain2DLaw"
            },
            "Variables"        : {
                "THICKNESS"                   : 1.0,
                "MATERIAL_POINTS_PER_ELEMENT" : 1,
                "DENSITY"                     : 7850.0,
                "YOUNG_MODULUS"               : 206900000000.0,
                "POISSON_RATIO"               : 0.29
            },
            "Tables"           : null
        }
    },{
        "model_part_name" : "Initial_MPM_Material.Material_Domain_2",
        "properties_id"   : 2,
        "Material"        : {
            "constitutive_law" : {
                "name" : "HenckyMCPlasticPlaneStrain2DLaw"
            },
            "Variables"        : {
                "THICKNESS"                   : 1.0,
                "MATERIAL_POINTS_PER_ELEMENT" : 1,
                "DENSITY"                     : 7850.0,
                "YOUNG_MODULUS"               : 206900000000.0,
                "POISSON_RATIO"               : 0.29,
                "COHESION"                    : 0.0,
                "INTERNAL_FRICTION_ANGLE"     : 0.5235987755982,
                "INTERNAL_DILATANCY_ANGLE"    : 0.0
            },
            "Tables"           : null
        }
    },{
        "model_part_name" : "Initial_MPM_Material.Material_Domain_3",
        "properties_id"   : 3,
        "Material"        : {
            "constitutive_law" : {
                "name" : "HyperElasticNeoHookeanPlaneStrain2DLaw"
            },
            "Variables"        : {
                "THICKNESS"                   : 1.0,
                "MATERIAL_POINTS_PER_ELEMENT" : 1,
                "DENSITY"                     : 7850.0,
                "YOUNG_MODULUS"               : 206900000000.0,
                "POISSON_RATIO"               : 0.29
            },
            "Tables"           : null
        }
    }]
}
```

The `property` field contains a list of one or more objects, each defining the constitutive law associated to a submodel part of the `Intitial_MPM_Material` model part.

For example, the `ParticleMaterials.json` file shown in the code-block above refers to a problem in which three continuum bodies are discretised by means of material points, each body being assigned a different constitutive law.

Each json object in the list assigned to the `property` field is defined through the following attributes.
* `model_part_name`: name of the submodel part to which the properties in that object are assigned;
* `properties_id`: integer identifying that specific property;
* `Material`: json object providing detailed information about the constitutive law. It contains the following name-value pairs.
    - `constitutive_law`: object containing a string identifying the [constitutive law](../Constitutive_Laws/constitutive_laws);
    - `Variables`: object defining some variables
        * `THICKNESS`: thickness for 2D problems
        * `MATERIAL_POINTS_PER_ELEMENT`: number of material points to be used on each element of the mesh discretising the submodel part indentified by `model_part_name`. The discretisation of the continuum to simulate by means of material points is discussed in more details [here](./mdpa). Depending on the element type, the following number of material points per element can be used:

            | Triangles | Tetrahedra | Quadrilaterals | Hexahedra |
            | :-------: | :--------: | :------------: | :-------: |
            |  1        |  1         |  1             |   1       |
            |  3        |  4         |  4             |   8       |
            |  4        |  8         |  9             |  27       |
            |  6        | 14         | 16             |  64       |
            | 12        | 24         | 25             | 125       |

        * other variables depending on the type of constitutive law that has been chosen.

A detailed discussion of the consitutive laws implemented within the `MPMApplication` and the list of variables to be added to the object assigned to the `Variables` attribute can be found [here](../Constitutive_Laws/constitutive_laws).

