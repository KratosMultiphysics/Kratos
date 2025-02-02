---
title: json files
keywords: mpm json
tags: [mpm json]
sidebar: mpm_application
summary: 
---

In this section we discuss the role played by the `ProjectParameters.json` and `ParticleMaterials.json` files, which define, respectively, the settings and the material properties of the problem to be solved.

`JSON` is an open-standard format that uses human-readable text to transmit data objects consisting of attribute–value pairs. Kratos uses a thin wrapper arround this syntax, the `Parameters` object.

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

    - `parallel_type`: a label that specifies the parallelization method used for solving the problem. Currently, the `MPMApplication` only supports shared memory parallelization, meaning the only valid value is `"OpenMP"`.

    - `echo_level`: an integer that controls the verbosity level of the simulation output. Higher values result in more messages printed to the standard output, providing greater detail on the simulation's status.
        * Minimum: `0` - silent
        * Maximum: `4` - detailed output

    -  `start_time`: a number defining the starting time of the simulation.

    - `end_time`: a number defining the ending time of the simulation.

* `solver_settings`, a set of attribute-value pairs defining the settings for the solvers, like analysis type, linear solver, etc.
* `processes`, a set of lists containing the processes to apply, e.g., boundary and initial conditions. A detailed description of the processes implemented within the `MPMApplication` can be found in [Processes: an Overview](../Processes/Overview.md).
* `output_processes`, a set of lists containing the output processes for the post-processing of the results. A detailed description of the output processes that can be used with the `MPMApplication` can be found in [Output Processes: an Overview](../Post_processing/Overview.md).

## `ParticleMaterials.json`

The input file `ParticleMaterials.json` has the following structure

```json
{
    "properties" : [{
        "model_part_name" : "Initial_MPM_Material.Parts_Material_domain",
        "properties_id"   : 1,
        "Material"        : {
            "constitutive_law" : {
                "name" : "LinearElasticIsotropicPlaneStrain2DLaw"
            },
            "Variables"        : {
                "THICKNESS"                   : 1.0,
                "MATERIAL_POINTS_PER_ELEMENT" : 6,
                "DENSITY"                     : 7850.0,
                "YOUNG_MODULUS"               : 206900000000.0,
                "POISSON_RATIO"               : 0.29
            },
            "Tables"           : null
        }
    }]
}
```

