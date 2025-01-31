---
title: json files
keywords: mpm json
tags: [mpm json]
sidebar: mpm_application
summary: 
---

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

### `analysis_stage`
A string that specifies the Python module containing the implementation of the `MpmAnalysis` class. This class extends the `AnalysisStage` class from the Kratos Core and defines the overall workflow of the simulation.

### `problem_data`

##### `problem_name`
A label that identifies the problem to be solved.

##### `parallel_type`
Specifies the parallelization method used for solving the problem. Currently, the MPMApplication only supports shared memory parallelization, meaning the only valid value is "OpenMP".

##### `echo_level`
An integer that controls the verbosity level of the simulation output.
Higher values result in more messages printed to the standard output, providing greater detail on the simulation's status.
* Minimum: 0 - silent
* Maximum: 4 - detailed output

##### `start_time`
The starting time of the simulation.

##### `end_time`
The ending time of the simulation.

### `solver_settings`
### `processes`
### `output_processes`

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

