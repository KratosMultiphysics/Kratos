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

## `ParticleMaterials.json`

The input file `ParticleMaterials.json` has the following structure

```json
{
    "properties" : [{
        "model_part_name" : "Initial_MPM_Material.Parts_Material_domain_Material_domain_Auto1",
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

