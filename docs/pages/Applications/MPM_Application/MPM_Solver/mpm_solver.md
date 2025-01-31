---
title: MPM Solver
keywords: mpm constitutive laws
tags: [mpm constitutive laws]
sidebar: mpm_application
summary: 
---

```json
{
    "solver_settings"  : {
        "time_stepping"                   : {
            "time_step" : 0.01
        },
        "domain_size"                     : 2,
        "echo_level"                      : 0,
        "solver_type"                     : "Dynamic",
        "time_integration_method"         : "implicit",
        "scheme_type"                     : "newmark",
        "analysis_type"                   : "non_linear",
        "model_part_name"                 : "MPM_Material",
        "grid_model_import_settings"      : {
            "input_type"     : "mdpa",
            "input_filename" : "file_name_Grid"
        },
        "model_import_settings"           : {
            "input_type"     : "mdpa",
            "input_filename" : "file_name_Body"
        },
        "material_import_settings"        : {
            "materials_filename" : "ParticleMaterials.json"
        },
        "convergence_criterion"           : "residual_criterion",
        "displacement_relative_tolerance" : 0.0001,
        "displacement_absolute_tolerance" : 1e-9,
        "residual_relative_tolerance"     : 0.0001,
        "residual_absolute_tolerance"     : 1e-9,
        "max_iteration"                   : 10,
        "pressure_dofs"                   : false,
        "linear_solver_settings"          : {
            "solver_type" : "LinearSolversApplication.sparse_lu"
        },
        "auxiliary_variables_list"        : ["NORMAL","IS_STRUCTURE"]
    }
}
```
