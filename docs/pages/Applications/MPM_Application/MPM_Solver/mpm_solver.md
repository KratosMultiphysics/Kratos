---
title: MPM Solver (Base)
keywords: mpm constitutive laws
tags: [mpm constitutive laws]
sidebar: mpm_application
summary: 
---

## Overview

The solver is responsible for handling the simulation physics and solving the problem. We should remark the following points, that apply to all Kratos Applications:
* solvers are always implemented at the Python level;
* all solvers extend the base `PythonSolver` class, which can be found in [this script](https://github.com/KratosMultiphysics/Kratos/kratos/python_scripts/python_solver.py). This base solver does not include any physics and its purpose is to define the solver class API.

Each Kratos Multiphysics application may implement multiple solvers, each tailored to different types of problems or simulation strategies. For example, the `MPMApplication` implements the following solvers:
* [`MPMImplicitDyanamicSolver`](./mpm_implicit_solver) - solves time-dependent problems by means of an implicit time scheme;
* `MPMExplicitSolver` - solves time-dependent problems using an explicit time scheme;
* `MPMQuasiStaticSolver` - solves quasi-static problems;
* `MPMStaticSolver` - solves static problems.

In Kratos, it is common practice to define a base solver for a specific physics-based application and then derive specialized solvers from it. For instance, all MPM solvers listed before extend the base `MPMSolver`, which manages common aspects of MPM-based simulations. The `MPMSolver` itself extendes the base `PythonSolver` class.

The solver is also responsible for creating and managing the solution strategy (e.g., the Newton-Raphson algorithm for non-linear problems) and the linear solver (if required) to solve the system of equations.

The solver is initialized through a custom `AnalysisStage` class, which uses (part of) the parameters from the `solver_settings` section of the `ProjectParameters.json` input file (introduced [here](../Input_Files/json#projectparametersjson)). The `solver_settings` also contain other parameters used by the solver to define the linear solver and the solution strategy.

Once initialised, the solver can be accessed from the `AnalysisStage` class by calling the `_GetSolver()` method.

## Parameters & Defaults

The most important parameters contained in the `solver_settings` section of the `ProjectParameters.json` input file and required for the initialisation and configuration of MPM solvers extending the `MPMSolver` class are the following.

```json
{
    "solver_settings"  : {
        "model_part_name"                 : "MPM_Material",
        "echo_level"                      : 0,
        "domain_size"                     : 2,
        "solver_type"                     : "Dynamic",
        "time_stepping"                   : {
            "time_step" : 0.01
        },
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
        "analysis_type"                   : "non_linear",
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

* `model_part_name`: string with the name to be assigned to the material point model part.
* `echo_level`: an integer controlling the verbosity level of the simulation output. Higher values result in more messages printed to the standard output.
* `domain_size`: size of the domain in which the problem is set; admissible values are `2` (two-dimensional problems) and `3` (three-dimensional problems).
* `solver_type`: string used to initialize the appropriate MPM solver class depending on the type of problem (dynamic, quasi-static, static) to be solved. Admissible values are:
    * `"dynamic"`: initializes either the `MPMImplicitDyanamicSolver` or the `MPMExplicitSolver`, depending on whether the additional parameter `time_integration_method` is set to `"implicit"` or `"explicit"`, respectively;
    * `"quasi-static"`: initializes the `MPMQuasiStaticSolver`;
    * `"static"`: initializes the `MPMStaticSolver`.

  The following table summarizes the available solvers and the parameters used to initialise them

| MPM Solver                  | `solver_type`    | `time_integration_method`  |
|-----------------------------| :--------------: | :------------------------: |
| `MPMImplicitDyanamicSolver` | `"dynamic"`      | `"implicit"`               |
| `MPMExplicitSolver`         | `"dynamic"`      | `"explicit"`               |
| `MPMQuasiStaticSolver`      | `"quasi-static"` | *not required*             |
| `MPMStaticSolver`           | `"static"`       | *not required*             |

* `time_stepping`: object defining the time step (if required).
* `grid_model_import_settings`: object containing the following name-value pairs:
    * `input_filename`: name (without extension) of the file defining the background grid
    * `input_type`: input file type.
* `model_import_settings`: object containing the following name-value pairs:
    * `input_filename`: name (without extension) of the file defining the mesh discretising the body to be simulated by means of the MPM
    * `input_type`: input file type.
* `material_import_settings`: name of the json file (typically `ParticleMaterials.json`) containing the information about the physical properties of the materials.
* `analysis_type`: admissible values are `"linear"` and `"non_linear"`. If set to `"non_linear"`, a linearization procedure based on the Newton-Raphson strategy is used. In this case, the following additional parameters defining the convergence criterion of the Newton-Raphson algorithm must be defined.
    * `convergence_criterion`: type of convergence criterion to be used; available options are `"displacement_criterion"` and `"residual_criterion"`.
    * `displacement_relative_tolerance`: relative tolerance value for displacement-based convergence criterion.
    * `displacement_absolute_tolerance`: absolute tolerance value for displacement-based convergence criterion.
    * `residual_relative_tolerance`: relative tolerance value for residual-based convergence criterion.
    * `residual_absolute_tolerance`: absolute tolerance value for residual-based convergence criterion.
    * `max_iteration`: maximum number of non-linear iterations.
* `pressure_dofs`: boolean variable identifying the type of MPM formulation to be used (displacement based if `false`, mixed displacement-pressure if `true`).
* `linear_solver_settings`: object defining the type of solver to be used for solving the linear system.
* `auxiliary_variables_list`: list of strings, each one representing an additional variable required in the simulation process.

## Source Code

[<i class="fa fa-github"></i> `MPMApplication/python_scripts/mpm_solver.py`](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_solver.py)
