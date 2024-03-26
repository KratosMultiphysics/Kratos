---
title: Dirichlet Boundary Conditions
keywords: json
tags: [DirichletBoundaryConditions.md, DBC, Dirichlet]
sidebar: structural_mechanics_application
summary:
---
# Overview
The set of different Dirichlet boundary conditions are explained in this document. Each condition requires a `python_module` for corresponding process. These are explained in the file [Python Modules for Processes](./PythonModulesForProcesses.html).

## Displacement:

- Variable Name: `DISPLACEMENT`
- It specifies the displacement values or constraints at certain points or along certain edges or faces of a computational domain. These prescribed displacements can be either fixed (known values) or constrained to follow a specific function or relationship.
- The required parameters for the `variable_name: DISPLACEMENT` are `constrained` and `value`. Parameter `constrained` explains if the corresponding DOF is constrained in the directions [`X, Y, Z`]. Parameter `value` can be used to prescribe the fixed value.

### Example:
```json
"constraints_process_list" : [{
            "python_module" : "assign_vector_variable_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "AssignVectorVariableProcess",
            "Parameters"    : {
                "model_part_name" : "Structure.DISPLACEMENT_base",
                "variable_name"   : "DISPLACEMENT",
                "interval"        : [0.0,"End"],
                "constrained"     : [true,true,true],
                "value"           : [0.0,0.0,0.0]
            }
        }]
```
- You can find this constraint applied in the following example: [Line Load Example](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/Examples/line_load.zip)

- The displacement constraint in the corresponding direction can be directly applied by utilizing a `master_variable_name`. These are as follows: `DISPLACEMENT_X`, `DISPLACEMENT_Y`, and `DISPLACEMENT_Z`. (TODO: Needs to be verified!)

### Example:
```json
{
    "python_module" : "impose_rigid_movement_process",
    "kratos_module" : "KratosMultiphysics.StructuralMechanicsApplication",
    "Parameters"    : {
        "model_part_name"      : "Parts_Parts_Auto1",
        "master_variable_name" : "DISPLACEMENT_Z",
        "interval"             : [0.0,"End"]
    }
}
```
## Rotation:

- Variable Name: `ROTATION`
- Rotation constraint is used to fix the rotation of a certain point, line, or a surface of a structure, ensuring it remains rigid or constrained in its rotational behavior. It dictates the exact rotational values or constraints that the solution must satisfy along specified portions of the boundary.
- The required parameters for the `variable_name: ROTATION` are `constrained` and `value`. Parameter `constrained` explains if the corresponding DOF is constrained in the directions [`X, Y, Z`]. Parameter `value` can be used to prescribe the fixed value.

### Example:
```json
{
    "python_module" : "assign_vector_variable_process",
    "kratos_module" : "KratosMultiphysics",
    "process_name"  : "AssignVectorVariableProcess",
    "Parameters"    : {
        "model_part_name" : "Structure.ROTATION_fix",
        "variable_name"   : "ROTATION",
        "interval"        : [0.0,"End"],
        "constrained"     : [true,true,true],
        "value"           : [0.0,0.0,0.0]
    }
}
```
- You can find this constraint applied in the following example: [Point Moment Example](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/Examples/point_moment.zip)