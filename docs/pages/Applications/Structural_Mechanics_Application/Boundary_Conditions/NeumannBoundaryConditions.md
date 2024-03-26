---
title: Neumann Boundary Conditions
keywords: json
tags: [NeumannBoundaryConditions.md, NBC, Neumann]
sidebar: structural_mechanics_application
summary:
---
# Overview
The set of different Neumann boundary conditions are explained in this document. Each loading condition requires a `python_module` for corresponding process. These are explained in the file [Python Modules for Processes](./PythonModulesForProcesses.html).

## Point Load:
- Variable Name: `POINT_LOAD`
- A point load is a concentrated force applied at a specific location on a structure or object. It is represented as a single force vector acting at a particular point.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/images/point_load.png)

### Example:
```json
"loads_process_list"       : [{
        "python_module" : "assign_vector_by_direction_to_condition_process",
        "kratos_module" : "KratosMultiphysics",
        "check"         : "DirectorVectorNonZero direction",
        "process_name"  : "AssignVectorByDirectionToConditionProcess",
        "Parameters"    : {
            "model_part_name" : "Structure.PointLoad2D_load",
            "variable_name"   : "POINT_LOAD",
            "interval"        : [0.0,"End"],
            "modulus"         : 1.0,
            "direction"       : [0.0,-1,0.0]
        }
    }],
```
- You can find the above example here: [Point Load Example](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/Examples/point_load.zip)

## Point Moment:
- Variable name: `POINT_MOMENT`
- A point moment is a measure of rotational force around a specific point. It quantifies the tendency of a force to rotate an object about a given point and is calculated as the product of the force applied and the perpendicular distance from the point of rotation to the line of action of the force.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/images/point_moment.png)

### Example:
```json
"loads_process_list" : [{
        "python_module" : "assign_vector_by_direction_to_condition_process",
        "kratos_module" : "KratosMultiphysics",
        "check"         : "DirectorVectorNonZero direction",
        "process_name"  : "AssignVectorByDirectionToConditionProcess",
        "Parameters"    : {
            "model_part_name" : "Structure.PointMoment3D_moment",
            "variable_name"   : "POINT_MOMENT",
            "interval"        : [0.0,"End"],
            "modulus"         : 1.0,
            "direction"       : [0.0,0.0,1]
        }
    }],
```
- You can find the above example here: [Point Moment Example](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/Examples/point_moment.zip)

## Line Load:
- Variable name: `LINE_LOAD`
- A line load is a uniformly distributed load applied over a line, rather than concentrated at a single point.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/images/line_load.png)

### Example:
```json
"loads_process_list"       : [{
        "python_module" : "assign_vector_by_direction_to_condition_process",
        "kratos_module" : "KratosMultiphysics",
        "check"         : "DirectorVectorNonZero direction",
        "process_name"  : "AssignVectorByDirectionToConditionProcess",
        "Parameters"    : {
            "model_part_name" : "Structure.LineLoad2D_load",
            "variable_name"   : "LINE_LOAD",
            "interval"        : [0.0,"End"],
            "modulus"         : 1.0,
            "direction"       : [0.0,-1.0,0.0]
        }
    }],
```
- You can find the above example here: [Line Load Example](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/Examples/line_load.zip)

## Surface Load:
- Variable name: `SURFACE_LOAD`
- A surface load is a distributed force applied over a specific area, rather than concentrated at a single point or along a line.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/images/surface_load.png)

### Example:
```json
"loads_process_list"       : [{
        "python_module" : "assign_vector_by_direction_to_condition_process",
        "kratos_module" : "KratosMultiphysics",
        "check"         : "DirectorVectorNonZero direction",
        "process_name"  : "AssignVectorByDirectionToConditionProcess",
        "Parameters"    : {
            "model_part_name" : "Structure.SurfaceLoad3D_surface_load",
            "variable_name"   : "SURFACE_LOAD",
            "interval"        : [0.0,"End"],
            "modulus"         : 1000.0,
            "direction"       : [0.0,-1,0.0]
        }
    },
```
- You can find the above example here: [Surface Load Example](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/Examples/surface_load.zip)

## Pressure:
- Pressure can be described as Positive or Negative face pressure.
- Positive Face Pressure:
    - Variable name: `POSITIVE_FACE_PRESSURE`
    - This means applying pressure on the side of the surface that points outward. When applying positive face pressure, the direction of the pressure force is opposite to the surface normal. It is represented by the vector (1, -1), indicating a force acting away from the surface.
- Negative Face Pressure:
    - Variable name: `NEGATIVE_FACE_PRESSURE`
    - This refers to applying pressure on the side of the surface that points inward. When applying negative face pressure, the direction of the pressure force aligns with the surface normal. It is represented by the vector (-1, 1), indicating a force acting in the same direction as the surface normal, effectively pushing into the surface.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/images/pressure_negative_positive.png)

### Example:
```json
"loads_process_list"       : [{
        "python_module" : "assign_scalar_variable_to_conditions_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"  : "AssignScalarVariableToConditionProcess",
        "Parameters"    : {
            "model_part_name" : "Structure.SurfacePressure3D_pressure_surface",
            "variable_name"   : "POSITIVE_FACE_PRESSURE",
            "interval"        : [0.0,"End"],
            "value"           : 30000000.0
        }
    }]

"loads_process_list"       : [{
        "python_module" : "assign_scalar_variable_to_conditions_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"  : "AssignScalarVariableToConditionProcess",
        "Parameters"    : {
            "model_part_name" : "Structure.SurfacePressure3D_pressure_surface",
            "variable_name"   : "NEGATIVE_FACE_PRESSURE",
            "interval"        : [0.0,"End"],
            "value"           : 30000000.0
        }
    }]
```
- You can find the above examples here:
    - [Positive Face Pressure Example](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/Examples/positive_face_pressure.zip)
    - [Negative Face Pressure Example](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/Examples/negative_face_pressure.zip)