---
title: Neumann Boundary Conditions
keywords: json
tags: [NeumannBoundaryConditions.md]
sidebar: structural_mechanics_application
summary:
---
# Overview
The set of different Neumann boundary conditions are explained in this document. Each loading condition requires a `process_name`. These are explained in the process name's documentation.

## Point Load:
- Variable Name: `POINT_LOAD`
- A point load is a concentrated force applied at a specific location on a structure or object. It is represented as a single force vector acting at a particular point.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/point_load.png)

### Example:
```json
"loads_process_list"       : [{
        "python_module" : "assign_vector_by_direction_to_condition_process",
        "kratos_module" : "KratosMultiphysics",
        "check"         : "DirectorVectorNonZero direction",
        "process_name"  : "AssignVectorByDirectionToConditionProcess",
        "Parameters"    : {
            "model_part_name" : "Structure.PointLoad3D_BC_NEUMANN",
            "variable_name"   : "POINT_LOAD",
            "modulus"         : 10.0,
            "direction"       : [0.0,-1,0.0],
            "interval"        : [0.0,"End"]
        }
    }],
```

## Point Moment:
- Variable name: `POINT_MOMENT`
- A point moment is a measure of rotational force around a specific point. It quantifies the tendency of a force to rotate an object about a given point and is calculated as the product of the force applied and the perpendicular distance from the point of rotation to the line of action of the force.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/point_moment.png)

### Example:
```json
"loads_process_list" : [{
        "python_module" : "assign_vector_variable_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"  : "AssignVectorVariableProcess",
        "Parameters"    : {
                "mesh_id"         : 0,
                "model_part_name" : "Structure.PointLoad3D_neumann",
                "constrained"     : false,
                "variable_name"   : "POINT_MOMENT",
                "value"           : [0.0,0.0,"25000*t"]
        }
    }],
```

## Line Load:
- Variable name: `LINE_LOAD`
- A line load is a uniformly distributed load applied over a line, rather than concentrated at a single point.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/line_load.png)

### Example:
```json
"loads_process_list"       : [{
        "python_module" : "assign_vector_by_direction_to_condition_process",
        "kratos_module" : "KratosMultiphysics",
        "help"          : "This process sets a vector variable value over a condition",
        "check"         : "DirectorVectorNonZero direction",
        "process_name"  : "AssignVectorByDirectionToConditionProcess",
        "Parameters"    : {
            "mesh_id"         : 0,
            "model_part_name" : "Structure.LineLoad2D_Load_on_lines_Auto1",
            "variable_name"   : "LINE_LOAD",
            "modulus"         : 1.0,
            "direction"       : [1.0,0.0,0.0],
            "interval"        : [0.0,"End"]
        }
    }],
```

## Surface Load:
- Variable name: `SURFACE_LOAD`
- A surface load is a distributed force applied over a specific area, rather than concentrated at a single point or along a line.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/surface_load.png)

### Example:
```json
"loads_process_list"       : [{
        "python_module" : "assign_vector_by_direction_to_condition_process",
        "kratos_module" : "KratosMultiphysics",
        "help"          : "This process sets a vector variable value over a condition",
        "check"         : "DirectorVectorNonZero direction",
        "process_name"  : "AssignVectorByDirectionToConditionProcess",
        "Parameters"    : {
            "mesh_id"         : 0,
            "model_part_name" : "Structure.SurfaceLoad3D_Load_on_surfaces_Auto2",
            "variable_name"   : "SURFACE_LOAD",
            "modulus"         : 1.0,
            "direction"       : [0.0,1.0,0.0],
            "interval"        : [0.0,"End"]
        }
    },
```

## Pressure:
- Pressure can be described as Positive or Negative face pressure.
- Positive Face Pressure:
    - Variable name: `POSITIVE_FACE_PRESSURE`
    - This means applying pressure on the side of the surface that points outward. When applying positive face pressure, the direction of the pressure force is opposite to the surface normal. It is represented by the vector (1, -1), indicating a force acting away from the surface.
- Negative Face Pressure:
    - Variable name: `NEGATIVE_FACE_PRESSURE`
    - This refers to applying pressure on the side of the surface that points inward. When applying negative face pressure, the direction of the pressure force aligns with the surface normal. It is represented by the vector (-1, 1), indicating a force acting in the same direction as the surface normal, effectively pushing into the surface.

![StructuralMechanicsApplication](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/StructuralMechanicsApplication/pressure_negative_positive.png)

### Example:
```json
"loads_process_list"       : [{
        "python_module" : "assign_scalar_variable_to_conditions_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"  : "AssignScalarVariableToConditionProcess",
        "Parameters"    : {
            "mesh_id"         : 0,
            "model_part_name" : "Structure.SurfacePressure3D_Surface_Q4_thick",
            "variable_name"   : "POSITIVE_FACE_PRESSURE",
            "value"           : 2040.37,
            "interval"        : [0.0,"End"]
        }
    },{
        "python_module" : "assign_scalar_variable_to_conditions_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"  : "AssignScalarToConditionsProcess",
        "Parameters"    : {
            "mesh_id"         : 0,
            "model_part_name" : "Structure.SurfacePressure3D_surface",
            "variable_name"   : "NEGATIVE_FACE_PRESSURE",
            "value"           : -2040.37,
            "interval"        : [0.0,"End"]
        }
    }]
```