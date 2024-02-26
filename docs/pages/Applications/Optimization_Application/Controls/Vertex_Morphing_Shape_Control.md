---
title: Vertex Morphing Shape Control
keywords: 
tags: [Vertex_Morphing_Shape_Control.md]
sidebar: optimization_application
summary: 
---

## Introduction

```VertexMorphingShapeControl``` is used to control the shape of a specified domain (given by a model part). This control can use either [Explicit Vertex Morphing](../Filtering/Explicit_Vertex_Morphing.html) or Implicit Vertex Morphing as the filtering method. Therefore, the physical space gradients will be transformed to smoothened vertex morphed control space gradients.

This has an inbuilt mesh motion solver to solve for the mesh once the shape is changed. Therefore, it requires to have mesh motion settings as well.

## Json settings

### Explicit vertex morphing
Following json snippet illustrates one use case of explicit vertex morphing
```json
{
    "name": "explicit_shape_control",
    "type": "shape.vertex_morphing_shape_control",
    "module" : "KratosMultiphysics.OptimizationApplication.controls",
    "settings":{
        "controlled_model_part_names": [],
        "filter_settings": {
            "type" : "surface_explicit",
            "filter_function_type" : "linear",
            "damping_function_type" : "sigmoidal",
            "radius": 0.000000000001,
            "max_nodes_in_filter_radius": 1000
        },
        "mesh_motion_settings" : {},
        "output_all_fields": false,
        "fixed_model_parts": {},
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| name | A unique string name |
| type  | "shape.vertex_morphing_shape_control"  |
| module  | "KratosMultiphysics.OptimizationApplication.model_part_controllers"  |
| controlled_model_part_names | List of model part names of which the shape should be controlled. |
| filter_settings::type | "surface_explicit" |
| filter_settings::filter_function_type | Types of filter functions to be used. ["gaussian"](../Filtering/Explicit_Vertex_Morphing.html#gaussian-filter-function), ["linear"](../Filtering/Explicit_Vertex_Morphing.html#linear-filter-function), ["constant"](../Filtering/Explicit_Vertex_Morphing.html#constant-filter-function), ["cosine"](../Filtering/Explicit_Vertex_Morphing.html#cosine-filter-function), ["quartic"](../Filtering/Explicit_Vertex_Morphing.html#quartic-filter-function) are supported.|
| filter_settings::damping_function_type | Type of the damping function to be used. ["gaussian"](../Filtering/Explicit_Vertex_Morphing.md#gaussian-filter-function-1), ["linear"](../Filtering/Explicit_Vertex_Morphing.md#linear-filter-function-1), ["constant"](../Filtering/Explicit_Vertex_Morphing.md#constant-filter-function-1), ["cosine"](../Filtering/Explicit_Vertex_Morphing.md#constant-filter-function-1), ["quartic"](../Filtering/Explicit_Vertex_Morphing.md#quartic-filter-function-1) "sigmoidal" are supported. |
| filter_settings::radius| Filter radius |
| filter_settings::max_nodes_in_filter_radius| Number of max nodes to be found in the filter radius. This specifies maximum number of neighbours will be searched for. If the radii is higher, or mesh is refined, then this number should be increased. |
| mesh_motion_settings | Mesh motion solver settings |
| output_all_fields | Output intermediate results also to the ```OptimizationProblem``` data container. |
| fixed_model_parts | List of model part names to be dampened |

### Implicit vertex morphing
Following json snippet illustrates one use case of explicit vertex morphing
```json
{
    "name": "implicit_shape_control",
    "type": "shape.vertex_morphing_shape_control",
    "module" : "KratosMultiphysics.OptimizationApplication.controls",
    "settings":{
        "controlled_model_part_names": [],
            "filter_settings": {
                "type" : "bulk_surface_implicit",
                "radius": 0.000000000001,
                "linear_solver_settings" : {}
            },
            "mesh_motion_settings" : {},
            "output_all_fields": false,
            "fixed_model_parts": {},
    }
}
```

| Option | Allowed values |
| ------------- | ------------- |
| name | A unique string name |
| type  | "shape.vertex_morphing_shape_control"  |
| module  | "KratosMultiphysics.OptimizationApplication.model_part_controllers"  |
| controlled_model_part_names | List of model part names of which the shape should be controlled. |
| filter_settings::type | "bulk_surface_implicit" |
| filter_settings::radius| Filter radius |
| filter_settings::linear_solver_settings | Linear solver to be used in the implicit vertex morphing solver |
| mesh_motion_settings | Mesh motion solver settings |
| output_all_fields | Output intermediate results also to the ```OptimizationProblem``` data container. |
| fixed_model_parts | List of model part names to be dampened |

## Source files
* [Doxygen](TODO) TODO
* [applications/OptimizationApplication/python_scripts/controls/shape/vertex_morphing_shape_control.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/controls/shape/vertex_morphing_shape_control.py)
* [applications/OptimizationApplication/python_scripts/filtering/helmholtz_analysis.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/OptimizationApplication/python_scripts/filtering/helmholtz_analysis.py)


