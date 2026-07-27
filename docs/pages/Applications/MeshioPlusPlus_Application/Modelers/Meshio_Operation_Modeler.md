---
title: Meshio Operation Modeler
keywords: meshio meshioplusplus modeler operation clean refine partition split
tags: [Meshio_Operation_Modeler.md]
sidebar: meshioplusplus_application
summary: Applies a meshio++ mesh operation to a model part from JSON settings.
---

## Overview

`MeshioOperationModeler` (`modelers.meshio_operation_modeler.MeshioOperationModeler`) is the data-driven, JSON-configurable entry point to the operations exposed by `MeshioPlusPlusMeshOperations` (see [Utilities](../Utilities/Mesh_Operations.html) for the full list and the settings each operation reads).

The operation is selected with `"operation"`, and its own settings are passed through `"operation_settings"`. Query what a build supports with `KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.GetSupportedOperations()`, and the available settings with `GetDefaultParameters()`.

## Usage

```json
{
    "modelers" : [{
        "modeler_name" : "KratosMultiphysics.MeshioPlusPlusApplication.modelers.meshio_operation_modeler.MeshioOperationModeler",
        "Parameters"   : {
            "operation"              : "clean",
            "input_model_part_name"  : "Structure",
            "output_model_part_name" : "StructureClean",
            "operation_settings"     : {
                "weld"      : true,
                "tolerance" : 1e-6
            }
        }
    }]
}
```

## Result

The result is written into `"output_model_part_name"`. The multi-output operations (`"split"`, `"partition"`) instead create one model part per piece, named `"<output_model_part_name>_<operation>_<index>"`; their names are listed in the report (`echo_level > 0` prints it).

## Parameters

| Setting | Default | Description |
|---|---|---|
| `operation` | `"clean"` | The operation to apply; see [Utilities](../Utilities/Mesh_Operations.html). |
| `input_model_part_name` | `""` *(required)* | The model part to operate on. |
| `output_model_part_name` | `""` *(required)* | The model part receiving the result; created when the modeler is instantiated. |
| `operation_settings` | `{}` | The operation's own settings, forwarded as-is to `MeshioPlusPlusMeshOperations.Execute`. |
| `echo_level` | `0` | Verbosity; `> 0` prints the operation's report. |
