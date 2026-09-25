---
title: Meshio Interpolate Modeler
keywords: meshio meshioplusplus modeler interpolate sample nearest barycentric
sidebar: meshioplusplus_application
tags: [Meshio_Interpolate_Modeler.md]
summary: Samples one model part's field data onto another's geometry from JSON settings.
---

## Overview

`MeshioInterpolateModeler` (`modelers.meshio_interpolate_modeler.MeshioInterpolateModeler`) is the data-driven, JSON-configurable entry point to `MeshioPlusPlusMeshOperations.Interpolate` (see [Utilities](../Utilities/Mesh_Operations.html#additional-entry-points)).

Unlike the [Meshio Operation Modeler](Meshio_Operation_Modeler.html), this modeler needs two independent sources: `"source_model_part_name"` (carrying the field data to sample) and `"target_model_part_name"` (contributing only its geometry) — which is why `"interpolate"` is a separate modeler rather than another `"operation"` entry.

Field data is selected with the same settings `MeshioPlusPlusIO` uses (see [Field data](../Utilities/Mesh_Operations.html#field-data)), applied to **both** the source and the target — the usual case is that only the source carries the named data, but a name present on both is exactly what `"on_conflict"` resolves.

## Usage

```json
{
    "modelers" : [{
        "modeler_name" : "KratosMultiphysics.MeshioPlusPlusApplication.modelers.meshio_interpolate_modeler.MeshioInterpolateModeler",
        "Parameters"   : {
            "source_model_part_name" : "Structure",
            "target_model_part_name" : "Target",
            "output_model_part_name" : "StructureOnTarget",
            "operation_settings"     : {
                "nodal_solution_step_data_variables" : ["TEMPERATURE"],
                "method"      : "barycentric",
                "extrapolate" : true
            }
        }
    }]
}
```

## Result

The result — the target's topology plus the interpolated arrays — is written into `"output_model_part_name"`.

## Parameters

| Setting | Default | Description |
|---|---|---|
| `source_model_part_name` | `""` *(required)* | The model part carrying the field data to sample. |
| `target_model_part_name` | `""` *(required)* | The model part contributing the geometry to sample onto. |
| `output_model_part_name` | `""` *(required)* | The model part receiving the result; created when the modeler is instantiated. |
| `operation_settings` | `{}` | Field-data selection plus `"method"` (`"nearest"`/`"barycentric"`), `"names"` (empty = every array the field-data settings collect), `"extrapolate"`, `"default_value"`, `"on_conflict"` (`"error"`/`"overwrite"`/`"suffix"`); forwarded as-is to `MeshioPlusPlusMeshOperations.Interpolate`. |
| `echo_level` | `0` | Verbosity; `> 0` prints the operation's report. |
