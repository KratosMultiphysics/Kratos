---
title: Meshio Input Modeler
keywords: meshio meshioplusplus modeler input mdpa vtu gmsh med
tags: [Meshio_Input_Modeler.md]
sidebar: meshioplusplus_application
summary: Imports a model part from any meshio++-supported mesh format.
---

## Overview

`MeshioInputModeler` (`modelers.meshio_input_modeler.MeshioInputModeler`) imports a model part from any format `MeshioPlusPlusIO` supports (see the [Overview](../General/Overview.html#supported-formats)).

The format is taken from `"input_format"`, or resolved from the extension of `"input_filename"` when it is `"auto"`.

## Usage

```json
{
    "modelers" : [{
        "modeler_name" : "KratosMultiphysics.MeshioPlusPlusApplication.modelers.meshio_input_modeler.MeshioInputModeler",
        "Parameters"   : {
            "model_part_name" : "Structure",
            "input_filename"  : "bracket.msh"
        }
    }]
}
```

## Parameters

| Setting | Default | Description |
|---|---|---|
| `model_part_name` | `""` *(required)* | The model part to create and fill; it is created when the modeler is instantiated, before the solver adds its variables. |
| `input_filename` | `""` *(required)* | The file to read. |
| `input_format` | `"auto"` | Explicit format name (see [`GetSupportedReadFormats`](../General/Overview.html)), or `"auto"` to resolve it from the extension. |
| `time_step` | `0` | Selects one step of a multi-step file for the formats meshio++ reads selectively (mdpa, med, exodus, gmsh, tecplot, ensight, cgns, openfoam); negative counts from the end. Any other format ignores it and is always read whole. |
| `lenient` | `false` | Downgrades an mdpa/med construct the reader cannot represent from an error to a warning and a skip, instead of throwing. |
| `openfoam_region` | `""` | OPENFOAM only: selects one region of a case with no single `constant/polyMesh`, but `constant/<region>/polyMesh` per region. |
| `echo_level` | `0` | Verbosity. |
