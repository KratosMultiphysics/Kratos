---
title: Mesh Operations
keywords: meshio meshioplusplus clean transform split refine decimate smooth partition quality stats diff data_calc data_condition data_manage interpolate
tags: [Mesh_Operations.md]
sidebar: meshioplusplus_application
summary: The meshio++ mesh and data operations exposed as a Kratos utility.
---

## Overview

`MeshioPlusPlusMeshOperations` (`custom_utilities/meshioplusplus_mesh_operations.h`) exposes the meshio++ operations layer to Kratos. Every operation is data-driven through `Parameters`, keyed by an `"operation"` name that mirrors the meshio++ command line verbs, so a new upstream operation is a table entry rather than a new class.

It is reachable three ways:

- directly from C++, `MeshioPlusPlusMeshOperations::Execute(...)`;
- directly from Python, `KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(...)`;
- data-driven from JSON, through the [Meshio Operation Modeler](../Modelers/Meshio_Operation_Modeler.html).

Operations fall into two groups:

- **mesh-producing** (`clean`, `transform`, `refine`, ...): the result is written into the destination model part, and the returned report carries the operation's own counts (welded points, bandwidth before/after, ...);
- **report-only** (`stats`, `quality`): the destination is left untouched and the whole result is in the returned report.

`"split"` and `"partition"` are a third case: **multi-output**. They create one model part per piece — named `<destination>_<operation>_<index>` and registered as siblings in the same `Model` — because a Kratos sub model part is a view into its root's entity containers while each piece is independently renumbered from 1, so the pieces cannot be sub model parts of a shared destination.

```python
import KratosMultiphysics as KM
import KratosMultiphysics.MeshioPlusPlusApplication as KratosMeshioPlusPlus

settings = KM.Parameters("""{
    "operation" : "clean",
    "weld"      : true,
    "tolerance" : 1e-6
}""")
destination = model.CreateModelPart("Clean")
report = KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.Execute(source, settings, destination)
```

> meshio++ is serial: these operations do not support distributed model parts. The intended distributed workflow is `"partition"` with ghost layers, feeding an MPI assembly.

## Operations

Every setting below is read from the same flat `Parameters` object as `"operation"` itself (directly, or under `"operation_settings"` when driven by the [Operation Modeler](../Modelers/Meshio_Operation_Modeler.html)).

| Operation | Description | Settings used |
|---|---|---|
| `clean` | Weld coincident nodes, drop degenerate/duplicate cells, remove orphan nodes. | `weld`, `tolerance`, `remove_orphans`, `drop_degenerate`, `drop_duplicate_cells` |
| `transform` | Affine transform (scale, then rotate, then translate), optionally rotating vector field data. | `scale`, `rotation_axis`, `rotation_angle`, `translation`, `rotate_vector_data` |
| `convert_cells` | Element-type conversion: `"linearize"` (drop mid nodes), `"simplexify"` (decompose into simplices), `"elevate"` (linear → quadratic). | `mode`, `record_parent_ids` |
| `refine` | Uniform subdivision. | `levels`, `record_parent_ids` |
| `decimate` | Surface decimation by quadric-error edge collapse. | `target_ratio`, `target_faces`, `max_error`, `preserve_boundary`, `preserve_features`, `feature_angle` |
| `smooth` | Coordinate smoothing (Laplacian/Taubin). | `method`, `iterations`, `lambda`, `mu`, `fix_boundary`, `preserve_features`, `feature_angle`, `guard_inversion` |
| `reorder` | Node renumbering (`"rcm"`, `"morton"`, `"hilbert"`); the report carries the bandwidth before and after. | `method` |
| `extract_surface` | Boundary extraction (3D → faces, 2D → edges). | `record_parent_ids` |
| `extract_skin` | Boundary skin of a volume mesh (the `SkinDetectionProcess` algorithm). | `linearize` |
| `crop_bbox` | Subset by axis-aligned bounding box. | `box_min`, `box_max`, `keep_partial_cells`, `record_parent_ids` |
| `crop_halfspace` | Subset by half-space `(p - origin)·normal >= 0`. | `origin`, `normal`, `keep_partial_cells`, `record_parent_ids` |
| `slice` | Planar cross-section (volume → surface, surface → lines). | `origin`, `normal`, `record_parent_ids` |
| `isosurface` | Level set of a scalar `point_data` field, one dimension below the cut cells. Needs `"array_name"` staged as field data (see below), otherwise it throws. | `array_name`, `isovalues`, `record_parent_ids` |
| `attach_quality` | Attaches per-cell quality metrics as cell data, under meshio++'s own array names (`"quality:scaled_jacobian"`, ...). These do not match a registered `Variable`, so the metrics are in the returned report but are **not** written back onto the destination model part — see the write-back note below. | — |
| `data_calc` | Evaluates an expression (array names, numbers, `+-*/`, `abs`/`sqrt`/`min`/`max`/`norm`) into a new array. | `expression`, `output`, `output_overwrite`, `location` |
| `data_condition` | Conditions an array's values in place: `"clamp"` (min/max), `"normalize"` (rescale to `[lo, hi]`), `"standardize"` (zero mean, unit variance). | `mode`, `scope`, `names`, `lo`, `hi`, `output_suffix`, `preserve_dtype`, `location` |
| `data_manage` | Keeps, drops or renames arrays (in that fixed order). | `keep`, `drop`, `rename`, `ignore_missing` |
| `data_info` | *(report-only)* Reports every staged array's name, shape, dtype and value statistics. | — |
| `point_data_to_cell_data` | Averages nodal data onto cells (mean over each cell's own nodes). | `names`, `weight`, `prefix`, `output_suffix`, `overwrite`, `nan_policy`, `nan_replacement` |
| `cell_data_to_point_data` | Averages cell data onto nodes (mean, optionally measure-weighted, over the incident cells). | `names`, `weight`, `prefix`, `output_suffix`, `overwrite`, `nan_policy`, `nan_replacement` |
| `split` | *(multi-output)* Partition into submeshes by cell type, connected component, or integer tag. | `split_by`, `tag_name` |
| `partition` | *(multi-output)* Domain decomposition (SFC or KaHIP), with an optional shared-node ghost layer for MPI halos. | `number_of_parts`, `imbalance`, `seed`, `ghost_layers`, `weights_variable`, `record_parent_ids` |
| `stats` | *(report-only)* Bounding box, per-cell-type counts, total area, signed/unsigned volume, inverted-cell count. | — |
| `quality` | *(report-only)* Per-cell quality metrics summarized per metric, plus inverted/degenerate counts. | — |

All operations also read the shared settings `"entity_type"` (`"elements"`, `"conditions"`, or `"automatic"` for both) and `"use_deformed_configuration"` (current vs. initial node coordinates).

## Field data

Nodal, elemental and conditional data reaches every operation through the same settings `MeshioPlusPlusIO` uses, read once per call and staged into the meshio++ mesh before the operation runs:

| Setting | Description |
|---|---|
| `nodal_solution_step_data_variables` | Historical nodal variables to stage as point data. |
| `nodal_data_value_variables` | Non-historical nodal variables to stage as point data. |
| `nodal_flags` | Node `Flags` to stage as point data. |
| `element_data_value_variables` | Non-historical elemental variables to stage as cell data. |
| `element_flags` | Element `Flags` to stage as cell data. |
| `condition_data_value_variables` | Non-historical conditional variables to stage as cell data. |
| `condition_flags` | Condition `Flags` to stage as cell data. |
| `gauss_point_variables_in_elements` | Elemental Gauss-point values to stage as cell data. |
| `write_ids` | Also stage Kratos node/entity ids as an array. |

All default to empty (nothing staged) — an operation that needs a specific array, like `isosurface`'s `"array_name"`, has that array present in the mesh only when its variable is listed in one of the settings above.

> **The write-back constraint.** A resulting array is written back onto the destination model part only when its name matches a registered `Variable<T>` with the right component count — Kratos stores non-historical/historical data keyed by `Variable` objects, not arbitrary strings. An operation's own invented names, such as `attach_quality`'s `"quality:scaled_jacobian"`, are computed and appear in the report but cannot be retrieved from the model part afterwards. Point an operation's own naming setting (`data_calc`'s `"output"`, `data_manage`'s `"rename"`, ...) at an existing variable name to get the result back onto the mesh.

Query the exact set a build supports:

```python
KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.GetSupportedOperations()
```

## Additional entry points

Beyond `Execute`, `MeshioPlusPlusMeshOperations` exposes operations that do not fit the single-source, single-`"operation"`-setting shape above:

| Method | Description |
|---|---|
| `Merge(sources, settings, destination)` | Merges several model parts into one, with `"weld"`/`"tolerance"` (coincident-node welding), `"drop_duplicate_cells"`. |
| `Interpolate(source, target, settings, destination)` | Samples `source`'s field data onto `target`'s geometry (`"method"`: `"nearest"`/`"barycentric"`; `"names"`, `"extrapolate"`, `"default_value"`, `"on_conflict"`), filling `destination` with the target's topology plus the interpolated arrays. Needs two independent meshes, so it is not reachable through `Execute` — use it directly, or through the dedicated [Meshio Interpolate Modeler](../Modelers/Meshio_Interpolate_Modeler.html). |
| `Diff(first, second, settings)` | Compares two model parts with `"absolute_tolerance"`/`"relative_tolerance"`; returns `{"equal": bool}`. |
| `ComputeStatistics(model_part)` | The `"stats"` report, callable directly without a destination. |
| `ComputeQuality(model_part)` | The `"quality"` report, callable directly without a destination. |
| `ComputeBandwidth(model_part)` | The bandwidth of the node adjacency graph — what `"reorder"` reduces. |

## Relation to existing Kratos utilities

Where Kratos already has an equivalent (`SkinDetectionProcess`, `ReorderAndOptimizeModelPartProcess`, `IntegrationValuesExtrapolationToNodesProcess`), the meshio++ operations are a format-interop alternative, not a replacement — use them when you are already moving a mesh through meshio++ (reading/writing a non-native format, converting between formats) and want the operation applied in the same pass.
