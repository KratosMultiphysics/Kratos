---
title: Mesh Operations
keywords: meshio meshioplusplus clean transform split refine decimate smooth partition quality stats diff
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
| `isosurface` | Level set of a scalar `point_data` field, one dimension below the cut cells. | `array_name`, `isovalues`, `record_parent_ids` |
| `attach_quality` | Attaches per-cell quality metrics as cell data. | — |
| `split` | *(multi-output)* Partition into submeshes by cell type, connected component, or integer tag. | `split_by`, `tag_name` |
| `partition` | *(multi-output)* Domain decomposition (SFC or KaHIP), with an optional shared-node ghost layer for MPI halos. | `number_of_parts`, `imbalance`, `seed`, `ghost_layers`, `weights_variable`, `record_parent_ids` |
| `stats` | *(report-only)* Bounding box, per-cell-type counts, total area, signed/unsigned volume, inverted-cell count. | — |
| `quality` | *(report-only)* Per-cell quality metrics summarized per metric, plus inverted/degenerate counts. | — |

All operations also read the shared settings `"entity_type"` (`"elements"`, `"conditions"`, or `"automatic"` for both) and `"use_deformed_configuration"` (current vs. initial node coordinates).

Query the exact set a build supports:

```python
KratosMeshioPlusPlus.MeshioPlusPlusMeshOperations.GetSupportedOperations()
```

## Additional entry points

Beyond `Execute`, `MeshioPlusPlusMeshOperations` exposes four operations that do not fit the single-source, single-`"operation"`-setting shape above:

| Method | Description |
|---|---|
| `Merge(sources, settings, destination)` | Merges several model parts into one, with `"weld"`/`"tolerance"` (coincident-node welding), `"drop_duplicate_cells"`. |
| `Diff(first, second, settings)` | Compares two model parts with `"absolute_tolerance"`/`"relative_tolerance"`; returns `{"equal": bool}`. |
| `ComputeStatistics(model_part)` | The `"stats"` report, callable directly without a destination. |
| `ComputeQuality(model_part)` | The `"quality"` report, callable directly without a destination. |
| `ComputeBandwidth(model_part)` | The bandwidth of the node adjacency graph — what `"reorder"` reduces. |

## Relation to existing Kratos utilities

Where Kratos already has an equivalent (`SkinDetectionProcess`, `ReorderAndOptimizeModelPartProcess`, `IntegrationValuesExtrapolationToNodesProcess`), the meshio++ operations are a format-interop alternative, not a replacement — use them when you are already moving a mesh through meshio++ (reading/writing a non-native format, converting between formats) and want the operation applied in the same pass.
