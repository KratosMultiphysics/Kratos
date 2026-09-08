---
title: Mesh Operations
keywords: meshio meshioplusplus clean transform split refine subdivide agglomerate undo_green decimate decimate_volume smooth remesh remesh_volume optimize_volume partition quality stats diff data_calc data_condition data_manage data_integrate interpolate conservative_interpolate gradient hessian estimate_error voxelize sdf crop_predicate grid
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
- **report-only** (`stats`, `quality`, `data_info`, `data_integrate`): the destination is left untouched and the whole result is in the returned report.

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
| `refine` | Uniform subdivision, or selective/adaptive subdivision of a chosen subset — see below. `record_hierarchy` attaches the persistent `refine:cell_id`/`refine:parent_id` arrays a multigrid caller resolves across the sequence of meshes it keeps, and which `UndoGreen` needs to identify a green group later. | `levels`, `record_parent_ids`, `cells`, `region`, `predicate_array`, `predicate_op`, `predicate_value`, `closure`, `record_levels`, `record_hierarchy` |
| `subdivide` | **Polyhedral** refinement: one polyhedral child per face. Unlike `refine` it handles an arbitrary polyhedron and needs no closure, since a shared face is never touched and no hanging node can appear. Its output is polyhedral, which no Kratos `Element` can hold — `simplexify_result` (default `true`) decomposes it into tetrahedra before the write-back; with it `false` an empty destination is refused by name rather than silently returned. | `record_parent_ids`, `simplexify_result` |
| `agglomerate` | **Polyhedral** coarsening: the many-to-one counterpart of `subdivide`, merging face-adjacent cells into one polyhedron per group. Same `simplexify_result` handling. | `target_group_size`, `simplexify_result` |
| `decimate` | **Surface** decimation by quadric-error edge collapse. | `target_ratio`, `target_faces`, `max_error`, `preserve_boundary`, `preserve_features`, `feature_angle` |
| `decimate_volume` | The **tetrahedral** counterpart of `decimate`: the same quadric-collapse idea, but the budget is a cell count rather than a face count, hence `target_cells`. | `target_ratio`, `target_cells`, `max_error`, `placement`, `preserve_boundary`, `preserve_features`, `feature_angle`, `frozen` |
| `remesh` | **Surface** remeshing of a triangle mesh by ACVD clustering, with an isotropic, quadric or anisotropic metric, curvature gradation and boundary preservation. | `num_clusters`, `subdivide_level`, `subsample_ratio`, `max_subdivide`, `max_iterations`, `max_repair_passes`, `metric`, `gradation`, `preserve_boundary`, `max_anisotropy` |
| `remesh_volume` | A **closed surface in, a tetrahedral volume out**: a lattice cut against the surface's signed distance, so it shares the whole lattice/SDF settings block with `voxelize` and `compute_sdf` rather than declaring its own. | `resolution` or `cell_size`, `bounds`, `padding`, `padding_relative`, `max_cells`, `max_tets`, `warp_fraction`, plus the surface-distance settings |
| `optimize_volume` | Improves a **tetrahedral** mesh without changing its cell budget: vertex relocation plus 2-3/3-2 flips. The worst element's quality is monotone — the report carries it before and after. Uses `optimize_iterations` rather than `max_iterations` because upstream's default here (10) differs from `remesh`'s (100), and one key cannot carry two defaults. | `optimize_iterations`, `relocate`, `flip`, `preserve_boundary`, `min_improvement`, `frozen` |
| `smooth` | Coordinate smoothing (Laplacian, Taubin or ODT). | `method`, `iterations`, `lambda`, `mu`, `fix_boundary`, `preserve_features`, `feature_angle`, `guard_inversion` |
| `reorder` | Node renumbering (`"rcm"`, `"morton"`, `"hilbert"`); the report carries the bandwidth before and after. | `method` |
| `extract_surface` | Boundary extraction (3D → faces, 2D → edges). | `record_parent_ids` |
| `extract_skin` | Boundary skin of a volume mesh (the `SkinDetectionProcess` algorithm). | `linearize` |
| `crop_bbox` | Subset by axis-aligned bounding box. | `box_min`, `box_max`, `keep_partial_cells`, `record_parent_ids` |
| `crop_halfspace` | Subset by half-space `(p - origin)·normal >= 0`. | `origin`, `normal`, `keep_partial_cells`, `record_parent_ids` |
| `crop_predicate` | Keep the cells whose scalar `cell_data` value satisfies a comparison. Shares `refine`'s comparison vocabulary *and its evaluator*, so a non-finite cell value never matches — including under `!=`, where IEEE would say it does. `keep_partial_cells` is deliberately absent: a cell-data predicate is already one value per cell and has nothing to reduce. | `predicate_array`, `predicate_op`, `predicate_value`, `record_parent_ids` |
| `gradient` | Gradient, divergence or curl of a `point_data` field. Green-Gauss is exact for a linear field on any cell. A `cell_data` input raises by name — a piecewise-constant field has no derivative. | `array_name`, `gradient_operator`, `gradient_method`, `location`, `output`, `output_overwrite`, `component` |
| `hessian` | The second-derivative sibling of `gradient`. | `array_name`, `gradient_method`, `location`, `output`, `output_overwrite` |
| `estimate_error` | Zienkiewicz-Zhu recovery-based error indicator, with `"absolute"`, `"fraction"` or `"dorfler"` marking (or `"none"` to estimate only). The marked array is exactly what selective `refine`'s `predicate_array` was built to consume, which closes the adaptive loop. Both outputs are colon-namespaced (`error:zz`, `error:marked`), so point `output`/`marked_name` at registered `Variable` names to read them back — see the write-back note below. | `array_name`, `error_method`, `marking`, `marking_value`, `output`, `marked_name`, `output_overwrite` |
| `voxelize` | A regular hexahedron lattice around the mesh: the whole bounding box (`"all"`), only the cells a surface passes through (`"surface"`), or only those inside it (`"inside"`). The latter two measure against a *surface* — run `extract_skin` first. | `resolution` or `cell_size` (exactly one), `bounds`, `padding`, `padding_relative`, `fill`, `attach_occupancy`, `max_cells`, `output`, plus the surface-distance settings |
| `compute_sdf` | Generates a grid *and* fills it with the signed distance to the source surface. `"octree"` refines adaptively near the surface and reaches the same zero contour as a uniform grid of the same finest resolution from far fewer cells. | `structure`, `resolution` or `cell_size`, `bounds`, `padding`, `padding_relative`, `root_resolution`, `max_depth`, `band_cells`, `record_levels`, `max_cells`, `output`, plus the surface-distance settings |
| `slice` | Planar cross-section (volume → surface, surface → lines). | `origin`, `normal`, `record_parent_ids` |
| `isosurface` | Level set of a scalar `point_data` field, one dimension below the cut cells. Needs `"array_name"` staged as field data (see below), otherwise it throws. | `array_name`, `isovalues`, `record_parent_ids` |
| `attach_quality` | Attaches per-cell quality metrics as cell data, under meshio++'s own array names (`"quality:scaled_jacobian"`, ...). These do not match a registered `Variable`, so the metrics are in the returned report but are **not** written back onto the destination model part — see the write-back note below. | — |
| `data_calc` | Evaluates an expression (array names, numbers, `+-*/`, `abs`/`sqrt`/`min`/`max`/`norm`) into a new array. | `expression`, `output`, `output_overwrite`, `location` |
| `data_condition` | Conditions an array's values in place: `"clamp"` (min/max), `"normalize"` (rescale to `[lo, hi]`), `"standardize"` (zero mean, unit variance). | `mode`, `scope`, `names`, `lo`, `hi`, `output_suffix`, `preserve_dtype`, `location` |
| `data_manage` | Keeps, drops or renames arrays (in that fixed order). | `keep`, `drop`, `rename`, `ignore_missing` |
| `data_info` | *(report-only)* Reports every staged array's name, shape, dtype and value statistics. | — |
| `data_integrate` | *(report-only)* The cell-measure-weighted integral and mean of each array, over the whole mesh and over every named cell region. Per-component throughout, so a vector array reports one entry per component. | `names` |
| `point_data_to_cell_data` | Averages nodal data onto cells (mean over each cell's own nodes). | `names`, `weight`, `prefix`, `output_suffix`, `overwrite`, `nan_policy`, `nan_replacement` |
| `cell_data_to_point_data` | Averages cell data onto nodes (mean, optionally measure-weighted, over the incident cells). | `names`, `weight`, `prefix`, `output_suffix`, `overwrite`, `nan_policy`, `nan_replacement` |
| `split` | *(multi-output)* Partition into submeshes by cell type, connected component, or integer tag. | `split_by`, `tag_name` |
| `partition` | *(multi-output)* Domain decomposition (SFC or KaHIP), with an optional shared-node ghost layer for MPI halos. | `number_of_parts`, `imbalance`, `seed`, `ghost_layers`, `weights_variable`, `record_parent_ids` |
| `stats` | *(report-only)* Bounding box, per-cell-type counts, total area, signed/unsigned volume, inverted-cell count. | — |
| `quality` | *(report-only)* Per-cell quality metrics summarized per metric, plus inverted/degenerate counts. | — |

All operations also read the shared settings `"entity_type"` (`"elements"`, `"conditions"`, or `"automatic"` for both), `"use_deformed_configuration"` (current vs. initial node coordinates), and `"record_provenance"` (default `true`) — whether this operation joins the provenance chain a later `MeshioPlusPlusIO` write records into the file. It is deliberately not the same setting as the IO's own `"provenance"`, which selects a *mode* (`off`/`best_effort`/`required`) for a write; these are different questions.

## Selective refinement (`refine`)

By default `refine` subdivides every cell (byte-identical output when none of the settings below are set). Setting exactly one of the following switches to selective mode; setting more than one is an error:

- `"cells"` — explicit, block-major cell indices to refine.
- `"region"` — the name of an existing sub model part of the source model part (nested parts are addressed `/`-joined, e.g. `"Outer/Inner"`). A part built only from nodes selects every cell touching any of them; one carrying elements/conditions selects those cells directly.
- `"predicate_array"` / `"predicate_op"` (`"<"`, `"<="`, `">"`, `">="`, `"=="`, `"!="`) / `"predicate_value"` — threshold a staged scalar `cell_data` array, e.g. composing with `attach_quality`'s output.

`"closure"` (`"redgreen"` default, `"propagate"`, `"balanced"`) controls how the hanging nodes a partial refinement leaves behind are resolved:

- `"redgreen"` promotes an affected neighbor's split to the smallest conforming superset — keeps the extra refinement local.
- `"propagate"` promotes any affected neighbor straight to a full split — always conforming, but not local (converges to uniform refinement of the whole connected component).
- `"balanced"` does not close at all: it keeps genuine hanging nodes and only enforces 2:1 balance, the classic adaptive-mesh-refinement meaning of "propagate" — the only mode whose cost is bounded by the selection rather than by the mesh. It is materially better than it once was: meshio++ v9.23.0 fixed a second balanced pass leaving the two sides of a hanging interface referencing distinct but exactly coincident nodes (a *torn* mesh, not merely a 1-irregular one), and fixed `refine:hanging` under-reporting the constrained nodes it promises to mark exactly.

`"record_levels"` attaches meshio++'s own `Int64` `refine:level` bookkeeping array (cumulative refinement depth per cell); `"balanced"` always attaches an `Int64` `refine:hanging` array flagging constrained nodes. Both are computed but, like `record_parent_ids`'s `*:parent_cell` arrays on every operation that has it, are **not** written back onto the destination model part through `Execute()` — their names are not registered Kratos `Variable`s (see the write-back constraint below).

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
| `write_ids` | Also stage Kratos node/entity ids as `KRATOS_*_ID` arrays, for inspection. |
| `write_mdpa_ids` | Stage those ids as meshio++'s `"mdpa:id"` instead — a *format contract*, not a diagnostic: the `mdpa` writer renumbers entities to `1..n` unless the mesh carries it, so a model part whose ids have gaps (routine after a sub model part extraction or an entity removal) loses them silently on write without this. |

All default to empty (nothing staged) — an operation that needs a specific array, like `isosurface`'s `"array_name"`, has that array present in the mesh only when its variable is listed in one of the settings above.

## Surface distance and grids

`voxelize` and `compute_sdf` share one block of surface-distance settings, because meshio++ embeds the same options struct in both:

| Setting | Default | Description |
|---|---|---|
| `sdf_sign` | `"pseudonormal"` | How the distance is signed. |
| `sdf_weight` | `"angle"` | How incident faces are weighted at a vertex (`pseudonormal` only). |
| `sdf_location` | `"corner"` | Evaluate at cell corners or centres. |
| `band` | `0.0` | Distances beyond this are clamped and marked; non-positive computes the full field. |
| `record_closest_cell` / `record_inside` | `false` | Attach the extra diagnostic arrays. |
| `watertight_check` | `"warn"` | What to do about a surface that is not closed. |
| `grid_cell_size` | `0.0` | Accelerator bucket size; `0` derives one. It cannot change the answer — only the time. |
| `max_winding_work` | `2.0e9` | Refuses a winding-number query above this rather than running for an hour. |

> **`"padding_relative"` carries one default for operations upstream gives two.** meshio++'s own default is `0.0` for `voxelize` but `0.1` for `compute_sdf` and `remesh_volume`; a single Kratos setting can only have one, and it carries `voxelize`'s `0.0`. So the latter two are padded *less* here than a direct meshio++ call unless you set it explicitly. It is kept at `0.0` rather than "corrected" because changing it would silently change what every existing `voxelize` call produces.

> **`"output"` is what makes these usable from Kratos.** meshio++ names its result `sdf:distance` (and `voxel:occupancy`), which no Kratos `Variable` can ever be — they are colon-namespaced and Kratos variables are not. Chaining a `data_manage` rename afterwards does not help either, because the array is already gone at the model part boundary. Setting `"output"` renames the array *before* the write-back, so pointing it at `DISTANCE` gives a level-set initialisation in one call.

> **The write-back constraint.** A resulting array is written back onto the destination model part only when it is `Float64` or `Int64` **and** its name matches a registered `Variable<T>` with the right component count — Kratos stores non-historical/historical data keyed by `Variable` objects, not arbitrary strings. An operation's own invented names, such as `attach_quality`'s `"quality:scaled_jacobian"` or `refine`'s `"refine:level"`/`"refine:hanging"`, are computed and appear in the report but cannot be retrieved from the model part afterwards, regardless of dtype. Point an operation's own naming setting (`data_calc`'s `"output"`, `data_manage`'s `"rename"`, ...) at an existing variable name to get the result back onto the mesh.

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
| `ConservativeInterpolate(source, target, settings, destination)` | The conservative sibling of `Interpolate`: redistributes by intersected cell measure so the **summed quantity is preserved**, rather than sampling point values. The right choice for an extensive field (a mass, a heat load) and the wrong one for an intensive one (a temperature). Settings: `"names"`, `"default_value"`, `"on_conflict"`. |
| `Diff(first, second, settings)` | Compares two model parts with `"absolute_tolerance"`/`"relative_tolerance"`; returns `{"equal": bool}`. |
| `ComputeStatistics(model_part)` | The `"stats"` report, callable directly without a destination. |
| `ComputeQuality(model_part)` | The `"quality"` report, callable directly without a destination. |
| `ComputeBandwidth(model_part)` | The bandwidth of the node adjacency graph — what `"reorder"` reduces. |
| `Grid(settings, destination)` | Builds a regular hexahedron lattice from **no source model part at all** — the library's only generator. Settings: `dims`, `origin`, `spacing`, `max_cells`. |
| `DistanceToSurface(query, surface, settings, destination)` | Attaches the signed distance from `query` to `surface`. Two independent meshes, like `Interpolate`; honours `"output"` exactly as `compute_sdf` does. |
| `CheckSurfaceWatertight(surface)` | *(report-only)* `watertight` plus the boundary-edge, non-manifold-edge, inconsistent-pair and degenerate-triangle counts behind the verdict. |

meshio++'s `undo_green` is deliberately **not** exposed either, for a different reason: it resolves a refinement's green closure through the `refine:cell_id`/`refine:parent_id` arrays that `refine(record_hierarchy)` attaches, and those names are colon-namespaced — so by the write-back constraint above they can never reach a model part. The hierarchy is already gone by the time a refined mesh *is* a `ModelPart`, which makes a ModelPart-taking wrapper for it impossible rather than merely awkward. It stays usable from meshio++ itself, where the mesh never leaves the library.

meshio++'s `sample_distance` is deliberately **not** exposed: it takes a raw `(n, 3)` point array rather than a model part, and `DistanceToSurface` covers the same ground for any mesh.

## Relation to existing Kratos utilities

Where Kratos already has an equivalent (`SkinDetectionProcess`, `ReorderAndOptimizeModelPartProcess`, `IntegrationValuesExtrapolationToNodesProcess`), the meshio++ operations are a format-interop alternative, not a replacement — use them when you are already moving a mesh through meshio++ (reading/writing a non-native format, converting between formats) and want the operation applied in the same pass.
