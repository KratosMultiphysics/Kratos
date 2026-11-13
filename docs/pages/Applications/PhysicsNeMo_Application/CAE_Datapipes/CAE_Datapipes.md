---
title: CAE Datapipes
keywords: domino transolver cae datapipe stl export
tags: [CAE_Datapipes.md]
sidebar: physicsnemo_application
summary: The CAE exporter writing per-case .npz files in physicsnemo.datapipes.cae's exact layout, the DoMINO and Transolver pipe factories, and in-loop DoMINO deployment with pretrained checkpoints.
---

# CAE datapipes: DoMINO and Transolver on Kratos data

physicsnemo's large aerodynamic models come with config-driven datapipes (`physicsnemo.datapipes.cae`) that read per-case dictionaries. `CaeDatasetExportProcess` writes Kratos results in exactly that layout — one `.npz` per case with the **superset** of the keys `DoMINODataPipe` and `TransolverDataPipe` consume (each pipe selects its subset; extra keys are ignored) — so both pipes run on Kratos data out of the box. FIGConvNet remains bring-your-own-datapipe.

<p align="center">
    <img src="images/npz_layout.svg" alt="The arrays of one CAE case file and which pipelines consume them"/>
</p>
<p align="center">Figure 1: One case file, its arrays and shapes, and the two pipelines that read it.</p>

## The exporter

```json
{
    "python_module" : "cae_dataset_export_process",
    "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.export",
    "Parameters"    : {
        "model_part_name"         : "FluidModelPart",
        "surface_model_part_name" : "FluidModelPart.Body",
        "surface_fields"          : [ { "variable_name" : "PRESSURE",  "data_location" : "node_historical" } ],
        "volume_fields"           : [ { "variable_name" : "VELOCITY",  "data_location" : "node_historical" } ],
        "global_params"           : { "stream_velocity" : 30.0, "air_density" : 1.226 },
        "global_params_order"     : ["stream_velocity", "air_density"],
        "output_path"             : "cae_cases",
        "case_id"                 : 0
    }
}
```

- `model_part_name`: the volume part — its **nodes** become `volume_mesh_centers`, its nodal `volume_fields` the targets.
- `surface_model_part_name` + `surface_source_container` ("Conditions" default): tessellated into the STL through the mesh bridge (`BuildProvenance` — watertight smallest-id triangulation, `tessellation_mode`/`higher_order_mode` knobs apply). Triangle centers/areas/normals are computed from the triangulation; nodal surface fields become per-triangle vertex means, entity/Gauss fields are replicated/mean-collapsed onto the sub-triangles.
- `global_params`: flat dict of doubles. Every key is written both as an individual shape-(1,) array (Transolver's `stream_velocity`/`air_density` convention) and stacked into `global_params_values`/`global_params_reference` (k, 1) for DoMINO. **Kratos Parameters do not preserve insertion order** — set `global_params_order` explicitly for DoMINO, which is order-sensitive. `global_params_reference` (optional dict) defaults to the values.
- `case_id >= 0` names the file `<prefix>_<case_id>.npz` (one file per run, for steady parameter sweeps); the default `-1` uses the STEP.
- Empty field lists omit the corresponding keys (DoMINO then loads in inference mode).
- MPI-aware: topology and fields are gathered (`GatherModelPartToRank0` for the surface, node/field gathers for the volume) and rank 0 writes the serial layout.

## npz layout (the superset)

| Key | Shape | Notes |
|---|---|---|
| `stl_coordinates` | (P, 3) f32 | triangulated surface points |
| `stl_faces` | (3T,) int32 | **flattened** — DoMINO expects it flat |
| `stl_centers`, `surface_mesh_centers` | (T, 3) f32 | triangle centers |
| `stl_areas`, `surface_areas` | (T,) f32 | triangle areas |
| `surface_normals` | (T, 3) f32 | unit normals |
| `surface_fields` | (T, F_s) f32 | when configured |
| `volume_mesh_centers` | (N, 3) f32 | the volume part's nodes |
| `volume_fields` | (N, F_v) f32 | when configured |
| `global_params_values/_reference` | (k, 1) f32 | DoMINO style, `global_params_order` order |
| `<each global key>` | (1,) f32 | Transolver style |
| `TIME`, `STEP` | (1,) | never 0-d (the npz reader rejects 0-d) |

## Loading through the pipes

`torch_dataset` provides lazy factories plus the per-pipe minimal key sets as constants (`DOMINO_SURFACE_KEYS`, `DOMINO_VOLUME_KEYS`, `TRANSOLVER_SURFACE_KEYS`, `TRANSOLVER_VOLUME_KEYS`):

```python
from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import (
    CreateDoMINODataPipe, CreateTransolverDataPipe, CreateCaeDataset)

pipe = CreateDoMINODataPipe(
    "cae_cases", "surface",
    bounding_box=((0, 0, 0), (1, 1, 1)))    # REQUIRED - DoMINO's preprocessing
sample = pipe[0]                            # surface neighbors etc. computed by the pipe (kNN)

pipe = CreateTransolverDataPipe("cae_cases", "volume")
sample = pipe[0]                            # embeddings / fields / fx (globals)
```

Handled gotchas (all verified against physicsnemo 2.2.0):

- DoMINO needs **both** volume and surface bounding boxes even for surface-only models; `bounding_box_surface` defaults to `bounding_box`.
- Surface-neighbor keys are computed internally by kNN — never exported.
- Transolver's default `resolution=200000` breaks small meshes → the factory defaults to `None` (all points); `volume_sample_from_disk` is forced off (the npz reader cannot subsample on disk).
- CPU is the safe default (`device=None`); pass a CUDA device to enable GPU preprocessing.

## DoMINO deployment (`DominoInferenceProcess`)

DoMINO consumes DoMINODataPipe's preprocessed sample dict (geometry encodings, SDF grids, kNN surface neighborhoods), so deployment reuses the exact training pipeline: per execution the current state is exported as a single-case `.npz` into a scratch directory (through this page's export machinery), preprocessed by a `DoMINODataPipe` (`sampling=False`: node/cell order preserved, samples come pre-batched), and the model's predictions are written back — per-node **volume** outputs onto the volume part's nodes, per-triangle **surface** outputs collapsed onto their parent conditions/elements by mean via the mesh-bridge provenance (surface output fields therefore use `"condition"`/`"element"` data locations).

```json
{
    "python_module" : "domino_inference_process",
    "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.inference",
    "Parameters"    : {
        "volume_model_part_name"  : "Main",
        "surface_model_part_name" : "Main.Skin",
        "model_settings"          : { "checkpoint_file" : "domino.mdlus", "checkpoint_type" : "physicsnemo" },
        "model_type"              : "surface",
        "bounding_box"            : [-0.1, -0.1, -0.1, 1.1, 1.1, 1.1],
        "global_params"           : { "stream_velocity" : 30.0, "air_density" : 1.226 },
        "global_params_order"     : ["stream_velocity", "air_density"],
        "datapipe_overrides"      : { "grid_resolution" : [64, 64, 64], "num_surface_neighbors" : 7 },
        "output_fields_surface"   : [ { "variable_name" : "PRESSURE", "data_location" : "condition" } ]
    }
}
```

Consistency requirements (validated the hard way): the **model's** `interp_res` must equal the **datapipe's** `grid_resolution`, and its `num_neighbors_surface` the datapipe's `num_surface_neighbors` — they describe the same background grid and kNN stencils. Both bounding boxes are mandatory for DoMINO's preprocessing (`bounding_box_surface` defaults to `bounding_box`). Pretrained checkpoints (DoMINO-Automotive-Aero) come from NGC/NIM externally; any `.mdlus` DoMINO deploys, model cards included.

### Driving a *pretrained* checkpoint

The public `nvidia/domino_drivaerml` checkpoints are ungated and load through `model_registry` with `"checkpoint_type": "physicsnemo"` — they are exactly the `.mdlus` form this process takes. Two things must be right, and only one of them fails loudly.

**Match the datapipe to the checkpoint.** The shipped `DoMINODataConfig` defaults do not match a pretrained model; the settings below do, and are all expressible today:

```json
"datapipe_overrides" : { "grid_resolution"       : [128, 64, 64],
                         "num_surface_neighbors" : 7,
                         "normalize_coordinates" : true },
"bounding_box"         : [-3.5, -2.25, -0.32,  8.5, 2.25, 3.00],
"bounding_box_surface" : [-1.5, -1.4,  -0.32,  5.0, 1.4,  1.4 ],
"global_params"        : { "stream_velocity" : 30.0, "air_density" : 1.205 },
"global_params_order"  : ["stream_velocity", "air_density"]
```

A wrong `grid_resolution` used to surface as an opaque reshape failure deep inside the geometry encoder; it is now checked against the checkpoint's own `grid_resolution` at load and named. A wrong `normalize_coordinates` or a defaulted `bounding_box_surface` produce **no error at all** — only wrong numbers — so both are warned about. Note `num_surface_neighbors: 7` yields `surface_mesh_neighbors` with K = 6: the kNN includes the point itself and the pipe drops it, exactly as upstream does.

**De-normalize the output.** This is the one that silently corrupts results. A pretrained DoMINO predicts **dimensionless, normalized** fields; the raw tensor is not a pressure. Upstream applies the inverse of the training normalization and then redimensionalizes by `U²ρ`:

```json
"scaling_factors_file" : ".../domino_drivaerml_surface_checkpoint/scaling_factors.pkl",
"normalization"        : "min_max_scaling",
"redimensionalize"     : true
```

For this checkpoint that turns a raw `0.1386` into **−609 Pa** — a Cp of about −1.1 against the 542 Pa dynamic pressure at 30 m/s. Without it the value written into Kratos is wrong by roughly three orders of magnitude *and* shifted. Both settings default to off, so configurations written against raw output are unaffected.

`scaling_factors.pkl` is read through physicsnemo-cfd's restricted unpickler: the file references a `utils.ScalingFactors` class that is not an installed module, so a plain `pickle.load` raises `ModuleNotFoundError`, and the public `ScalingFactors.load` does exactly that. Note also that `global_stats.json`, shipped alongside, is **not** what de-normalization uses — it is informational.

