---
title: Process reference
keywords: process reference settings json factory export inference validation remesh adjoint
tags: [Process_Reference.md]
sidebar: physicsnemo_application
summary: Every process the application ships - the JSON that attaches it, every setting with its default and meaning - grouped by what the process does to the solve.
---

# Process reference

Everything that can be attached to a solve from `ProjectParameters.json` lives under `processes/`, and nothing else in the application has a `Factory`. This page lists all of them with their settings. Defaults are quoted from the code; a value of `PLEASE_SPECIFY_...` means the key is mandatory.

<p align="center">
    <img src="images/execution_points.svg" alt="Where each family of processes fires inside the AnalysisStage"/>
</p>
<p align="center">Figure 1: When each process runs. See <a href="Architecture.html">Architecture</a> for the reasoning.</p>

## Which process

| One sample looks like | Export it with | Run a model on it with |
|---|---|---|
| nodal / elemental / condition fields as flat vectors | `dataset_export_process`, `streaming_dataset_export_process` | `inference_process`, `hybrid_initialization_process`, `onnx_inference_process`, `triton_inference_process` |
| the same over time | `dataset_export_process` per step | `time_series_inference_process` |
| a regular voxel grid | `grid_dataset_export_process` | `grid_inference_process`, `superresolution_process`, `sequence_inference_process`, `diffusion_inference_process` |
| the mesh with its connectivity | `mesh_export_process` (`.pmsh`), `curator_export_process` | `graph_inference_process` |
| the nodes as an unordered cloud | `dataset_export_process` | `point_cloud_inference_process` |
| a CAD surface plus a volume | `cae_dataset_export_process` | `domino_inference_process`, `nim_inference_process` |
| particles with trajectories | `dataset_export_process` per step | `particle_inference_process` |
| a few parameters, a whole field back | - | `rom_surrogate_process` |
| no model, only the PDE and its boundary data | - | `pinn_solve_process` |
| a scrubbable digital twin | `usd_export_process` | - |
| judging, adapting, differentiating | - | `validation_metrics_process`, `adaptive_remesh_process`, `adjoint_sensitivity_process` |

## Blocks shared by several processes

**A field entry.** Wherever a list of fields appears, each entry names a variable and a [data location](Kratos_Concepts_For_ML.html):

```json
{ "variable_name" : "VELOCITY", "data_location" : "node_historical" }
```

Locations: `node_historical`, `node_non_historical`, `element`, `condition`, `element_gauss_point`, `condition_gauss_point` (the last two read-only).

**`model_settings`** on every torch-backed inference process, validated by `model_registry.LoadModel`:

| Key | Default | Meaning |
|---|---|---|
| `checkpoint_file` | mandatory | path to the checkpoint; the model card is looked up as `<checkpoint_file>.card.json` |
| `checkpoint_type` | `"torchscript"` | `"torchscript"` (`torch.jit.load`) or `"physicsnemo"` (`Module.from_checkpoint`, `.mdlus`) |
| `device` | `"auto"` | `"auto"` (CUDA when available), `"cpu"`, `"cuda"`, `"cuda:N"` |
| `model_card_policy` | `"advisory"` | what a field-list mismatch with the card does: `"advisory"` warns, `"strict"` raises, `"ignore"` skips the check (de-normalization is applied regardless) |
| `torch_compile` | `false` | wrap a physicsnemo checkpoint in `torch.compile(fullgraph=True)`; refused for TorchScript |
| `nvtx_ranges` | `false` | NVTX ranges around gather, forward and scatter for Nsight Systems |

**`uncertainty`** on `inference_process` and everything inheriting it (see [Uncertainty](../Uncertainty/Uncertainty.html)):

| Key | Default | Meaning |
|---|---|---|
| `method` | `"none"` | `"mc_dropout"`, `"ensemble"` or `"gp"` |
| `num_samples` | `16` | stochastic passes (dropout) |
| `seed` | `-1` | RNG seed for the passes; `-1` leaves it alone |
| `gp_head_file` | `""` | the `.gp_head.pt` sidecar for `"gp"` |
| `gp_feature_fields` | `[]` | features the GP head reads instead of the model inputs |
| `retain_ensemble` | `false` | keep the `(M, ...)` member stack as `last_ensemble` for CRPS |
| `uncertainty_fields` | `[]` | one field entry per output field; receives the standard deviation |

**`ood_guard`** on the same processes: `{"guard_file": "<checkpoint>.ood_guard", "policy": "advisory"}` - the guard written by `TrainModel`, and what a flagged input does (`"advisory"` warns, `"strict"` raises, `"ignore"` disables). Inputs are normalized per the card *before* the check. See [Uncertainty](../Uncertainty/Uncertainty.html).

**`geometry_guard`** on the same processes: `{"guard_file": "<checkpoint>.geometry_guard.npz", "policy": "advisory"}` - the same policy vocabulary asked of the *shape* instead of the field values. The model part's boundary surface is scored against a density model of the training family, checked once at the first execution and again whenever the node count changes. See [Uncertainty](../Uncertainty/Uncertainty.html).

**`execution_point`** (`"finalize_solution_step"` by default, or `"initialize_solution_step"`) and **`output_interval`** (`1`) appear on every per-step process.

**Tessellation keys** on the mesh-based exporters and DoMINO: `tessellation_mode` (`"smallest_id_diagonal"` or the legacy `"fan"`), `higher_order_mode` (`"reduce"`, `"subdivide"`, `"curved"`), `curved_refinement_levels` (`2`), `source_container` (`"Elements"` or `"Conditions"`). See [Mesh Bridge](../Mesh_Bridge/Mesh_Bridge.html).

## Export processes

Package `KratosMultiphysics.PhysicsNeMoApplication.processes.export`. All write in `ExecuteFinalizeSolutionStep`, every `output_interval` steps, and need no torch.

### `dataset_export_process`

One `.npz` per sampled step with one array per field (keys `"<VARIABLE>__<location>"`) plus `TIME` and `STEP`. MPI-aware: ghost-free gathers, rank 0 writes the serial layout.

```json
{ "python_module" : "dataset_export_process",
  "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.export",
  "Parameters"    : { "model_part_name" : "ThermalModelPart",
                      "list_of_fields"  : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                      "output_path"     : "physics_nemo_dataset", "file_prefix" : "sample", "output_interval" : 1 } }
```

| Key | Default | Meaning |
|---|---|---|
| `model_part_name` | mandatory | the part to read |
| `list_of_fields` | one placeholder entry | field entries |
| `output_path` | `"physics_nemo_dataset"` | folder |
| `file_prefix` | `"sample"` | `<prefix>_<step>.npz` |
| `output_interval` | `1` | |

### `grid_dataset_export_process`

The same fields resampled onto a regular voxel grid, `(C, D, H, W)` per file - the input of the grid, sequence and diffusion processes.

```json
{ "python_module" : "grid_dataset_export_process",
  "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.export",
  "Parameters"    : { "model_part_name" : "ThermalModelPart", "grid_shape" : [32, 32, 2], "bounding_box" : [],
                      "list_of_fields"  : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ] } }
```

| Key | Default | Meaning |
|---|---|---|
| `grid_shape` | `[8, 8, 8]` | voxels per axis; use a thin axis of 2 for planar cases and squeeze it at deployment |
| `bounding_box` | `[]` | `[xmin, ymin, zmin, xmax, ymax, zmax]`; empty means the model part's own box |
| `output_path`, `file_prefix`, `output_interval` | `"physics_nemo_grid_dataset"`, `"grid"`, `1` | |

### `mesh_export_process`

The tessellated model part with its fields as a `.pmsh` series - the layout physicsnemo's `MeshReader` consumes. MPI-aware through a rank-0 shadow part.

```json
{ "python_module" : "mesh_export_process",
  "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.export",
  "Parameters"    : { "model_part_name" : "FluidModelPart", "output_path" : "physics_nemo_meshes", "zero_pad_steps" : 4,
                      "list_of_fields"  : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ] } }
```

| Key | Default | Meaning |
|---|---|---|
| `source_container` | `"Elements"` | or `"Conditions"` for a surface series |
| tessellation keys | see above | |
| `output_path`, `file_prefix`, `output_interval` | `"physics_nemo_meshes"`, `"mesh"`, `1` | files are `<prefix>_<step>.pmsh` - the suffix is the reader's glob |
| `zero_pad_steps` | `0` | pad the step number; the reader sorts lexicographically, so `mesh_10` sorts before `mesh_2` without it |

### `cae_dataset_export_process`

Per-case `.npz` in the exact layout of `physicsnemo.datapipes.cae` (triangulated STL surface, surface and volume fields, ordered global parameters) for DoMINO and Transolver pipelines.

```json
{ "python_module" : "cae_dataset_export_process",
  "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.export",
  "Parameters"    : { "model_part_name" : "FluidModelPart", "surface_model_part_name" : "FluidModelPart.Body",
                      "surface_fields" : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ],
                      "volume_fields"  : [ { "variable_name" : "VELOCITY", "data_location" : "node_historical" } ],
                      "global_params_order" : ["STREAM_VELOCITY", "AIR_DENSITY"], "case_id" : 7 } }
```

| Key | Default | Meaning |
|---|---|---|
| `surface_model_part_name` | mandatory | the skin |
| `surface_source_container` | `"Conditions"` | |
| `surface_fields`, `volume_fields` | `[]` | field entries |
| `global_params_order` | `[]` | explicit order of the global parameters - `Parameters.keys()` is alphabetical, so the order must be written down |
| tessellation keys | see above | |
| `output_path`, `file_prefix` | `"physics_nemo_cae_dataset"`, `"case"` | |
| `case_id` | `-1` | the case index in the file name; `-1` uses the step |

### `curator_export_process`

A tessellated mesh series as an AI-ready Zarr store or a VTU series. Zarr is written by physicsnemo itself by default; only `"zarr_backend" : "curator"` or `"sink" : "vtu"` needs the git-only `physicsnemo_curator`.

| Key | Default | Meaning |
|---|---|---|
| `sink` | `"zarr"` | or `"vtu"` |
| `zarr_backend` | `"physicsnemo"` | `"physicsnemo"` (`mesh.io.to_zarr`, one store per step) or `"curator"` |
| `compression_level`, `chunk_size_mb` | `3`, `1.0` | Zarr chunking |
| `source_container`, tessellation keys, `output_path` (`"physics_nemo_curated"`), `file_prefix` (`"mesh"`), `output_interval` | | |

### `usd_export_process`

A time-sampled OpenUSD stage of the running solve, predicted and uncertainty fields included, readable by usdview, Omniverse or Blender. Needs `usd-core`.

| Key | Default | Meaning |
|---|---|---|
| `kind` | `"mesh"` | or `"points"` for a particle cloud (torch-free) |
| `output_file` | `"physics_nemo_twin.usda"` | recreated if it exists |
| `prim_path` | `""` | root prim; empty derives it from the model part name |
| `up_axis`, `meters_per_unit` | `"Z"`, `1.0` | stage metadata |
| `time_source` | `"step"` | or `"time"` - what becomes the USD time code |
| `time_codes_per_second` | `1.0` | |
| `source_container`, tessellation keys, `output_interval` | | topology is re-sampled only on steps where it changed |

### `streaming_dataset_export_process`

No files: each sampled step is pushed into a `LiveSampleQueue` that a training loop drains while the solve is still running. Samples are byte-identical to `dataset_export_process`'s.

| Key | Default | Meaning |
|---|---|---|
| `max_queue_size` | `0` | bound on the queue (`0` unbounded); the producer blocks when full |
| `close_on_finalize` | `true` | close the queue in `ExecuteFinalize` so the consumer learns the solve ended |

## Inference processes

Package `KratosMultiphysics.PhysicsNeMoApplication.processes.inference`. Every one gathers Kratos data into tensors, runs a model, and writes the prediction onto Kratos variables, de-normalized by the model card when one exists.

### `inference_process`

The base contract: input fields concatenated into one `(n_entities, width)` tensor, one forward pass, the output split back into the output fields. Everything else on this page that takes a trained model follows the same gather/forward/scatter shape.

```json
{ "python_module" : "inference_process",
  "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.inference",
  "Parameters"    : { "model_part_name" : "ThermalModelPart",
                      "model_settings"  : { "checkpoint_file" : "surrogate.pt", "device" : "auto" },
                      "input_fields"    : [ { "variable_name" : "HEAT_FLUX",   "data_location" : "node_historical" } ],
                      "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
                      "execution_point" : "finalize_solution_step" } }
```

| Key | Default | Meaning |
|---|---|---|
| `model_part_name`, `model_settings` | mandatory | |
| `input_fields`, `output_fields` | placeholders | entities are those of the first location; outputs default to `node_non_historical` |
| `execution_point`, `output_interval` | `"finalize_solution_step"`, `1` | |
| `ood_guard`, `geometry_guard`, `uncertainty` | `{}` | the shared blocks |

### `hybrid_initialization_process`

`inference_process` run once, in `ExecuteBeforeSolutionLoop`, to seed the solver instead of replacing it: Newton starts from the prediction and keeps its convergence guarantee. Same settings, except `execution_point` and `output_interval` are refused.

### `onnx_inference_process`

The `inference_process` contract with the forward pass through a cached ONNX Runtime session; torch stays only as the array bridge. `model_settings` has its own schema:

| Key | Default | Meaning |
|---|---|---|
| `onnx_file` | mandatory | the card is `<onnx_file>.card.json` |
| `device` | `"cpu"` | `"cuda"` or `"cuda:N"` needs `onnxruntime-gpu` |
| `require_device` | `false` | turn ONNX Runtime's silent CPU fallback into an error |
| `model_card_policy` | `"advisory"` | |

### `triton_inference_process`

The `inference_process` contract against a running Triton Inference Server; the solver host needs neither weights nor a GPU. `model_settings`:

| Key | Default | Meaning |
|---|---|---|
| `url`, `protocol` | `"localhost:8000"`, `"http"` | or `"grpc"` |
| `model_name`, `model_version` | `"kratos_surrogate"`, `""` | as in the repository written by `deployment.triton_export` |
| `input_name`, `output_name` | `"input"`, `"output"` | tensor names in `config.pbtxt` |
| `timeout` | `0.0` | |
| `card_file`, `model_card_policy` | `""`, `"advisory"` | the card is not next to a checkpoint here, so it is named explicitly |

### `point_cloud_inference_process`

The nodes as an unordered cloud: coordinates plus features in, per-node fields out, for Transolver, GeoTransolver, FLARE, FIGConvNet and generic point models - plus two OPERATOR interfaces that map a whole case to a field instead, xDeepONet from case parameters and GLOBE from boundary data. `inference_process` settings plus:

| Key | Default | Meaning |
|---|---|---|
| `model_interface` | `"generic"` | `"generic"` (coordinates prepended to the features, a batch axis added), `"transolver"`, `"geotransolver"`, `"flare"`, `"figconvnet"`, `"deeponet"`, `"globe"` |
| `normalize_coordinates` | `true` | min-max per model part |
| `pass_geometry` | `true` | hand the geometry tensor to models built with a geometry input |
| `subsampling` | `{}` | a token budget - keys below |
| `branch_input` | `{}` | `"deeponet"` only - the case parameters, keys below |
| `trunk_dimension` | `0` | `"deeponet"`/`"globe"` only - how many coordinate columns the model takes; `0` follows `DOMAIN_SIZE` |
| `globe` | `{}` | `"globe"` only - the boundaries, keys below |

`subsampling` keys: `method` (`"none"`, `"farthest_point"`, `"uniform"`), `num_points` (0 = no budget), `bounding_box_min`/`bounding_box_max` (`[]`, in the model part's own coordinates), `seed` (-1). Unselected nodes take their nearest selected node's prediction.

`branch_input` keys: `process_info_variables` (`[]`), `properties_variables` (`[]`), `constants` (`[]`), concatenated in that order.

`globe` keys: `boundary_sub_model_parts` (`[]`), `boundary_fields` (`[]`, entries `variable_name`/`data_location`/`rank`), `reference_lengths` (`{}`), `output_names` (`[]`), `source_container` (`"Elements"`).

### `rom_surrogate_process`

Case parameters in, modal coefficients out, full field reconstructed through a `RomApplication` POD basis (`u = Phi q`). `inference_process` settings plus:

| Key | Default | Meaning |
|---|---|---|
| `rom_basis_folder`, `rom_basis_name` | `"rom_data"`, `"RomParameters"` | the `CalculateRomBasisOutputProcess` output |
| `input_reduction` | `"mean"` | how per-node inputs become the parameter vector |
| `output_fields` | refused | derived from the basis's nodal unknowns |

### `graph_inference_process`

MeshGraphNet-family models on the mesh's own element-edge graph, extracted once. Optional multiscale hierarchy and proximity "world" edges.

```json
{ "python_module" : "graph_inference_process",
  "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.inference",
  "Parameters"    : { "model_part_name" : "FluidModelPart", "model_interface" : "meshgraphnet",
                      "model_settings"  : { "checkpoint_file" : "mgn.mdlus", "checkpoint_type" : "physicsnemo" },
                      "input_fields"    : [ { "variable_name" : "VELOCITY", "data_location" : "node_historical" } ],
                      "output_fields"   : [ { "variable_name" : "PRESSURE", "data_location" : "node_historical" } ] } }
```

| Key | Default | Meaning |
|---|---|---|
| `model_interface` | `"meshgraphnet"` | `"meshgraphnet"`, `"meshgraphkan"`, `"bistride"`, `"hybrid"` |
| `source_container` | `"Elements"` | |
| `update_edge_features` | `false` | recompute relative positions each step (moving meshes) |
| `multiscale_levels` | `1` | levels of the bistride hierarchy |
| `world_edges` | `{"type": "radius", "radius": 0.1, "max_neighbors": 16, "backend": "auto"}` | the proximity edges of the hybrid interface |

### `time_series_inference_process`

An autoregressive surrogate over nodal states: the last `history_size` states (oldest first) in, the next state out, written back so the next step sees it.

| Key | Default | Meaning |
|---|---|---|
| `history_size` | `2` | states in the rolling window; the first steps are warm-up |

### `grid_inference_process`

Same-resolution grid-to-grid models (FNO, AFNO, UNet, ModAFNO, GraphCast through the squeeze idiom): fields sampled on a grid, forward, scattered back to the nodes. A thin subclass of `superresolution_process` with one model part and one grid.

| Key | Default | Meaning |
|---|---|---|
| `model_part_name`, `grid_shape` | mandatory, `[16, 16, 16]` | replace the coarse/fine pair, which is refused |
| `bounding_box`, `squeeze_axis`, `model_interface` | `[]`, `-1`, `"grid"` | as in `superresolution_process`; `"modafno"` feeds the solver `TIME` as the timestep input |

### `superresolution_process`

Coarse-mesh solve in, fine-mesh field out: the coarse part sampled onto a grid, an SRResNet-style model, the fine grid scattered onto the fine part.

| Key | Default | Meaning |
|---|---|---|
| `coarse_model_part_name`, `fine_model_part_name` | mandatory | |
| `coarse_grid_shape` | `[8, 8, 8]` | the model decides the fine shape |
| `bounding_box` | `[]` | shared by both grids |
| `squeeze_axis` | `-1` | axis of size 2 to average away for planar cases (`-1` none) |
| `model_interface` | `"grid"` | |

### `sequence_inference_process`

One-to-many grid sequence models (`One2ManyRNN`, or FNO with a fourth dimension): seeded once from the solver's grid, then rolled forward.

| Key | Default | Meaning |
|---|---|---|
| `grid_shape`, `bounding_box`, `squeeze_axis` | `[8, 8, 8]`, `[]`, `-1` | |
| `window_size` | `2` | states per window |
| `window_as_time_axis` | `false` | present the window as a fourth spatial axis (FNO dimension 4) |

`model_interface` (`"grid"`, or `"dpot"` for `physicsnemo.models.dpot.DPOTNet`, whose `(B, H, W, T, C)` layout the process permutes to and from).

### `diffusion_inference_process`

Conditional diffusion on grids: the condition sampled from the model part, an ensemble generated, its mean written to the output fields and its standard deviation to the uncertainty fields, optionally onto a second (finer) model part.

| Key | Default | Meaning |
|---|---|---|
| `output_model_part_name` | `""` | write onto another part (downscaling); empty means the input part |
| `uncertainty_fields` | `[]` | one entry per output field, receives the ensemble std |
| `grid_shape`, `bounding_box`, `squeeze_axis` | `[8, 8, 2]`, `[]`, `-1` | |
| `denoiser_interface` | `"edm"` | `"edm"`, `"dit"`, `"unet3d"`, `"topodiff"` (a model already speaking the protocol contract uses `"protocol"` in `sampler_settings`) |
| `sampler_settings` | `{}` | ensemble size, steps, seeds, solver, and the `"guidance"` and `"patching"` blocks below |

`sampler_settings` keys: `num_samples` (8), `num_steps` (18, at least 2), `solver` (`"heun"`, `"euler"`, `"edm_stochastic_euler"`, `"edm_stochastic_heun"`), `solver_options` (`{}`, reaching the solver's constructor), `output_channels` (0 = the model's own), `seed` (-1), `api` (`"auto"`), `denoiser_interface` (`"edm"`), `preconditioner` (`"none"`), `sigma_min` (0.002), `sigma_max` (80.0), `rho` (7.0), `sigma_data` (0.5), `guidance` (`{}`), `patching` (`{}`). The separate `"regression_settings"` block, alongside `sampler_settings`, adds a CorrDiff regression mean.

`guidance` keys: `type` (`"none"`, `"data_consistency"`, `"model_consistency"`), `std_y` (0.1), `gamma` (0.0), `norm` (2), `observation_fields` (`[]`), `mask_field`, `operator` (`"kratos_residual"`), `residual_fields` (`[]`), `residual_model_part_name` (`""`), `use_stored_fixed_values` (`true`).

`patching` keys: `patch_shape` (`[]` = off), `overlap_pix` (0), `boundary_pix` (0), `chunk_size` (0). 2-D latents only.

### `particle_inference_process`

Learning-to-Simulate on a particle cloud: velocity history and node-type one-hots in, per-particle acceleration out, semi-implicit Euler moving the nodes, the proximity graph rebuilt every step.

| Key | Default | Meaning |
|---|---|---|
| `model_interface` | `"meshgraphnet"` | |
| `connectivity` | `{}` | `"type"` (`"radius"` or `"knn"`), `"radius"` (default `0.015` - match it to the particle spacing or the graph is empty), `"max_neighbors"`, `"backend"` (`"auto"`, `"numpy"`, `"warp"`, `"cupy"`), and `"box_size"` (`[]`, `[L]` or `[Lx, Ly, Lz]`) for a periodic minimum-image search |
| `history_size` | `2` | velocity states kept |
| `node_type_variable`, `num_node_types` | `""`, `0` | an integer variable one-hot encoded into the features |

### `domino_inference_process`

A pretrained or fine-tuned DoMINO on the current state: the surface and volume exported as one CAE datapipe case per execution, the per-node and per-triangle predictions written back onto the parent entities, de-normalized through the checkpoint's scaling factors.

| Key | Default | Meaning |
|---|---|---|
| `surface_model_part_name` | mandatory | the skin; `surface_source_container` `"Conditions"` |
| `volume_model_part_name` | `""` | needed for `model_type` `"volume"` or `"combined"` |
| `model_type` | `"surface"` | |
| `bounding_box`, `bounding_box_surface` | `[]` | must match the checkpoint's; the surface box is *not* the volume box for the public checkpoints |
| `global_params_order` | `[]` | see the CAE exporter |
| `output_fields_surface`, `output_fields_volume` | `[]` | field entries (surface defaults to `condition` location) |
| `scaling_factors_file`, `normalization`, `redimensionalize` | `""`, `"none"`, `false` | the pretrained checkpoints emit dimensionless min-max fields; set the file, `"min_max"` and `true` to get pascals back |
| `scratch_directory` | `""` | where the per-step case is written |
| tessellation keys | | |

### `nim_inference_process`

The same no-weights-no-GPU property against NVIDIA's packaged DoMINO-Automotive-Aero microservice: the skin goes out as an STL upload, point clouds come back and are scattered by nearest neighbour.

| Key | Default | Meaning |
|---|---|---|
| `base_url`, `endpoint` | `"http://localhost:8000"`, `"infer"` | |
| `api_key`, `api_key_env` | `""` | bearer token for hosted gateways; local containers need none |
| `timeout` | `120.0` | seconds |
| `stream_velocity`, `stencil_size`, `point_cloud_size` | `30.0`, `1`, `500000` | the documented form fields |
| `surface_output_fields`, `volume_output_fields` | `[]` | entries with a `nim_key` naming the returned array plus the usual variable and location |

### `pinn_solve_process`

No checkpoint: the process *is* the solve. A network is trained on the model part's nodes against a SymPy PDE and the Dirichlet fixities, once, in `ExecuteBeforeSolutionLoop`; inverse mode recovers a coefficient from observation fields.

| Key | Default | Meaning |
|---|---|---|
| `mode` | `"forward"` | or `"inverse"` |
| `physics` | `{}` | the PDE (`"builtin:diffusion"` and friends, or a module path) and its coefficients |
| `fields` | `[{"name": "u", "width": 1}]` | the network's outputs |
| `solution_fields`, `output_fields` | `TEMPERATURE` at `node_historical` | where fixities are read and where the solution is written |
| `observation_fields` | `[]` | inverse mode data |
| `network` | `{"layer_size": 64, "num_layers": 4, "activation_fn": "silu"}` | |
| `training` | epochs `500`, learning rate `1e-3`, `collocation_points` `0`, weights physics `1.0` / boundary `10.0` / data `1.0`, `seed` `0` | |
| `device`, `normalize_coordinates` | `"auto"`, `true` | normalization happens inside the model so autodiff stays physical |

## Judging, adapting, differentiating

Package `KratosMultiphysics.PhysicsNeMoApplication.processes`.

### `validation_metrics_process`

Compares predicted fields against reference fields every `output_interval` steps and writes a JSON report in `ExecuteFinalize`.

```json
{ "python_module" : "validation_metrics_process",
  "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes",
  "Parameters"    : { "model_part_name" : "ThermalModelPart", "output_file" : "validation_metrics.json",
                      "list_of_comparisons" : [ { "predicted_variable" : "TEMPERATURE", "predicted_location" : "node_historical",
                                                  "reference_variable" : "TEMPERATURE_REFERENCE", "reference_location" : "node_historical",
                                                  "metrics" : ["mse", "rmse", "max_abs_error", "relative_l2"] } ] } }
```

| Key | Default | Meaning |
|---|---|---|
| `list_of_comparisons` | one placeholder | entries with `predicted_variable`/`predicted_location`, `reference_variable`/`reference_location`, optional `weight_variable`/`weight_location`, and `metrics` from `mse`, `rmse`, `max_abs_error`, `wasserstein`, `relative_l2`, `relative_mse`, `weighted_mse`, `weighted_rmse`, `histogram_l1`, `entropy_difference` |
| `cfd_metrics` | `[]` | entries `{"name", "domain"}` from physicsnemo-cfd's registry plus a free-form `fields` block |
| `uncertainty_comparisons` | `[]` | mean and std variables against a reference: `coverage`, `nll`, `sharpness`, `calibration_error`, with `confidence_z` |
| `ensemble_comparisons` | `[]` | explicitly named `member_variables` (at least two): `crps`, `kcrps` |
| `output_interval`, `output_file` | `1`, `"validation_metrics.json"` | |

`list_of_comparisons` entries also take `bins` (32) for the distribution metrics. A `spectral_comparisons` block compares power spectra, with entries: `predicted_variable`/`predicted_location`, `reference_variable`/`reference_location`, `grid_shape` (`[32, 32, 2]`), `bounding_box` (`[]`), `squeeze_axis` (2), `high_wavenumber_fraction` (0.5) and `metrics` (`power_spectrum_relative_l2`, `high_wavenumber_energy_ratio`). Both default to `[]`, so an existing configuration gains nothing it did not ask for.

### `adaptive_remesh_process`

Closes the residual-scoring to mesh-adaptation loop: the real solver residual of the surrogate's field becomes an equidistributed size field, handed to `MeshingApplication`'s MMG.

| Key | Default | Meaning |
|---|---|---|
| `remesh_interval` | `1` | steps between remeshes |
| `size_settings` | `{}` | `"target_error"` (`1e-3`), `"exponent"` (`0.5`), `"minimal_size"` (`1e-4`), `"maximal_size"` (`1.0`), plus the MMG block (`discretization_type`, `framework`, interpolation flags) - see [Mesh Bridge](../Mesh_Bridge/Mesh_Bridge.html) |
| `echo_level` | `0` | |

### `adjoint_sensitivity_process`

Puts a sensitivity field (dJ/dX by default) onto the model part as an ordinary nodal variable, so every exporter carries it and a dataset can hold gradient targets.

| Key | Default | Meaning |
|---|---|---|
| `sensitivity_source` | `"shipped"` | the application's adjoint, or `"response_function"` for Kratos's own adjoint stack through `bridges.adjoint_bridge` |
| `objective` | `{}` | the objective (`weighted_sum` over fields, ...) for the shipped source |
| `dof_fields` | `[]` | the state fields the residual is assembled from |
| `fd_step` | `1e-6` | step of the per-entity local finite difference |
| `design_sub_model_part_name` | `""` | restrict the design nodes |
| `output_variable`, `output_location` | `"SHAPE_SENSITIVITY"`, `"node_non_historical"` | where the field goes |
| `response_settings` | `{}` | the Kratos response definition for the `"response_function"` source |
| `execution_point` | `"finalize_solution_step"` | also `"initialize_solution_step"` or `"before_output_step"` |
| `output_interval`, `echo_level` | `1`, `0` | |

Next: [Where things live](Module_Map.html) for the non-process modules, or the topic pages linked from each section.
