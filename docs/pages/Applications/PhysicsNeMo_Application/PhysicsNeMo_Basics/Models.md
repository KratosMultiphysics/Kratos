---
title: Models
keywords: physicsnemo models architectures fno meshgraphnet transolver domino
tags: [Models.md]
sidebar: physicsnemo_application
summary: The 25 architecture families under physicsnemo.models, and which ones this application deploys.
---

# Models

`physicsnemo.models` holds 25 architecture families. You do not need to know them all; you need to know **which shape of data you have**, because that is what picks the family.

## Pick by data shape

| Your data is | Use | Deployed here by |
|---|---|---|
| A regular grid, same resolution in and out | FNO, AFNO, UNet (DPOT fits mechanically, not yet pinned) | `grid_inference_process` |
| A coarse grid in, a fine grid out | SRResNet (`srrn`) | `superresolution_process` |
| An unstructured mesh with connectivity | MeshGraphNet family | `graph_inference_process` |
| An unordered point cloud | Transolver, FIGConvNet, FLARE | `point_cloud_inference_process` |
| A CAD surface plus a volume (external aero) | DoMINO | `domino_inference_process` |
| A time series of states | RNN, `mesh_reduced` | `sequence_inference_process`, `time_series_inference_process` |
| Particles with trajectories | MeshGraphNet, VFGN | `particle_inference_process` |
| A distribution, not a single answer | diffusion U-Nets, DiT | `diffusion_inference_process` |
| The globe | GraphCast, Pangu, FengWu, DLWP | `grid_inference_process` (GraphCast) |

The same choice as a chart - start from the shape of one sample and follow the arrows to the process that deploys it:

<div class="mermaid">
flowchart TD
    start([What does one sample look like?]) --> grid{a regular grid?}
    grid -->|same resolution in and out| fno[FNO, AFNO, UNet]
    grid -->|coarse in, fine out| srrn[SRResNet]
    grid -->|a time series of grids| rnn[One2ManyRNN, FNO dimension 4]
    grid -->|a distribution, not one answer| diff[diffusion U-Nets, DiT, DiffusionUNet3D]
    start --> mesh{an unstructured mesh?}
    mesh -->|use the connectivity| mgn[MeshGraphNet, BiStride, Hybrid, KAN]
    mesh -->|ignore it, points only| pc[Transolver, GeoTransolver, FLARE, FIGConvUNet]
    mesh -->|a CAD surface plus a volume| domino[DoMINO]
    mesh -->|a time series of nodal states| ts[any per-node model, rolled forward]
    start --> other{something else?}
    other -->|particles with trajectories| part[MeshGraphNet on a proximity graph, VFGN]
    other -->|a few parameters to a whole field| mlp[FullyConnected, or a POD basis plus a small net]
    other -->|no data, only the PDE| pinn[a PINN network]
    fno --> p1[grid_inference_process]
    srrn --> p2[superresolution_process]
    rnn --> p3[sequence_inference_process]
    diff --> p4[diffusion_inference_process]
    mgn --> p5[graph_inference_process]
    pc --> p6[point_cloud_inference_process]
    domino --> p7[domino_inference_process]
    ts --> p8[time_series_inference_process]
    part --> p9[particle_inference_process]
    mlp --> p10[inference_process, rom_surrogate_process]
    pinn --> p11[pinn_solve_process]
</div>

## The families

### Neural operators on grids

| Family | Classes | Notes |
|---|---|---|
| `fno` | `FNO` | The reference neural operator. `dimension=1..4`; the 4-D case is a spatio-temporal block operator |
| `afno` | `AFNO`, `ModAFNO` | Fourier operator with an adaptive token mixer. `ModAFNO` takes a *timestep* input, so the solver's TIME conditions it |
| `dpot` | `DPOTNet` | Denoising pre-trained operator transformer, 2-D and 3-D |
| `unet` | `UNet` | A plain 3-D convolutional U-Net |
| `pix2pix` | `Pix2Pix` | Image-to-image translation |
| `srrn` | `SRResNet` | Super-resolution residual network — coarse in, fine out |

### Graphs and meshes

| Family | Classes | Notes |
|---|---|---|
| `meshgraphnet` | `MeshGraphNet`, `BiStrideMeshGraphNet`, `HybridMeshGraphNet`, `MeshGraphKAN` | Encode–process–decode on an edge graph. The variants add a multiscale hierarchy, proximity "world" edges, and KAN layers |
| `mesh_reduced` | `Mesh_Reduced`, `Sequence_Model` | Reduce a mesh, then learn the dynamics in the reduced space with temporal attention |
| `graphcast` | `GraphCastNet` | The weather architecture, on an icosahedral grid |
| `vfgn` | `VFGNLearnedSimulator` and its encode/process/decode parts | Virtual Foundry GraphNet, for sintering and additive manufacturing |

### Point clouds and transformers

| Family | Classes | Notes |
|---|---|---|
| `transolver` | `Transolver` | Attention over learned "physics slices" of a point cloud |
| `geotransolver` | `GeoTransolver` | Geometry-aware successor |
| `figconvnet` | `FIGConvUNet` | Factorized implicit grids; per-point fields plus a scalar (drag-style) head |
| `flare` | `FLARE` | Experimental attention successor (also under `experimental.models`) |
| `domino` | `DoMINO` | External aerodynamics: CAD surface plus volume, with pretrained checkpoints |

### Generative

| Family | Classes | Notes |
|---|---|---|
| `diffusion_unets` | `SongUNet`, `SongUNetPosEmbd`, `DhariwalUNet`, `CorrDiffRegressionUNet`, `UNet`, `StormCastUNet` | The denoiser backbones. **2-D image oriented** |
| `dit` | `DiT` | Diffusion transformer — the denoiser this application wraps for non-image data |
| `topodiff` | `TopoDiff` | Diffusion for topology optimization — deployed through `denoiser_interface: "topodiff"`, its constraint channels passed natively as the model's second argument |

A **volumetric** 3-D denoiser exists too, under `physicsnemo.experimental.models.diffusion_unets.DiffusionUNet3D` — this application deploys it through `denoiser_interface: "unet3d"` (see the [Diffusion](../Diffusion/Diffusion.html) page and [Companion packages](Companion_Packages.html)).

### Sequences and weather

`rnn` (`One2ManyRNN`, `Seq2SeqRNN`), `dlwp` (`DLWP`), `dlwp_healpix` (`HEALPixUNet`, `HEALPixRecUNet`), `pangu` (`Pangu`), `fengwu` (`Fengwu`), `swinvrnn` (`SwinRNN`).

### The plain one

`mlp` (`FullyConnected`) — a multilayer perceptron. It is the right first model far more often than it looks: if your input is a handful of case parameters and your output is a field, you want this, not a neural operator.

### In physicsnemo 2.2 but not deployed here

Five families left this table and are now deployed: `dpot` (`DPOTNet`) through `SequenceInferenceProcess`'s `"model_interface": "dpot"`, `topodiff` through `denoiser_interface: "topodiff"`, `experimental.xdeeponet` and `experimental.globe` through the point-cloud process's `"deeponet"` and `"globe"` interfaces, and `experimental.aerojepa` through `training.aerojepa_pretraining`. What remains below has no Kratos counterpart deployed; only `VideoHealDA` is a roadmap row, the rest are listed so the gap stays visible rather than because a bridge is planned.

| Family | Class | What it would bring to Kratos |
|---|---|---|
| `pix2pix` | `Pix2Pix`, `Pix2PixUnet` | a plain convolutional image-to-image translator; fits the grid process mechanically |
| `experimental.strata`, `experimental.healda` | `Strata`, `VideoHealDA` | weather emulation on the sphere and HEALPix data assimilation; the assimilation idea matters for digital twins, the API is calendar-shaped |
| `pangu`, `fengwu`, `swinvrnn`, `dlwp_healpix` | as named | global weather architectures with no Kratos counterpart |

## Building blocks

`physicsnemo.nn` holds the layers the models are made of — around 150 of them. Two are worth knowing by name because this application uses them directly:

- `ConcreteDropout` — dropout with a *learned* rate, which makes MC-dropout uncertainty estimates meaningfully calibrated instead of arbitrary;
- `physicsnemo.nn.functional.derivatives` — differential operators used by the physics-informed path.

## What this application uses it for

| Architecture | Kratos-side entry | Status here |
|---|---|---|
| `FNO`, `AFNO`, `ModAFNO` | `processes.inference.grid_inference_process` | test-pinned |
| `SRResNet` | `processes.inference.superresolution_process` | test-pinned |
| `MeshGraphNet` + all three variants | `processes.inference.graph_inference_process` | test-pinned, one on a real fluid solve |
| `Transolver`, `GeoTransolver`, `FLARE`, `FIGConvUNet` | `processes.inference.point_cloud_inference_process` | test-pinned |
| `DoMINO` | `processes.inference.domino_inference_process` | test-pinned, pretrained checkpoints de-normalized |
| `One2ManyRNN` | `processes.inference.sequence_inference_process` | test-pinned |
| `Sequence_Model` | `training.rom_temporal` | test-pinned |
| `VFGNLearnedSimulator` | `bridges.vfgn_bridge` | driven by a real sintering solve |
| `DiT`, `CorrDiffRegressionUNet` | `processes.inference.diffusion_inference_process` | test-pinned |
| `GraphCastNet` | `processes.inference.grid_inference_process` | documented, needs `torch_sparse` or `dgl` |
| `FullyConnected` | anything | the default in most examples |

User-supplied 2-D UNets go through the grid process mechanically, by the same thin-axis squeeze idiom (physicsnemo's own `UNet` is 3-D and takes the volumetric path); `DLWP` is test-pinned on the volumetric path with its six cubed-sphere faces as the depth axis.

Next: [Data and datapipes](Data_And_Datapipes.html) — feeding them.
