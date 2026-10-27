# PhysicsNeMo Application

|      **Application**      |                                                                              **Description**                                                                              |                                       **Status**                                        |                       **Authors**                        |
|:--------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------:|:-----------------------------------------------------------:|
| `PhysicsNeMoApplication` | The *PhysicsNeMo Application* bridges *Kratos Multiphysics* with [NVIDIA PhysicsNeMo](https://docs.nvidia.com/physicsnemo/latest/index.html), a PyTorch-based framework for physics-ML: surrogate models, neural operators, and active learning with an external ground-truth solver. | <img src="https://img.shields.io/badge/Status-%F0%9F%9A%80%20Actively%20developed-Green"  width="200px"> | [*Vicente Mataix Ferrándiz*](mailto:vicente.mataix-ferrandiz@siemens.com) |

The application includes tests to check the proper functioning of the application.

<p align="center">
  <img src="../../docs/pages/Applications/PhysicsNeMo_Application/General/images/architecture.svg" alt="Kratos on the left, NVIDIA PhysicsNeMo on the right, the application's packages in between, and the artifacts that travel between the steps along the bottom" width="900"/>
</p>

## 🧭 New here? Start with these three:

The feature list below is long because the application is. These are the three
pages that make it navigable:

| | |
|---|---|
| **[PhysicsNeMo Basics](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Overview.html)** | **NVIDIA PhysicsNeMo itself** — what it is, what every one of its modules does, and which part of this application uses it. Nine pages: [core and checkpoints](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Core_And_Checkpoints.html), [models](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Models.html), [data and datapipes](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Data_And_Datapipes.html), [mesh and geometry](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Mesh_And_Geometry.html), [symbolic and physics](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Symbolic_And_Physics.html), [diffusion and deployment](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Diffusion_And_Deployment.html), [distributed and scale](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Distributed_And_Scale.html), [companion packages](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/PhysicsNeMo_Basics/Companion_Packages.html). Read this if the rest of the documentation seems to assume things you were never told. |
| **[Where things live](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Module_Map.html)** | Every module by folder with a one-line purpose, and an *I want to X → use Y* index. |
| **[From scratch](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/From_Scratch.html)** | One surrogate end to end — export from a real solve, train, save with a model card, deploy in the loop, validate — naming the exact module at each step. |

### Documentation map

| Section | Read it for |
|---|---|
| [Architecture](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Architecture.html) | the three columns of the figure above, where processes hook into the solution loop, and which artifact each step writes |
| [Kratos concepts for ML readers](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Kratos_Concepts_For_ML.html) | ModelPart, Variables and data locations, Processes and execution points, TensorAdaptors - and what each becomes on the PhysicsNeMo side |
| [Process reference](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Process_Reference.html) | all twenty-six processes with the JSON that attaches them and every setting's default |
| [Installation and environment](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Installation_And_Environment.html), [Testing and contributing](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Testing_And_Contributing.html) | the dependency matrix with verified versions, the suites and the guards |
| [Performance](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Performance.html), [Troubleshooting and traps](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Troubleshooting_And_Traps.html), [Glossary](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Glossary.html) | what the bridge costs, every failure the history produced as symptom-cause-fix, and the vocabulary |
| [Examples](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/Examples/Examples.html) | the gallery of twenty-six use cases and the twenty notebooks |
| one page per topic | tensor and mesh bridges, training, inference, active learning, superresolution, graphs, particles, sequences, point clouds, CAE datapipes, ROM, adjoint, diffusion, uncertainty, physics-informed, co-simulation, distributed - each with a diagram |

### Layout

`python_scripts/` is a tree of packages, and the folder a module sits in says what kind of thing it is. Each package's `__init__.py` says what belongs in it; a process is attached the usual way, with `kratos_module` naming the package.

<details>
<summary>Package tree and the process-attachment snippet</summary>

```
python_scripts/
├── processes/          everything you attach to a solve - and only those have a Factory
│   ├── inference/          run a trained model in the solution loop      (16)
│   ├── export/             write solver data out as training data         (7)
│   └── ...                 adaptive remeshing, validation, adjoint dJ/dX  (3)
├── bridges/            Kratos data <-> physicsnemo data     (11 + mesh_bridge/)
├── training/           loops, datasets, schemes                           (9)
├── physics/            residuals, PINN machinery, sensitivities           (4)
├── deployment/         checkpoints, cards, ONNX/Triton/NIM/USD, uncertainty (9)
├── distributed/        MPI, halo-partitioned graphs, ShardTensor          (3)
├── active_learning/    Kratos as the labeling oracle
└── utilities/          small shared helpers
```

```json
{
    "python_module" : "inference_process",
    "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.inference"
}
```

</details>

## 😎 Features:

- **Tensor bridge**: zero-copy `ModelPart` data ↔ `torch.Tensor` built on the core `Kratos.TensorAdaptors`, plus a dataset-export process producing `.npz` training samples from any combination of nodal, elemental and Gauss-point fields, and a `torch.utils.data.Dataset` factory over the exported samples for direct training — plus a **streaming path** (`streaming_dataset`) that trains directly out of a running solve with no file round trip, yielding samples byte-identical to the dumped ones (asserted by running one case both ways), and **shrink-and-perturb warm restarts** in `TrainModel` for when Kratos data drifts to a new geometry family. Dataset curation for the `.pmsh` mesh series: **coherent** random rotate/scale/translate augmentations (`MakeMeshAugmentations`/`CreateAugmentedMeshDataset` — vector and rank-2 tensor fields transform *with* the coordinates, which upstream's defaults silently skip, plus the dtype cast its transforms require and reproducible per-epoch seeding) and `MultiDataset` mixing of several series (`CreateMultiMeshDataset`).

- **Optional CuPy acceleration** (`utilities/array_backend_utils`): a `"backend"` switch on the array-heavy bridge paths — the particle proximity graph, graph edge features, trilinear grid interpolation and the ROM basis projection — with **numpy the default and the fallback everywhere**, a measured per-site size threshold below which a CuPy request falls back on its own, and `AsTorchTensor`, which hands a device array to torch through DLPack rather than bouncing it through the host. `cupy` is optional exactly as `torch` is: nothing imports it at module scope, and `benchmarks/benchmark_bridges.py --backend both` reproduces every number claimed for it.

- **Mesh bridge**: tessellation of arbitrary Kratos meshes (hexahedra, prisms, pyramids, quadrilaterals, higher-order elements) into the simplicial representation required by `physicsnemo.mesh` — watertight on unstructured meshes via the smallest-node-id diagonal rule (Dompierre et al.), with optional subdivision of quadratic geometries through their real mid-side nodes and an opt-in **curved (isoparametric) mode** that samples the exact quadratic geometry on a configurable refinement lattice with synthetic points (interpolated on gather, dropped on scatter-back), watertight across curved neighbours — plus a field-provenance map that round-trips predictions back onto the original Kratos entities; `DomainMesh` export with named boundaries from sub-model-parts, and native save/load of PhysicsNeMo's memory-mapped mesh format. On top of it:
    - **Discrete calculus** (`calculus_bridge`) — LSQ/DEC gradient, divergence, curl, Laplacian and integrals on the tessellated mesh, autograd-differentiable, with the backend-validity guards and boundary masks the upstream operators lack.
    - **Mesh generation from implicit geometry** (`mesh_bridge/generate`) — SDF primitives and combinators into a meshed volume, marching-cubes level sets, quality-guaranteed 2D loop filling, **tetrahedral filling of watertight 3D surfaces** (winding-number-carved, so non-convex solids fill correctly, and self-validating against the input's own volume and boundary area; an opt-in **`"method": "tetgen"`** backend adds **exact boundary recovery** — input facets kept verbatim, Schönhardt-class solids filled via interior Steiner points — through the optional AGPL-licensed `tetgen` package, never chosen by `"auto"`), differentiable refitting, and `PopulateModelPartFromMesh`, the direction the bridge never had: generated geometry materialized as real Kratos entities and handed straight to MMG or a solver.
    - **Exact NURBS sampling** (`mesh_bridge/nurbs_sampling`) — the isogeometric analogue of the curved mode: `IgaApplication`-style NURBS surfaces and volumes (core `NurbsSurfaceGeometry3D`/`NurbsVolumeGeometry`) sampled on a parametric lattice with every point evaluated by the C++ geometry itself (rational arcs land on the circle to machine precision), tessellated watertight (Kuhn 6-tet split), with the **isogeometric gather** — control-point fields interpolated through the geometry's own NURBS basis via `CreateQuadraturePointGeometries`, since a control point is generally not on the geometry.
    - **Grid-stencil counterparts** in `grid_bridge.ComputeGridDerivatives`/`ComputeGridVectorOperator` — gradients, divergence, curl and Laplacian, with the float64-precision backend choice upstream gets wrong by default.
    - **Signed distance fields as features** (`mesh_bridge/spatial`) — written into an ordinary nodal variable so every existing gather picks them up, on a boundary surface this bridge re-orients because upstream's extraction winds it inconsistently.
    - **Surrogate-error-driven adaptive remeshing** (`adaptive_remeshing`/`AdaptiveRemeshProcess`) — residual scoring → equidistributed size field → `MeshingApplication` MMG scalar-metric remesh, plus error-weighted `partition_cells` surface clustering and a Warp-backed `remesh` wrapper.
    - **A differentiable shape-deformation layer** (`mesh_bridge/deformation`) — FFD/RBF/morph/displace parameterizations, mesh-quality energies including the inversion term that keeps an optimizer from tearing the mesh, and coordinate write-back through the position tensor adaptor; the matching `sensitivity_utils.ComputeShapeSensitivities`, which gives `dJ/d(control)` in one backward pass, lives with the other adjoints.
    - **Mesh operations** (`mesh_bridge/operations`) — curvature written onto the nodes as a feature beside the SDF, `repair_mesh` with its report, loop/butterfly/linear subdivision, Laplacian smoothing, the P0 `integrate_moment`, and `ExtrudeModelPart` sweeping a planar Kratos case into a 3-D one as real entities; `DeformMesh` adds 2.2's two *mesh*-level deformers, `shrinkwrap` and the Sobolev displacement filter.
    - **Torch-native sampling over the mesh's own BVH** (`mesh_bridge/sampling`), with a `backend` switch on `grid_bridge.SampleFieldsOnGrid`. The roadmap asked for a measurement before any switch, and the measurement says do not switch: the Kratos locator is 2-3x faster (795 ms against 1767 ms on a 64³ lattice over 196,608 tetrahedra), so the default stays `"kratos"` and the BVH path is kept for the gradient it carries, which the locator cannot.
    - **AI-ready Zarr without physicsnemo-curator** — `SaveMeshZarr`/`LoadMeshZarr` round-trip a `Mesh` or a `DomainMesh` with its boundary names through `mesh.io.to_zarr`, and `CuratorExportProcess` defaults to that backend; curator is git-only, builds a Rust toolchain at install and rejects a `DomainMesh` in its sinks.
    - **Kratos VTK output as training data** (`bridges/vtk_bridge`) — read-only and deliberately so, since a VTK file carries no provenance and nothing can be scattered back through it.

- **In-loop inference**: an `InferenceProcess` running a trained model (TorchScript or physicsnemo checkpoint) at a configurable point of the solution loop, writing predictions into existing Kratos variables, and a `HybridInitializationProcess` warm-starting solves from a model prediction. Deployment extras: **ONNX export and inference** (`ExportOnnxModel` via `physicsnemo.deploy.onnx` plus an `OnnxInferenceProcess` running a cached ONNX Runtime session — a portable deployment artifact needing neither torch-the-model-runtime nor physicsnemo), opt-in **`torch.compile(fullgraph=True)`** wrapping of loaded physicsnemo checkpoints, and opt-in **NVTX ranges** around every deployment process's gather/forward/scatter hot paths for Nsight Systems profiling. **Serving**: `triton_export` writes a Triton Inference Server model repository (dynamic-entity-axis ONNX or TorchScript plus a generated `config.pbtxt`, tensor names taken from the model card) and `TritonInferenceProcess` calls the running server from inside the solution loop, so the solver host needs neither weights nor a GPU. **NIM client** (`nim_client`/`NimInferenceProcess`): the same no-weights-no-GPU property against NVIDIA's packaged physics microservices — the model part's triangle skin goes out as the documented multipart STL upload (`design_stl` + `stream_velocity`/`stencil_size`/`point_cloud_size`), the returned `.npz` point clouds come back onto the Kratos nodes by nearest neighbour (the NIM samples its own points — no exact provenance exists on this path), scalars like `drag_force`/`lift_force` are collected; stdlib `urllib` + numpy only, payload contract pinned against a stub transport since running a NIM needs an NGC key and docker — the one part of this feature no test here has seen live.

- **Superresolution**: a grid bridge sampling unstructured-mesh fields onto regular voxel grids (FE shape-function interpolation, exact for linear fields) and a `SuperResolutionProcess` deploying `physicsnemo.models.srrn.SRResNet`-style models: coarse-mesh solve in, superresolved fine-mesh field out.

- **Graph neural networks on the mesh**: a graph bridge extracting the true element-edge graph of a `ModelPart` (bidirectional, MeshGraphNet edge encoding) and a `GraphInferenceProcess` deploying `physicsnemo.models.meshgraphnet.MeshGraphNet`-style models directly on unstructured meshes (requires the optional `torch_geometric`/`torch_scatter`), plus a `TimeSeriesInferenceProcess` for autoregressive transient surrogates.
    - The **scalable external-aero variants** deploy through the same process via `model_interface`: `BiStrideMeshGraphNet`'s multiscale U-Net with a scipy-based hierarchy builder replacing an upstream one blocked on `sparse_dot_mkl`, `HybridMeshGraphNet`'s mesh+proximity "world" edges, and `MeshGraphKAN` — validated on a real in-memory lid-driven-cavity FluidDynamics solve.
    - **Temporal training schemes** (`temporal_training`) — single-step, time-conditional and one-shot window datasets sharing the deployment window convention, plus `TrainAutoregressive` backpropagating through a self-fed rollout with optional per-step gradient checkpointing, driven by real transient Kratos solves.
    - The **GraphCast** grid surrogate recipe (`physicsnemo.models.graphcast.GraphCastNet` through `GridInferenceProcess`'s squeeze idiom, batch-size-1 training, a numpy-only shallow-water reference integrator generating the trajectories) is documented and test-pinned, gated on the extra `torch_sparse` (or `dgl`) dependency.

- **Grid-to-grid models and real-solver integration**: a `GridInferenceProcess` deploying FNO/UNet-style same-resolution grid models — 3D and, via the thin-axis squeeze idiom, the 2D operator zoo (`FNO(dimension=2)` and AFNO are test-pinned through the squeeze path, `DLWP` through the 3D path with its six cubed-sphere faces as the depth axis; physicsnemo's own `UNet` is 3D-only and takes the `(C, D, H, W)` path, a user-supplied 2D UNet the squeeze path) plus time-modulated ModAFNO (`model_interface: "modafno"` feeds the solver TIME as the timestep input) and spatiotemporal block operators (`FNO(dimension=4)` through `SequenceInferenceProcess`'s window-as-time-axis mode) — and integration tests/examples driven by **real Kratos solves** — ConvectionDiffusionApplication stationary heat conduction and a StructuralMechanicsApplication cantilever (crash/deformation-surrogate pattern) — conditioned on the compiled applications via the standard `CheckIfApplicationsAvailable` pattern.

- **Seven more compiled Kratos applications, each with a case**: a lid-driven-cavity PINN through `PinnSolveProcess`'s builtin incompressible Navier–Stokes (relative L2 about 0.5 against the VMS solve — the cavity's structure, not its detail), transient vortex shedding on a `FluidDynamicsApplication` cylinder, a learned `RANSApplication` closure written into the running loop (4.2e-4 against the mean predictor's 1.0e-1), an aerodynamic surrogate over a Mach sweep on compressible potential flow (3.4e-4 against 4.7e-2), `GeoMechanicsApplication` consolidation (9.8e-5 against persistence at 3.2e-2), a `ContactStructuralMechanicsApplication` patch test, and a real isogeometric solve feeding the isogeometric gather. The solver cases **reuse each owning application's own test data in place** and gate on its presence rather than copying it, which keeps a 216 kB cylinder mesh out of this application and stops the fixtures drifting from the cases they came from.

- **Non-matching transfer and MPI scale**: a `mapping_bridge` moving fields between the solver mesh and ML grids/meshes that match no tessellation, via MappingApplication's mappers (`nearest_element`, `barycentric`, ..., MPI-capable), with structured background-grid model parts and `(C, D, H, W)` gathering.
    - MPI-aware `DatasetExportProcess` **and** `MeshExportProcess` — ghost-free `DataCommunicator` gathers of fields, and full mesh **topology** reconstructed on a rank-0 shadow part via the reusable `GatherModelPartToRank0` primitive, with rank 0 writing the exact serial file layout.
    - Matched **process groups and device meshes** (`CreateMatchedProcessGroup(s)`, `InitializeDeviceMesh`) pairing physicsnemo subgroups with registered Kratos sub-communicators over the same ranks.
    - **Halo-partitioned distributed graph training** (`graph_partition_utils`) — per-rank subgraphs whose owned sets partition the global node set exactly and whose one-hop neighbourhoods match a serial run, fixing a real interface truncation in the plain graph bridge, plus gloo `DistributedDataParallel` gradient sync asserted bit-identical across ranks. `ExchangeHaloFeatures` moves the feature half of that exchange through `physicsnemo.distributed.autograd.indexed_all_to_all_v`, so a halo-subgraph loss backpropagates onto the owning rank — verified at np=2/3 against the eager exchange for values and an analytic expectation for gradients. Gloo cannot perform a ragged all-to-all and a halo is ragged by definition, so every block is padded to one size and sliced back, the padding rows repeating an existing index so their discarded copies carry exactly zero gradient.
    - **Domain parallelism over `ShardTensor`** (`domain_parallel_utils`) — one Kratos field or grid split across the ranks, uneven shards allowed, with physicsnemo's halo-exchanging convolution and mesh-wide reductions on top; asserted at np=2/3 over gloo (an uneven sharded field gathers to the serial layout, a pointwise forward and a 3×3 halo convolution equal the serial ones, and the backward of a sharded mean leaves the serial gradient on every rank — so no DDP-style all-reduce, which the test pins as wrong). Built on a torch `DeviceMesh` directly, since `DistributedManager.initialize_mesh` puts the mesh on CUDA whenever a GPU is visible; the NCCL transport stays untested here.
    - An MPI test suite (`tests/test_PhysicsNeMoApplication_mpi.py`, no Metis required).

- **Sequence, point-cloud and diffusion models**: three deployment families sharing the same gather/scatter contract.
    - **Sequences**: a grid-series exporter (`GridDatasetExportProcess`) with sequence/pair dataset factories feeding `physicsnemo.models.rnn` one-to-many surrogates, deployed by a seed-once/roll-forward `SequenceInferenceProcess` — which also carries **DPOT**, the PDE foundation-model architecture, through `"model_interface": "dpot"`, the adapter being a pair of permutations because DPOT reads `(B, H, W, T, C)` where every other grid model here reads `(B, C, T, H, W)`, a mismatch that transposes the domain silently rather than raising; and a **Virtual Foundry GraphNet** bridge for sintering/AM (`vfgn_bridge` — position-sequence rollout through the public `inference()`, driven by a **real thermo-mechanical sintering solve** (ConvectionDiffusion's `CoupledThermoMechanicalSolver`: cooling-driven thermal contraction — the test asserts the body actually contracts and that the shrinkage tracks the cooling rate, rather than merely running), a training path that routes around the upstream `forward()` being unusable in 2.2, and the normalization-statistics card round trip the particle path had been missing).
    - **Point clouds**: a `PointCloudInferenceProcess` running Transolver-style transformers, **GeoTransolver and FLARE** (the experimental FLARE-attention successors, also reachable through `physicsnemo-cfd`'s checkpoint-driven evaluation wrappers via `cfd_bridge`), FIGConvNet (per-point fields plus its drag-style scalar head) and generic point models directly on the mesh nodes. A `"subsampling"` block gives those transformers a **token budget** — a bounding-box filter first, in the model part's own coordinates, then farthest-point sampling or a seeded uniform draw — with every unselected node taking its nearest selected node's prediction, so the written field stays complete. Two further interfaces are not pointwise at all: **`"deeponet"`** maps case *parameters* to a field at the nodes, which is what `RomSurrogateProcess` does through a POD basis, without the basis; and **`"globe"`** maps *boundary* data to the interior the way an elliptic problem is determined by its boundary conditions, fed by `globe_bridge` and trained by `globe_training`'s own loop, because `TrainModel`'s DataLoader batches tensors and no collate function can stack meshes.
    - **Diffusion**: a bridge (`TrainDiffusionModel` with the conditional EDM loss, `DiffusionInferenceProcess`) doing CorrDiff-style downscaling with ensemble-mean predictions and per-node **uncertainty fields** — including the full **CorrDiff two-stage recipe** (`"regression"`/`"residual"` losses on upstream's `RegressionLoss`/`ResidualLoss`, `TrainCorrDiffPair`, and a `"regression_settings"` block on `DiffusionInferenceProcess` adding the regression mean to the generated ensemble) and a test-pinned **FWI-style subsurface inversion** recipe (sparse borehole observations + mask as conditioning, ensemble std as inversion uncertainty) — with **DiT denoisers** (`physicsnemo.models.dit.DiT`) plugging into the same sampler/loss machinery through `denoiser_interface: "dit"` / `WrapDenoiser` (conditioning concatenated into the input channels, sigma broadcast to timesteps) and **volumetric 3D U-Net denoisers** (`physicsnemo.experimental.models.diffusion_unets.DiffusionUNet3D`) through `denoiser_interface: "unet3d"` — conditioning passed natively as the model's `TensorDict` `"volume"` key rather than concatenated, full `(C, D, H, W)` grids with no thin-axis squeeze, trained through `MSEDSMLoss`, which is rank-agnostic. The bridge now runs on physicsnemo 2.2's **protocol API** throughout — `WrapDiffusionModel` adapting a denoiser to `DiffusionModel`, `MSEDSMLoss` over an `EDMNoiseScheduler`, `samplers.sample` and its four solvers — with `"api": "legacy"` still reaching the deprecated modules, because 2.2 ships no protocol equivalent of CorrDiff's `RegressionLoss`/`ResidualLoss`. Three capabilities follow from that migration: **DPS guidance** (`"guidance"` on `sampler_settings`) steering an already-trained denoiser toward masked measurements or any differentiable forward operator — the shipped one being the exact discrete FEM residual, so the solver's own physics grades the generator without retraining — **multi-diffusion patching** (`"patching"`) sampling grids beyond the training resolution, and **TopoDiff generative design** (`denoiser_interface: "topodiff"`), whose constraint channels go in natively as the model's second argument rather than concatenated, trained on a real SIMP compliance solve whose sampled designs go back to Kratos to be re-solved.

- **Pretrained DoMINO, correctly de-normalized, and fine-tuning on top of it**: the public ungated `nvidia/domino_drivaerml` checkpoints load straight through `model_registry`, and `DominoInferenceProcess` now applies the de-normalization they need — a pretrained DoMINO emits *dimensionless* fields, so writing its raw output onto Kratos entities is wrong by roughly three orders of magnitude (a raw `0.1386` is really −609 Pa).
    - It guards the three ways such a checkpoint is silently misfed: a `grid_resolution` disagreeing with the checkpoint's own is now named rather than surfacing as an opaque reshape error, and `normalize_coordinates`/`bounding_box_surface` defaults that produce no error at all are warned about.
    - `domino_finetune` ships both adaptation recipes — the predictor-corrector decomposition with the frozen predictor's output cached (so its cost is independent of the predictor's size) and LoRA adapters that merge back into an ordinary `.mdlus` needing no configuration change to deploy.
- **Sharded (FSDP2) checkpoints**: `SaveTrainedModel` gathers an FSDP2 model's DTensors to full tensors and writes them from rank 0, so a sharded model produces one ordinary `.mdlus` that `model_registry` loads. Saving the DTensors directly — what it used to do — yielded a checkpoint that reported success and could not be loaded, and across ranks each wrote only its own shard. Asserted at np=2 and np=3 over gloo, with the parameters confirmed genuinely split (`local != global`) before the save; only the multi-GPU NCCL transport is untested here.
- **Card-carried output de-normalization**: a model trained on normalized targets emits normalized predictions, and writing those onto Kratos variables as if they were physical is wrong by the training scaling — silently, since the values stay finite. The model card's `"output_normalization"` key travels that scaling with the checkpoint and every process that loads a checkpoint with a card inverts it: `InferenceProcess` and everything on its write path, `ParticleInferenceProcess` (where it matters most, because the prediction is integrated twice into node positions), the graph and time-series processes, the grid writers (per channel along axis 0), the co-simulation surrogate wrapper and `SurrogateResponseFunction` — the last on the written field *and* inside its autograd objective, so `dJ/dX` carries the training scale. Spreads are scaled but never shifted, an absent key is exactly the identity, and a degenerate scale or wrong channel count raises rather than producing silent NaNs. Each covered process is pinned by a mutation that removes the call; two of them (the co-simulation wrapper and the response function) had been documented as covered while passing no normalization at all.
- **GPU ONNX inference**: `"device" : "cuda"`/`"cuda:N"` on a CUDA build of ONNX Runtime, verified against `onnxruntime-gpu` 1.29.0 (a CUDA-13 build) — the CUDA provider genuinely instantiated, GPU and CPU agreeing to ~2.5e-4 relative under TF32. The bridge guards the two failures ONNX Runtime reports as success: a silent fall back to CPU (a missing CUDA build, or a nonexistent device index, both yield a working CPU session — `require_device` turns that into an error) and a dropped device index (`"cuda:1"` previously ran on device 0).
- **Vertex-morphing comparison**: the shipped deformation layer pinned against `ShapeOptimizationApplication`'s `MapperVertexMorphing`, whose reference fields were generated once in a wheel-only environment and committed as a fixture — no rebuild, and no mixing of a GCC-built wheel with this Clang-built core. The two are close relatives, not the same operator: Kratos normalizes by `sum(w)` (an exact partition of unity), physicsnemo by `1 + sum(a)` (a regularized compact Shepard field), so they agree in the dense-control limit and visibly differ when controls are sparse — and `morph` interpolates at a control point where morphing damps.
- **Exact shape gradients and shape optimization**: `sensitivity_utils.ComputeShapeSensitivityField` gives the discretely exact `dJ/dX` at *every* node from **one** pass over the mesh — perturbing each entity's own nodes and re-evaluating only that entity's local right-hand side, since moving a node perturbs only its adjacent entities. That replaces the per-parameter path's `6N` global assemblies with a cost linear in the mesh (measured ~15x faster at 3200 2-D triangles, ~100x at 24k 3-D tetrahedra, independent of the number of design parameters), and `ComputeControlSensitivities` then pushes the field back through the differentiable deformers to give `dJ/d(control)` for FFD/RBF/morph/displace parameterizations — the FEM-exact counterpart of the surrogate-side `ComputeShapeSensitivities`. Notebook 17 drives a gradient descent with it onto a target.
- **Adjoint cross-validation**: the shipped shape sensitivities (`sensitivity_utils` + `differentiable_residual`) are pinned against `StructuralMechanicsApplication`'s entirely separate adjoint stack — `AdjointFiniteDifferencing*` elements, `AdjointNodalDisplacementResponseFunction`, `ResidualBasedAdjointStaticScheme` and `SensitivityBuilder` producing `SHAPE_SENSITIVITY` — with full finite differences through real solves as a third opinion. All three agree to eight significant digits, and the deformation layer's boundary displacement is smoothed into the interior through `MeshMovingApplication` without folding elements.
- **Co-simulation surrogates**: a `cosim_surrogate_solver_wrapper` letting a trained model participate as a *solver* in Kratos co-simulation — referenced by module path in `CoSimulationApplication`'s `ProjectParameters`, exchanging interface fields through its data-transfer layer (MappingApplication mappers included) and convergence accelerators, with the flat `InferenceProcess` contract or any point-cloud interface, optional time ownership, and lazy model loading through the shared model registry (weak and strong Gauss-Seidel coupling both exercised in the tests). `"distributed" : true` runs the surrogate across MPI ranks — Metis-partitioned or, for element-free clouds Metis cannot graph, partition-free — inferring on owned nodes only and synchronizing ghosts, honouring `mpi_settings` so the wrapper can live on a subset of the ranks.

- **CAE datapipes**: a `CaeDatasetExportProcess` writing per-case `.npz` files in the exact layout `physicsnemo.datapipes.cae` consumes (triangulated STL from the mesh bridge, surface/volume fields, DoMINO- and Transolver-style global parameters), with `CreateDoMINODataPipe`/`CreateTransolverDataPipe`/`CreateCaeDataset` factories — DoMINO and Transolver pipelines run on Kratos data out of the box, MPI included. **DoMINO deploys in-loop** through `DominoInferenceProcess`: the current state is exported as a single preprocessed datapipe case per execution, and the model's per-node volume / per-triangle surface predictions are written back onto the Kratos entities (surface triangles collapse onto their parent conditions/elements via the mesh-bridge provenance).

- **Lagrangian particle surrogates**: a `particle_bridge` building proximity graphs (radius/kNN, warp-accelerated with an exact numpy fallback) over particle node clouds in `graph_bridge`'s exact contract, a `ParticleInferenceProcess` running the Learning-to-Simulate loop (velocity history + node-type one-hots in, per-particle acceleration out, semi-implicit Euler advancing and moving the nodes, graph rebuilt every step) for `physicsnemo.models.meshgraphnet.MeshGraphNet`-style models — the process reads any particle node cloud, so `MPMApplication`/`SPHApplication`/`DEMApplication`/`PfemFluidDynamicsApplication` fit the contract, though none of those applications is compiled in the reference environment and the tests drive synthetic clouds, and a `CreateParticleTrajectoryDataset` factory windowing position trajectories into training pairs with normalization statistics. The **Lennard-Jones molecular-dynamics recipe** runs on the same machinery: `"box_size"` makes the proximity search periodic (minimum image, scipy `cKDTree`), `utilities/lennard_jones_reference` is the numpy velocity-Verlet data source when neither OpenMM nor `DEMApplication` is at hand, and the test-pinned recipe trains the generic `MeshGraphNet` — NVIDIA's own Lennard-Jones example uses no molecular architecture either — to per-atom forces and potential energies and deploys it through `ParticleInferenceProcess` unchanged.

- **RomApplication interoperability**: `rom_bridge` consumes `CalculateRomBasisOutputProcess`'s numpy POD bases (the exact interleaved row-order contract; no compiled RomApplication needed to consume), a `RomSurrogateProcess` deploys **neural-augmented reduced bases** in the solution loop (case parameters → modal coefficients → full-field reconstruction `u = Φq`, with the model card validating against the basis's nodal unknowns automatically), and `rom_temporal` brings **mesh-reduced temporal attention in ROM space** — physicsnemo's decoder-only `Sequence_Model` trained teacher-forced on reduced trajectories with autoregressive rollout.

- **Adjoint integration — gradients both ways**: what `rom_bridge` is to RomApplication, `adjoint_bridge` is to Kratos's adjoint stack, and for the same reason: the contract it consumes (`ResponseFunctionInterface`) lives in the **core**, so it needs whichever application owns the response and never a compiled optimization application. `CreateResponseFunction` dispatches `structural_mechanics` / `convection_diffusion` / `compressible_potential_flow` (plus a module-path escape hatch), and `EvaluateResponse` turns the `{entity_id: value}` dicts Kratos hands out — unordered, and keyed by an *adjoint* model part whose iteration order need not be the primal's — into arrays in this application's own row order, the one shared with `graph_bridge`, `DofFieldMap` and every gather in the tree. The bridge is cross-validated node by node against the shipped `sensitivity_utils` adjoint on **both** compiled stacks: the cantilever to 1e-4 relative (Kratos's `semi_analytic` gradient takes a *forward* difference and sets the tolerance) and a diffusion square to 1e-6, out-of-plane row exactly zero. On top of it:
    - **`AdjointSensitivityProcess`** puts dJ/dX on the model part as an ordinary nodal variable and writes no files, deliberately — `DatasetExportProcess`, the streaming exporter, `MeshExportProcess` and the vtu outputs then carry it unchanged, and naming it in `CreateNpzDataset`'s `output_keys` yields gradient-carrying training targets with no new dataset code.
    - **Sobolev training** (`sobolev_training.MakeSensitivityLossTerm`) grades a surrogate on Kratos's exact dJ/dX, not only on field values — measured on a held-out set at the same seed and architecture, the gradient RMSE drops **4.8x** (0.376 → 0.079) and the value RMSE improves as a by-product. `TrainModel` gained the plumbing it needs: `extra_loss_terms` now resolves arity (a four-argument term receives the batch targets, three-argument terms are untouched) and `"target_channels"` keeps the data loss on the columns the model actually predicts — without it torch does not reject a 1-channel prediction against a 4-column target, it **broadcasts**, and trains against the mean of the value and the gradient.
    - **`SurrogateResponseFunction`** is the mirror of the co-simulation wrapper: a trained model where a *response function* goes rather than where a solver goes, with Kratos's own factory signature so any driver resolving responses by module path takes it unchanged. `"gradient_mode"` chooses autograd through the model (cheap, surrogate-accurate) or the FEM adjoint around the state it wrote (discretely exact *for that state* — a surrogate replaces the solve, not the sensitivity analysis).

- **Physics-informed residuals, PINN solves, exact residual losses and adjoints**: three residual notions.
    1. `solver_residuals.ResidualEvaluator` assembles the real PDE residual of any (ML-predicted) field through the solver's own builder machinery — the cheap non-differentiable *score* feeding query strategies, epoch callbacks and validation.
    2. `physics_informed`: SymPy-defined strong-form residuals evaluated by `physicsnemo.sym`'s `PhysicsInformer` (bundled — no extra install) as real **training loss terms** (`TrainModel(..., extra_loss_terms=[...])`; autodiff at point coordinates, least-squares on mesh/particle graphs, finite-difference/spectral on grids), with builtin diffusion, convection-diffusion, **linear elasticity** and **incompressible Navier–Stokes** PDEs (vector fields auto-split into per-component sympy functions), plus a `PinnSolveProcess` for pure-PINN forward solves from the model part's Dirichlet data and inverse coefficient recovery.
    3. `differentiable_residual`: the **exact discrete residual through the real FEM assembly** as a `torch.autograd.Function` — forward = `BuildRHS`, backward = the consistent tangent's transpose (a matvec, no solve), gradcheck-pinned on real ConvectionDiffusion and StructuralMechanics cases — with `MakeExactResidualLossTerm` making the physics' own verdict a gradient-carrying loss, now covering **transients** too (element-integrated time stepping unchanged; Bossak/Newmark displacement schemes via `scheme=` plus a per-step `InitializeSolutionStep`, assembling the effective tangent K + M(1-alpha)c0 + D c1).

    Alongside them, `sensitivity_utils` provides cheap surrogate `dJ/dx` (autograd through any point-cloud interface) plus **exact adjoint parameter sensitivities** (`Kᵀλ = ∂J/∂u`, one solve for all parameters, validated against full finite differences through real solves).

- **Training utilities and model governance**: a `Parameters`-driven training loop (`TrainModel`) and checkpoint saving (`SaveTrainedModel`, physicsnemo `.mdlus` or TorchScript) with optional **model cards** — JSON sidecars describing a checkpoint's fields, validated (advisorily) by every deployment process — plus autoregressive **rollout evaluation** (`EvaluateRollout`) exposing multi-step error growth of time-series surrogates. A `"performance"` block adds AMP, physicsnemo's `StaticCaptureTraining`, gradient clipping, the `Profiler` and `LaunchLogger`, learning-rate schedules, resumable checkpoints and `"optimizer": "muon"` — every one of them off by default, an empty block reproducing the previous history exactly, and a run resumed after two of four epochs returning the two epochs the uninterrupted run produced. **Self-supervised geometry pretraining** (`aerojepa_pretraining`) learns shape embeddings from a family of SDF-generated surfaces before any solve exists, the pretext task being signed-distance regression because upstream ships AeroJEPA's forward but neither a JEPA loss nor an EMA target encoder.

- **Validation, mesh datapipes and distributed alignment**: a `ValidationMetricsProcess` benchmarking predictions against reference fields with `physicsnemo.metrics` (JSON reports) — including `relative_mse` and the **distribution** metrics `histogram_l1`/`entropy_difference` on shared bin edges, plus a `"spectral_comparisons"` block comparing azimuthally averaged **power spectra**, because a surrogate can match a field on average while collapsing its extremes: the test pins a smooth prediction of a rough field scoring an RMSE below 0.2 while carrying under a tenth of the reference's fine-scale energy, extended by `cfd_bridge` with delegation to the optional `physicsnemo-cfd` package — its domain-aware CFD metric registry (relative-L2, drag/lift, physics residuals, UQ metrics) via the process's `cfd_metrics` block, plus Kratos↔pyvista `Flowfield` conversion and hybrid-initialization blending; a `MeshExportProcess` producing `.pmsh` mesh series directly consumable by `physicsnemo.datapipes.mesh_dataset` for training mesh-based models; a `curator_bridge` turning any Kratos solve into a **`physicsnemo-curator` ETL source** (the sinks ship upstream) with a `CuratorExportProcess` writing AI-ready Zarr stores or VTU grids per step straight from the solution loop; a **digital-twin USD export** (`usd_export`/`UsdExportProcess`, optional `usd-core` — no Omniverse install needed to write) turning a running solve, deployed-surrogate predictions and uncertainty fields included, into one time-sampled OpenUSD stage any USD viewer (usdview, Omniverse, Blender) scrubs directly — points and `primvars:` per step, topology re-sampled only on the steps it actually changes so adaptive-remeshing series stay valid, a `UsdGeomPoints` mode for particle clouds, MPI-gathered like the other exporters; and `distributed_utils` aligning `physicsnemo.distributed.DistributedManager` with Kratos's `DataCommunicator`, with a loud consistency check when the two disagree.

- **Uncertainty and governance**: guardrails, calibrated error bars and probabilistic validation across every deployment process.
    - **Out-of-distribution guardrails** (`ood_guard_utils` bridging `physicsnemo.experimental.guardrails.embedded.OODGuard`) — calibrated on the training inputs by `TrainModel`, saved as a sidecar, checked per inference under an advisory/strict/ignore policy.
    - **Geometry guardrails** (`geometry_guard_utils`) — the same question asked of the *shape* rather than the field values: a density model over descriptors of the outward-oriented boundary surface the mesh bridge already builds, so a surrogate handed a geometry outside its training family says so. The input-value guard cannot see that at all. A `"geometry_guard"` block on `InferenceProcess` reaches every deployment process at once.
    - **Predictive uncertainty for any deployed model** via the `"uncertainty"` block: MC dropout over dropout-like layers including learnable `physicsnemo.nn.ConcreteDropout`, checkpoint ensembles, or a **GP head** giving calibrated posterior variance whose *epistemic* term rises where the surrogate extrapolates — fitted post-training by `FitGpHead`, which implements the inducing-point/auxiliary-MSE/KL-ramp recipe upstream requires but does not ship, and saved as a sidecar since gpytorch modules cannot be TorchScripted. Mean goes to the output fields, per-node std to dedicated uncertainty fields, generalizing the diffusion bridge's pattern.
    - **Calibration metrics** (coverage, NLL, sharpness) answering whether those error bars are honest at all, and **ensemble scoring** (`crps`/`kcrps` over explicitly named members, with `retain_ensemble` keeping the members an inference process would otherwise reduce away).
    - **Probabilistic validation**: `relative_l2`/`weighted_mse`/`weighted_rmse` comparisons, ensemble `crps`/`kcrps` via `physicsnemo.metrics`, and multi-step rollout UQ in `EvaluateRollout` (per-step metric curves and MC-dropout std growth).

- **Active learning**: Kratos as the ground-truth solver in a `physicsnemo.active_learning` loop, via a `LabelStrategy` implementation with two execution backends — an in-process backend for small problems and a subprocess backend (recommended) that keeps Kratos MPI ranks and `torch.distributed` ranks in separate OS processes, fans labeling batches out over `max_parallel_jobs` concurrent solves, and supports HPC job submission. Completed by a query-strategy library (ensemble disagreement, predictive entropy, and **solver-residual scoring** — the real PDE residual assembled by Kratos ranks the surrogate's weak spots), a `MetrologyStrategy` backed by the validation-metrics machinery, and optional strict model-card enforcement at deployment.

  ```python
  import KratosMultiphysics as Kratos
  from KratosMultiphysics.PhysicsNeMoApplication.active_learning.kratos_label_strategy import CreateKratosLabelStrategy
  from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.subprocess_backend import SubprocessBackend

  backend = SubprocessBackend(Kratos.Parameters("""{
      "template_directory": "my_case_template",
      "run_command": ["python3", "MainKratos.py"]
  }"""))
  label_strategy = CreateKratosLabelStrategy(backend, provides_fields={"VELOCITY__node_historical"})
  # pass label_strategy to physicsnemo.active_learning.config.StrategiesConfig
  ```

  The case template must run this application's `DatasetExportProcess` so results come back as `.npz` files (a `srun`/`sbatch --wait` prefix in `run_command` turns the subprocess backend into an HPC job submission).

## 🗺️ Roadmap:

Pending work, and nothing else — fifteen items across four categories, each gated on **hardware** (one GPU on the reference machine, no GPU CI runner), **upstream** (the API does not exist yet, raises `NotImplementedError`, or does not import), **external access** (a credential, service or download this environment lacks), **build** (a Kratos application the reference build does not compile), or **open** (nothing blocks it — it needs data, or someone's time). What has shipped is described under Features above and in the documentation; how items left this list, and what each round learned doing so, is kept on the [Retrospectives](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Retrospectives.html) page. Contributions and prioritization requests are welcome — please open an issue on the Kratos repository mentioning `PhysicsNeMoApplication`.

<details>
<summary>Full roadmap, by category</summary>

Every entry is verified live against **physicsnemo 2.2.0** — 2.2.1 is a fix-only release and the upstream 2.3.0 changelog is still empty, so 2.2 *is* the current surface — and against the reference build, rather than assumed. Items marked *experimental* live under `physicsnemo.experimental`, whose API may change between releases without notice.

### Verification gates on shipped machinery

| # | Item | Gate | What is missing, specifically |
|---|---|---|---|
| 1 | GPU CI for the device-dependent paths | hardware | the application builds and runs its torch-free tests in the Linux CI; every ML-dependent test self-skips there. Needs a GPU runner. |
| 2 | Multi-GPU (NCCL) validation of the shipped distributed paths — halo-partitioned graph training, FSDP2 checkpoint round-trips, `ShardTensor` domain parallelism | hardware | all three are asserted at np=2/3 over CPU/gloo; NCCL rejects two ranks on the single GPU here, so the accelerator transport is untested |
| 3 | The CUDA-only `ShardTensor` paths — ring attention, the sharded kNN/radius search, `scatter_tensor`, `sharding_shapes="infer"` | hardware | upstream runs them on CUDA streams (and `scatter_tensor` through the `DistributedManager`); none can run on the CPU mesh the shipped `domain_parallel_utils` builds |
| 4 | Upstream 3-D `fill_interior` | upstream | its `n = 3` raises `NotImplementedError` in 2.2 ("pending exact 3D boundary recovery"); the *need* is closed locally by the opt-in `"method": "tetgen"` backend, so this is about dropping that optional dependency |
| 5 | Live verification of the NIM client | external access | `nim_client`, `NimInferenceProcess` and the payload contract ship stub-pinned; running a NIM needs an NGC API key and docker. `TestNimLive` self-skips until `PHYSICSNEMO_NIM_URL` points at one |

### Recipes that ship test-pinned but lack data or a worked example

| # | Item | Gate | What is missing, specifically |
|---|---|---|---|
| 6 | Richer crash architectures | open | modeling work on the shipped structural-surrogate pattern; needs a crash dataset, not new machinery |
| 7 | Virtual Foundry / additive-manufacturing datasets | open | `vfgn_bridge` and the real sintering case ship; furnace and AM trajectories at scale do not |

### physicsnemo 2.2 capabilities not yet bridged

Read off the installed package's module tree and the 2.2.0 changelog; each signature below was checked with `inspect.signature` against 2.2.0.

| # | Item | Gate | What is missing, specifically |
|---|---|---|---|
| 8 | DPOT's pretrained weights | external access | the *architecture* ships - `"model_interface": "dpot"` on `SequenceInferenceProcess`, with LoRA fine-tuning through `domino_finetune`'s generic pattern - and a rollout beats persistence on Kratos grid series. The authors' pretrained checkpoints are not verified from here, so "foundation model" stays a claim about someone else's weights |
| 9 | Kratos HDF5 output through `HDF5Reader` | build | the VTK half shipped (`bridges/vtk_bridge`, training-only: a VTK file carries no provenance, so there is no scatter-back); `physicsnemo.datapipes.readers.HDF5Reader` over `HDF5Application` output is gated on that application, which the reference build does not compile |
| 10 | physicsnemo-cfd evaluation benchmarks and reports on Kratos data | upstream | `physicsnemo.cfd.evaluation.benchmarks` declares `run_benchmark`/`write_report` and `reports.register_visual`, but the installed 0.0.3a0 source clone raises `ImportError` on the import; revisit on the next physicsnemo-cfd release |
| 11 | Equivariant surrogates instead of augmentation | upstream, *experimental* | attempted and refused on measured evidence rather than left open: `physicsnemo.experimental.nn.symmetry`'s `rotate_grid_coefficients` does not act as an SO(3) rotation in this environment - with the grid mask applied correctly it moves the l=0 **scalar** part (2.6), fails to preserve the l=1 norm (7.0), round-trips under neither inverse convention and does not commute with a z-rotation. A rotation that moves the scalar part is not a rotation, and four probes did not find the convention that makes it one. Needs an upstream rotation helper that demonstrably is a group action, or a worked upstream example to copy |
| 12 | HealDA-style data assimilation for transient twins | upstream | `physicsnemo.experimental.models.healda.VideoHealDA` maps sparse observations plus a prior state to an assimilated field, which is what a sensor-fed digital twin wants - but its API is HEALPix-pixel and calendar-shaped (`npix`, `second_of_day`, `day_of_year`). Watch item until a domain-agnostic variant exists |

### Kratos-side cases still missing

What is missing here is on the Kratos side - a cross-validation, a driver loop, or an application the reference build does not compile - not anything upstream.

| # | Item | Gate | What is missing, specifically |
|---|---|---|---|
| 13 | Adjoint cross-validation of the aerodynamic surrogate | open | the surrogate itself ships on the compiled `CompressiblePotentialFlowApplication` (a Mach sweep to airfoil-skin pressure, held-out MSE 3.4e-4 against the mean predictor's 4.7e-2), and `adjoint_bridge`'s `compressible_potential_flow` dispatch is pinned as resolving to the real factory. A *full* adjoint needs the separate adjoint solver and primal-restart files the owning application sets up for its own regression, which no case here builds yet |
| 14 | A closed-loop shape optimizer with cross-validation | build | `ComputeControlSensitivities` + `MeshMovingApplication` smoothing are shipped; a driver loop cross-validated against `ShapeOptimizationApplication`/`OptimizationApplication` needs those two compiled (the committed vertex-morphing fixture stands in for now) |
| 15 | Particle applications behind `ParticleInferenceProcess` | build | the process reads any particle node cloud; `MPMApplication`, `SPHApplication`, `DEMApplication` and `PfemFluidDynamicsApplication` are not compiled here, so the tests drive synthetic clouds |

</details>

## ⚠️ Dependency policy — read before contributing:

`torch` and `nvidia-physicsnemo` (and the extras: `nvidia-physicsnemo-cfd` with `pyvista`, `torch_geometric`, `cupy`, `usd-core`, `tetgen`, ...) are **optional, pure-Python runtime dependencies**:

- `import KratosMultiphysics.PhysicsNeMoApplication` always succeeds without them; only the specific submodules that need them fail — lazily, at call time, with an actionable error message.
- There is intentionally **no CMake gating** (no `find_package(Torch)`, no `USE_PHYSICSNEMO` option): nothing in the C++ core links against torch or CUDA. Do not add such a gate.
- New Python modules must keep all `import torch` / `import physicsnemo` / `import cupy` statements inside lazy `_TryImport*()` helpers, never at module scope of anything imported eagerly.

**Packaging policy** (wheels / the standard `KratosMultiphysics` distribution): the application ships like any other per-application Kratos wheel — `scripts/wheels/build_wheel.py` derives the app wheels automatically from whatever was compiled, so no wheel-side registration exists or is needed. `torch` / `nvidia-physicsnemo` are **never** declared as wheel dependencies (they are large, CUDA-variant-specific, and the user's choice of build); they remain the optional runtime `pip install`s above, and the lazy-import contract (`tests/test_import_contract.py`) is exactly what guarantees a wheel installed without them stays fully importable.

## ⚙️ Installation:

Compile like any other Kratos application (add `add_app ${KRATOS_APP_DIR}/PhysicsNeMoApplication` to your configure script). To use the ML features additionally install the runtime dependencies:

```bash
pip install torch             # CPU or CUDA build, your choice
pip install nvidia-physicsnemo
pip install git+https://github.com/NVIDIA/physicsnemo-cfd   # optional: hybrid-initialization recipes + CFD metric registry (source only - not on PyPI)
pip install cupy-cuda12x     # optional: the GPU array paths (match your CUDA toolkit)
pip install usd-core         # optional: the digital-twin USD export (UsdExportProcess) - no Omniverse install needed
pip install tetgen           # optional: exact boundary recovery in the tetrahedral fill ("method": "tetgen"; note TetGen itself is AGPL-licensed)
```

<details>
<summary>Before adding applications to an existing build — a trap worth knowing</summary>

Kratos reads `KRATOS_APPLICATIONS` from the *environment* and never writes it to the CMake cache, and `make install` never prunes the install tree. An install directory therefore accumulates every application ever built, while the build graph covers only the ones in the current configure script — and those two sets drift apart silently. Applications outside the graph are not rebuilt, so the next `libKratosCore.so` relink leaves them ABI-stale: they keep *existing*, and `CheckIfApplicationsAvailable` keeps reporting them as present, but importing them fails with an undefined symbol or segfaults. Because several Kratos applications import their optional companions on a *presence* check rather than a working one, one stale library can take down an application that was never touched.

So: reconfiguring to add an application is not a local operation. Restore the **full** historical application list in the same run, and treat the install tree as unreproducible from the cache alone — nothing records what produced it.

</details>

## ⚙️ Examples:

Twenty example notebooks covering the tensor bridge, mesh bridge, surrogate training/deployment, active learning, superresolution, transient rollouts, exact shape gradients, DoMINO fine-tuning, the adjoint integration and molecular dynamics can be found in [`examples/notebooks/`](examples/notebooks/).

They are **executed** by `tests/test_notebooks.py`, one test per notebook, so a
changed signature or return type breaks a test rather than rotting silently. Each
runs in a throwaway copy of the tree, because the notebooks write artifacts next to themselves and several of them load cases from `tests/kratos_solver_cases` and `tests/adjoint_cases`. At roughly five minutes for all twenty they sit on the **validation** suite, so they run only when that suite is asked for — not in CI, and not nightly.

Executing them proves no exception was raised, which is weaker than proving the right answer: the seven notebooks making a reproducible numeric claim (03, 05, 07, 17, 18, 19 and 20) assert it, and the unseeded ones deliberately assert nothing.

### Worked examples

Twenty-six fully documented use cases live in the [Examples repository](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application): self-contained scripts against real Kratos solves, with every figure regenerated by the code — the five newest being [Lennard-Jones molecular dynamics](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/lennard_jones_md) (a learned force field on periodic radius graphs, rolled out by `ParticleInferenceProcess`, with a log-axis plot of exactly where the self-driven rollout stops working rather than a claim that it does not), [Generative topology design](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/generative_topology_design) (Kratos in the loop twice — an optimality-criteria loop computes the training optima, and every design TopoDiff proposes is re-solved to score it), [Residual-guided diffusion](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/residual_guided_diffusion) (one trained denoiser steered at sampling time toward ten sparse sensors and toward the exact discrete FEM residual, with no retraining, and the `std_y` sweep that shows where guidance is silently off), [Geometry guardrail](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/geometry_guardrail) (a slab that normalizes to look exactly like the training cubes, so the field guard stays silent while the geometry guard rejects it at percentile 100) and [Boundary-to-interior operators](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/boundary_to_interior_operators) (GLOBE reading the actual boundary mesh and xDeepONet reading the case parameters, both beating the input-ignoring predictor on every held-out case, with no POD basis anywhere). Highlights (each image links to its case):

<details>
<summary>All highlights (thirteen cases, one image each)</summary>

**[Generative topology design](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/generative_topology_design)** — the solver computes the training set by running real topology optimization, and scores every proposed design by re-solving it; the sampled designs beat volume-matched random controls by 2.6x to 5.0x and stay clearly worse than the optimizer, which the case reports rather than hides:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/generative_topology_design"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/generative_topology_design/data/generated_designs.png" alt="Optimality-criteria optima, volume-matched random controls and sampled designs, every compliance on a log scale." width="560"/></a>
</p>

**[Thermal surrogate lifecycle](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/thermal_surrogate_lifecycle)** — dataset from real solves, FNO ensemble with model cards and an OOD guard, governed in-loop deployment, validation, ONNX distillation and a Triton serving layout:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/thermal_surrogate_lifecycle"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/thermal_surrogate_lifecycle/data/deployment_fields.png" alt="Solver vs FNO ensemble, error vs predicted error bar." width="560"/></a>
</p>

**[Transient thermo-mechanical surrogate](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/transient_thermomechanical_surrogate)** — coupled cooling-contraction ground truth, single-step vs BPTT training, self-fed rollout on an unseen cooling schedule:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/transient_thermomechanical_surrogate"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/transient_thermomechanical_surrogate/data/rollout.gif" alt="Coupled solver vs surrogate rollout." width="620"/></a>
</p>

**[Adjoint integration](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/adjoint_integration)** — Kratos's own adjoint stack read through `adjoint_bridge` on two physics (agreeing with the shipped adjoint to 7e-7 and 8e-9), a surrogate trained on those gradients (46x better `dJ/dθ` at held-out designs), and the descent that spends the difference — plus the model deployed back as a Kratos response function:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/adjoint_integration"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/adjoint_integration/data/gradient_enhanced_surrogate.png" alt="Training curves and the surrogate's dJ/dtheta against the exact adjoint." width="700"/></a>
</p>

**[GNN + exact shape optimization](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/gnn_and_exact_shape_optimization)** — MeshGraphNet on the mesh's own graph, exact adjoint shape gradients verified against re-solved finite differences to ~1e-10, an FFD optimization hitting its target exactly:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/gnn_and_exact_shape_optimization"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/gnn_and_exact_shape_optimization/data/shape_optimization.png" alt="Objective convergence and optimized domain." width="700"/></a>
</p>

**[Superresolution + CorrDiff diffusion](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/superresolution_and_diffusion)** — two-stage diffusion downscaling of genuinely under-resolved solves, with calibrated ensemble uncertainty:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/superresolution_and_diffusion"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/superresolution_and_diffusion/data/corrdiff_downscaling.png" alt="Coarse condition, regression mean, truth and ensemble spread." width="760"/></a>
</p>

**[Transient in-loop superresolution](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/transient_superresolution)** — the process attached to a running transient analysis, upscaling every step; the learned upsampler halves bilinear interpolation's error:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/transient_superresolution"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/transient_superresolution/data/transient_sr.gif" alt="Coarse, superresolved and fine transient fields." width="700"/></a>
</p>

**[Physics-refined superresolution](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/physics_informed_superresolution)** — the differentiable PDE residual grades the upsampler's output and refines it below the fine solve's own residual:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/physics_informed_superresolution"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/physics_informed_superresolution/data/pisr_metrics.png" alt="RMSE and PDE residual per field." width="640"/></a>
</p>

**[Lid-driven cavity MeshGraphNet](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/fluid_cavity_gnn)** — the first fluid case: Navier–Stokes cavity solves, a GNN propagating boundary conditions through the mesh graph to the recirculating flow it never sees in its input:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/fluid_cavity_gnn"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/fluid_cavity_gnn/data/cavity_fields.png" alt="VMS solver vs MeshGraphNet cavity flow." width="760"/></a>
</p>

**[Surrogate-driven adaptive remeshing](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/surrogate_driven_remeshing)** — the residual-scoring → MMG loop: elements concentrate exactly where the surrogate's physics error lives, with no reference solution needed:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/surrogate_driven_remeshing"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/surrogate_driven_remeshing/data/remeshing.png" alt="Surrogate error, driving residual, adapted mesh." width="760"/></a>
</p>

**[Lagrangian particle surrogate](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/lagrangian_particle_surrogate)** — meshless: a learned acceleration integrated into particle positions by the process, de-normalized by the model card:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/lagrangian_particle_surrogate"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/lagrangian_particle_surrogate/data/particles.gif" alt="Falling particle cloud, surrogate vs closed form." width="480"/></a>
</p>

**[PINN forward + inverse](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/pinn_forward_and_inverse)** — mesh-free solves against the FEM reference and inverse conductivity recovery from field observations:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/pinn_forward_and_inverse"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/pinn_forward_and_inverse/data/pinn_inverse.png" alt="PINN training and recovered conductivity." width="640"/></a>
</p>

**[Animated adaptive remeshing](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/animated_adaptive_remeshing)** — the mesh chases a moving heat source (surrogate-residual MMG adaptation), recorded by core Kratos' new `TransientPlotter` with the topology changing every frame:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/animated_adaptive_remeshing"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/animated_adaptive_remeshing/data/chasing_refinement.gif" alt="The adapted mesh following the moving source." width="520"/></a>
</p>

**[The sintering solve, recorded live](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/transient_thermomechanical_surrogate)** — the coupled thermo-mechanical transient rendered by `PyVistaAnimationOutputProcess` attached to the running analysis (temperature-colored, displacement-warped):

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/transient_thermomechanical_surrogate"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/transient_thermomechanical_surrogate/data/sintering_solve.gif" alt="The cooling, shrinking body recorded by the animation output process." width="440"/></a>
</p>

**[Co-simulation surrogate](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/cosim_surrogate)** — a trained checkpoint as a first-class CoSimulation solver, Aitken-accelerated to the learned operator's fixed point; **[uncertainty and trust](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/uncertainty_and_trust)** — spread vs actual error, calibration catching over-confident bars, the OOD guard; **[active learning](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/active_learning)** — Kratos as the labeler, with its honest result:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/uncertainty_and_trust"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/uncertainty_and_trust/data/uncertainty_sweep.png" alt="Ensemble spread vs actual error across the sweep." width="620"/></a>
</p>


**[Hybrid initialization](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/hybrid_initialization)** — the surrogate accelerating the solver itself: `HybridInitializationProcess` warm-starts Newton–Raphson on a geometrically nonlinear cantilever (25 % of the iterations saved, identical solutions), with the measured lesson that Newton rewards *smooth* seed accuracy — a rough seed closer to the solution loses to a cold start:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/hybrid_initialization"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/hybrid_initialization/data/newton_acceleration.png" alt="Newton iterations cold vs warm-started, and residual traces." width="640"/></a>
</p>

**[Implicit geometry + SDF surrogates](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/implicit_geometry_sdf_surrogate)** — meshes generated straight from signed distance functions (`mesh_bridge.generate`), solved by ConvectionDiffusion, and learned through a ladder of autograd-computed SDF encodings that generalizes across shapes:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/implicit_geometry_sdf_surrogate"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/implicit_geometry_sdf_surrogate/data/generated_family.png" alt="Six SDF-generated geometries with their solved fields." width="620"/></a>
</p>

**[Sparse-sensor field inversion](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/sparse_sensor_inversion)** — the FWI-style recipe: a mask-conditioned diffusion model reconstructs the full temperature field from ten sensor readings, ensemble std as the inversion uncertainty; **[ROM-space temporal surrogate](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/rom_temporal_surrogate)** — RomApplication's `CalculateRomBasisOutputProcess` writes the POD basis, physicsnemo's temporal attention learns the dynamics in eight coefficients:

<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/sparse_sensor_inversion"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/sparse_sensor_inversion/data/sensor_inversion.png" alt="Field inversion from ten sensors with ensemble uncertainty." width="700"/></a>
</p>
<p align="center">
  <a href="https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/rom_temporal_surrogate"><img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/physics_nemo_application/use_cases/rom_temporal_surrogate/data/coefficient_rollout.png" alt="Reduced coordinates, truth vs rollout at an unseen conductivity." width="700"/></a>
</p>

</details>

## 🗎 Documentation:

User documentation is available in the [Kratos documentation pages](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Overview.html) (sources under `docs/pages/Applications/PhysicsNeMo_Application/`). Further information regarding NVIDIA PhysicsNeMo can be found in the [official documentation](https://docs.nvidia.com/physicsnemo/latest/index.html).