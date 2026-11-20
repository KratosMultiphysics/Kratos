# `python_scripts/` — what is where

A module's folder says what kind of thing it is. If you are looking for
something, start here; the same map with an *I want to X → use Y* index is in
the [documentation](https://kratosmultiphysics.github.io/Kratos/pages/Applications/PhysicsNeMo_Application/General/Module_Map.html).

**The rule:** everything under `processes/` has a `Factory`, so it can be
attached from `ProjectParameters.json`. Nothing outside `processes/` does.

Each package's `__init__.py` carries a docstring saying what belongs in it.
Read it before adding a module — this directory was flat once.

```
processes/                          attach to a solve
├── inference/                      run a trained model in the solution loop
│   ├── inference_process           flat per-node/element feature vectors
│   ├── graph_inference_process     the element-edge graph
│   ├── grid_inference_process      a regular voxel grid
│   ├── point_cloud_inference_process   the nodes as an unordered cloud
│   ├── sequence_inference_process  a rolling window of grid states
│   ├── time_series_inference_process   a rolling window of nodal states
│   ├── particle_inference_process  a proximity graph over particles
│   ├── diffusion_inference_process an ensemble, with uncertainty fields
│   ├── domino_inference_process    a preprocessed DoMINO datapipe case
│   ├── onnx_inference_process      ONNX Runtime instead of torch
│   ├── triton_inference_process    a remote Triton server
│   ├── nim_inference_process       a running PhysicsNeMo NIM microservice
│   ├── superresolution_process     coarse grid in, fine grid out
│   ├── rom_surrogate_process       parameters in, u = phi q out
│   ├── hybrid_initialization_process   seeds the solver instead of replacing it
│   └── pinn_solve_process          no model: it *is* the solve
├── export/                         write solver data out as training data
│   ├── dataset_export_process      .npz per step
│   ├── grid_dataset_export_process the same, on a voxel grid
│   ├── cae_dataset_export_process  physicsnemo.datapipes.cae layout
│   ├── mesh_export_process         .pmsh mesh series
│   ├── curator_export_process      AI-ready Zarr / VTU
│   ├── usd_export_process          a time-sampled OpenUSD digital twin
│   └── streaming_dataset_export_process   a live queue, no files
├── adaptive_remesh_process         changes the mesh
├── adjoint_sensitivity_process     dJ/dX onto the model part, exporters carry it
└── validation_metrics_process      measures the result

bridges/                            Kratos data <-> physicsnemo data
├── torch_bridge                    fields <-> torch.Tensor, zero-copy
├── mesh_bridge/                    meshes (see its own __init__)
│   ├── tessellation                simplices, watertight across neighbours
│   ├── curved_tessellation         the isoparametric mode
│   ├── provenance                  predictions back onto the original entities
│   ├── domain_mesh_builder         DomainMesh with named boundaries
│   ├── generate                    geometry from implicit functions
│   ├── nurbs_sampling              exact NURBS (IGA) geometry on a lattice
│   ├── spatial                     signed distance fields as features
│   ├── deformation                 differentiable shape parameterizations
│   ├── operations                  curvature, repair, subdivision, extrusion
│   ├── sampling                    differentiable point sampling over a BVH
│   └── adaptive_remeshing          residual-driven MMG adaptation
├── graph_bridge                    the element-edge graph
├── grid_bridge                     unstructured fields <-> voxel grids
├── particle_bridge                 radius / kNN proximity graphs
├── mapping_bridge                  non-matching transfer via MappingApplication
├── calculus_bridge                 gradient, divergence, curl, Laplacian
├── rom_bridge                      RomApplication POD bases
├── adjoint_bridge                  Kratos adjoint gradients in row order
├── cfd_bridge                      pyvista and physicsnemo-cfd
├── globe_bridge                    boundary data -> the interior (GLOBE)
├── curator_bridge                  a solve as a physicsnemo-curator source
├── vtk_bridge                      Kratos VTK output as read-only training data
└── vfgn_bridge                     Virtual Foundry GraphNet, sintering / AM

training/                           loops, datasets, schemes
├── training_utils                  TrainModel, SaveTrainedModel, ExportOnnxModel
├── torch_dataset                   dataset and datapipe factories
├── streaming_dataset               train out of a running solve
├── temporal_training               window datasets, BPTT through a rollout
├── diffusion_utils                 diffusion, DPS guidance, multi-diffusion, CorrDiff
├── globe_training                  the boundary-driven operator's own loop
├── aerojepa_pretraining            geometry pretraining before any solve
├── sobolev_training                grade the surrogate on exact gradients too
├── domino_finetune                 predictor-corrector and LoRA adaptation
├── rom_temporal                    temporal attention in ROM space
└── rollout_utils                   multi-step error growth

physics/                            physics as a signal
├── solver_residuals                the real PDE residual - a score, not a gradient
├── physics_informed                SymPy strong-form residuals as a loss
├── differentiable_residual         the exact discrete residual, differentiable
├── diffusion_residual_operator     that residual as a diffusion guidance operator
└── sensitivity_utils               adjoints and exact shape gradients

deployment/                         checkpoint -> production
├── model_registry                  loading, model cards, de-normalization
├── onnx_utils                      the ONNX Runtime session and its devices
├── triton_export                   a Triton model repository
├── nim_client                      the documented physics-NIM HTTP contract
├── usd_export                      time-sampled OpenUSD stages (digital twins)
├── cosim_surrogate_solver_wrapper  a model as a CoSimulation solver
├── surrogate_response_function     a model as a Kratos response function
├── uncertainty_utils               MC dropout, ensembles, GP heads
├── geometry_guard_utils            OOD detection on the shape
└── ood_guard_utils                 out-of-distribution guardrails

distributed/                        MPI and multi-rank
├── distributed_utils               DistributedManager <-> DataCommunicator
├── graph_partition_utils           halo-partitioned graph training (data parallel)
└── domain_parallel_utils           ShardTensor over the Kratos ranks (domain parallel)

active_learning/                    Kratos as the labeling oracle

utilities/                          small shared helpers
├── tensor_adaptor_dataset_utils    the shared gather/scatter entry point
├── array_backend_utils             opt-in CuPy, with numpy the default
├── point_subsampling               point-cloud token budgets
├── nvtx_utils                      Nsight Systems ranges
├── shallow_water_reference         a numpy-only reference integrator (GraphCast recipe)
└── lennard_jones_reference         a numpy-only MD integrator (Lennard-Jones recipe)
```

## Adding a module

1. Pick the folder from the rule above; if none fits, that is worth a
   conversation before it becomes a tenth package.
2. Keep every `import torch`, `import physicsnemo` and `import cupy` inside a
   lazy `_TryImport*()` helper — `tests/test_import_contract.py` enforces it.
3. Add its tests to `tests/` and register the classes in
   `tests/test_PhysicsNeMoApplication.py`; `TestSuiteRegistration` fails
   otherwise.
4. Documented module paths are resolved against this tree by
   `TestDocumentedImportPathsResolve`, so a move that leaves prose behind fails
   the suite.
