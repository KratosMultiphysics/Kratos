# PhysicsNeMoApplication example notebooks

A learning path through the application, in order. Every notebook runs on CPU in a few
seconds to minutes with small analytic problems (07 and later use real solves), and writes its outputs to a local `./output/`
folder.

Every notebook is committed **with its outputs executed**, and every one renders what it
solves - the model part, the field a process wrote, the error - with pyvista through the
core `KratosMultiphysics.pyvista_utilities` bridge (mesh-aware renders, or point renders
for node-only parts), so the mesh and the fields are visible before any number is read.

| Notebook | What it shows |
|---|---|
| [01 — Tensor bridge and dataset export](01_tensor_bridge_and_dataset_export.ipynb) | `ModelPart` ↔ `torch.Tensor` (zero-copy in, checked copy out), `DatasetExportProcess`, `CreateNpzDataset` + `DataLoader` |
| [02 — Mesh bridge round trip](02_mesh_bridge_round_trip.ipynb) | Tessellation of a mixed hex+tet mesh with provenance, volume conservation, `DomainMesh` with a named boundary, exact nodal scatter-back, native save/load — both meshes rendered with the core `pyvista_utilities` bridge |
| [03 — Training and deploying a surrogate](03_training_a_surrogate.ipynb) | Parameter sweep → dataset export → MLP training → TorchScript checkpoint → `InferenceProcess` and `HybridInitializationProcess` deployment |
| [04 — Active learning with Kratos as the oracle](04_active_learning_with_kratos.ipynb) | `KratosALSample`, `InProcessBackend`, `CreateKratosLabelStrategy`, one real `physicsnemo.active_learning.Driver` step, and how `SubprocessBackend` scales this to HPC |
| [05 — Superresolution](05_superresolution.ipynb) | Coarse/fine meshes → voxel grids via `grid_bridge`, training `physicsnemo.models.srrn.SRResNet`, `.mdlus` checkpointing, and deployment with `SuperResolutionProcess` (slice comparison at an unseen parameter) |
| [06 — Mesh graph nets and transient surrogates](06_mesh_graph_nets_and_transient_surrogates.ipynb) | `graph_bridge` element-edge graph extraction, training `MeshGraphNet` directly on the mesh + `GraphInferenceProcess` (needs `torch_geometric`/`torch_scatter`), and an autoregressive next-state surrogate rolled forward by `TimeSeriesInferenceProcess` |
| [07 — Real-solver thermal surrogate](07_real_solver_thermal_surrogate.ipynb) | The Darcy-FNO recipe on **real Kratos solves**: `ConvectionDiffusionApplication` conductivity sweep, FNO training on `grid_bridge` grids, deployment with `GridInferenceProcess`, `ValidationMetricsProcess` numbers, mesh-aware field comparison via the core `pyvista_utilities` bridge (availability-gated on the compiled applications) |
| [08 — RNN grid sequences](08_rnn_grid_sequences.ipynb) | Transient field → grid series (`GridDatasetExportProcess`), `CreateGridSequenceDataset` + `One2ManyRNN` training, and seed-once/roll-forward deployment with `SequenceInferenceProcess` |
| [09 — Diffusion downscaling](09_diffusion_downscaling.ipynb) | CorrDiff-style conditional diffusion: matched (blurred, sharp) grid pairs, `TrainDiffusionModel` over `EDMPrecondSuperResolution`/`SongUNet`, and `DiffusionInferenceProcess` writing ensemble mean + per-node uncertainty |
| [10 — Structural deformation surrogate](10_structural_deformation_surrogate.ipynb) | The crash/deformation-surrogate pattern on **real StructuralMechanics solves**: cantilever load sweep, load → displacement-field surrogate via `TrainModel`, `InferenceProcess` deployment at an unseen load, warped-mesh comparison via the core `pyvista_utilities` bridge (availability-gated) |
| [11 — ROM surrogate and temporal attention](11_rom_surrogate_and_temporal_attention.ipynb) | Neural-augmented reduced bases on **real RomApplication POD output** (numpy-SVD fallback shown): `rom_bridge`, conductivity → modal-coefficients surrogate deployed by `RomSurrogateProcess`, reconstructed-field comparison via the core `pyvista_utilities` bridge, and `Sequence_Model` temporal attention over ROM trajectories |
| [12 — Co-simulation surrogate coupling](12_cosim_surrogate_coupling.ipynb) | A trained surrogate as a **solver** in Kratos co-simulation: `cosim_surrogate_solver_wrapper` referenced by module path, weak Gauss-Seidel coupling through `kratos_mapping` transfer, and an Aitken-accelerated strong-coupling fixed point (needs CoSimulationApplication + MappingApplication) |
| [13 — Lagrangian particle surrogate](13_lagrangian_particle_surrogate.ipynb) | Learning-to-Simulate on a Kratos particle cloud: `CreateParticleTrajectoryDataset` windows, training on velocity-history features, autoregressive deployment with `ParticleInferenceProcess` (proximity graph rebuilt every step, semi-implicit Euler advancing the nodes) |
| [14 — PINN forward and inverse solves](14_pinn_forward_and_inverse.ipynb) | `PinnSolveProcess` on `physicsnemo.sym` (bundled): a data-free Laplace solve from the model part's Dirichlet fixities alone, and inverse diffusion-coefficient recovery from observations (trainable PDE coefficient, detached field gradients) |
| [15 — CorrDiff two-stage downscaling](15_corrdiff_two_stage_downscaling.ipynb) | The regression + residual-diffusion split: `TrainCorrDiffPair` (upstream `RegressionLoss`/`ResidualLoss`, `SongUNetPosEmbd` sizing) and `DiffusionInferenceProcess` with `"regression_settings"` adding the regression mean to the ensemble — sharper mean, calibrated spread |
| [16 — Transient structural surrogate](16_transient_structural_surrogate.ipynb) | Temporal schemes on a real implicit-dynamic (Bossak) cantilever solve: `RunTransientAnalysis` trajectories, single-step window training, `TrainAutoregressive` BPTT with gradient checkpointing (2.7x lower rollout error), and deployment through `TimeSeriesInferenceProcess` |
| [17 — Exact shape gradients and shape optimization](17_exact_shape_gradients.ipynb) | The discretely exact `dJ/dX` of a real FEM solve from one element-local pass (`ComputeShapeSensitivityField`), the measured crossover against the per-parameter global path, `ComputeControlSensitivities` chaining it onto an FFD lattice, and a gradient descent driving the shape onto a target — validated against re-solve finite differences to ten digits |
| [18 — Fine-tuning a pretrained DoMINO](18_domino_finetuning.ipynb) | The public `nvidia/domino_drivaerml` checkpoint loaded through `model_registry`, the de-normalization a pretrained DoMINO needs, and both adaptation recipes of `domino_finetune`: the predictor-corrector decomposition with the frozen predictor's output cached, and LoRA adapters merged back into an ordinary `.mdlus` |
| [19 — Adjoint integration](19_adjoint_integration.ipynb) | Kratos's own adjoint stack read through `adjoint_bridge` (row order by id, never by iteration), `AdjointSensitivityProcess` putting dJ/dX on the model part so every exporter carries it, Sobolev training on the exact gradients, and `SurrogateResponseFunction` deploying the model where a Kratos response goes |
| [20 — Molecular dynamics with a MeshGraphNet](20_lennard_jones_molecular_dynamics.ipynb) | The Lennard-Jones recipe: a velocity-Verlet reference trajectory whose central-difference targets equal its forces to round-off, `particle_bridge`'s minimum-image radius graph in a **periodic box** (`"box_size"`), `MeshGraphNet` heads for per-atom force and potential energy, and `ParticleInferenceProcess` deployment asserting the model beats the mean-force predictor |

## Environment

```bash
export PYTHONPATH=/path/to/Kratos/bin/Release
export LD_LIBRARY_PATH=/path/to/Kratos/bin/Release/libs:$LD_LIBRARY_PATH
pip install torch nvidia-physicsnemo   # optional runtime dependencies
pip install jupyterlab matplotlib      # to run the notebooks
pip install pyvista                    # for the renders in 02, 07, 10, 11 (KratosMultiphysics.pyvista_utilities) and the point clouds of 13 and 20
```
