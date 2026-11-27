---
title: Glossary
keywords: glossary terms definitions surrogate neural operator datapipe mdlus pmsh provenance model card halo shardtensor crps epistemic
tags: [Glossary.md]
sidebar: physicsnemo_application
summary: The terms used across this documentation, from both sides of the bridge, each in a sentence or two with a pointer to where it matters.
---

# Glossary

**Adjoint.** The solve of the transposed tangent system, `K^T lambda = dJ/du`, which gives the sensitivity of an objective to every parameter or node coordinate at the cost of one solve. `physics.sensitivity_utils` ships one; `bridges.adjoint_bridge` reads Kratos's own. See [Adjoint](../Adjoint/Adjoint.html).

**Aleatoric versus epistemic.** Two components of predictive uncertainty: noise the data itself carries (aleatoric) and ignorance more data would remove (epistemic). The epistemic part is the out-of-distribution signal. See [Uncertainty and guardrails](../PhysicsNeMo_Basics/Uncertainty_And_Guardrails.html).

**BPTT.** Backpropagation through time: training an autoregressive surrogate by unrolling its own rollout and differentiating through it. `training.temporal_training` implements it with gradient checkpointing.

**Card, model card.** The JSON sidecar next to a checkpoint naming its input and output fields, their order and their normalization. Every deployment process validates against it and de-normalizes through it. See [Core and checkpoints](../PhysicsNeMo_Basics/Core_And_Checkpoints.html).

**Codimension.** Spatial dimension minus manifold dimension of a mesh: a surface in 3-D has codimension 1, a volume codimension 0. Some upstream operators are only valid at one codimension.

**Condition.** A Kratos boundary entity - a surface or line element carrying a boundary condition. Conditions of a sub-model-part become a named boundary mesh.

**CorrDiff.** The two-stage diffusion recipe: a regression model predicts the conditional mean, a residual diffusion model learns what is left. See [Diffusion](../Diffusion/Diffusion.html).

**CRPS.** Continuous ranked probability score: a proper scoring rule grading a whole ensemble against one observation, rewarding both accuracy and honest spread. Needs the members, not just their mean and standard deviation.

**Data location.** The string naming where a variable's values live on a Kratos entity: `node_historical`, `node_non_historical`, `element`, `condition`, `element_gauss_point`, `condition_gauss_point`. See [Kratos concepts](Kratos_Concepts_For_ML.html).

**Datapipe.** PhysicsNeMo's readers, transforms and datasets that turn files into batched tensors. See [Data and datapipes](../PhysicsNeMo_Basics/Data_And_Datapipes.html).

**De-normalization.** Inverting the scaling a model was trained under before writing its output onto a physical variable. Omitting it produces finite, plausible, wrong numbers - the most frequently rediscovered bug in this application's history.

**Denoiser, preconditioner, sampler.** The three replaceable parts of a diffusion model in `physicsnemo.diffusion`: the network, the noise scaling and conditioning, and the walk from noise to a sample.

**DomainMesh.** A PhysicsNeMo `Mesh` plus named boundary meshes plus global data - the natural image of a Kratos model part with sub-model-parts.

**DPS guidance.** Diffusion posterior sampling: steering a diffusion sampler toward observations or a forward model at sampling time, without retraining. Shipped as the `"guidance"` block of `sampler_settings`, whose `"kratos_residual"` operator is the exact discrete FEM residual.

**EDM.** The "elucidating the design space" diffusion formulation whose preconditioner and losses the bridge uses; its `sigma_data` sets the scale fields must be brought to.

**Execution point.** Which hook of the solution loop a process runs in - `initialize_solution_step` or `finalize_solution_step`. See [Architecture](Architecture.html).

**Ensemble.** Several trained models (or checkpoints) whose disagreement estimates uncertainty. Deep, snapshot, checkpoint and input ensembles differ in how the members are obtained.

**FSDP2.** Torch's fully sharded data parallelism: parameters split across ranks as `DTensor`s. `SaveTrainedModel` gathers them back into one checkpoint.

**Gauss point.** An integration point inside an element. Kratos can hand out values computed there but cannot write onto them, so Gauss-point fields are read-only for this application.

**Halo.** The layer of neighbouring entities a rank needs from other ranks so that a local operator (a convolution, a message-passing step) gives the same answer as the serial run. Kratos has a node halo (ghost nodes) but no element halo.

**Historical variable.** Nodal data kept in the solution-step buffer, one slot per time step - what the solver reads and writes. Non-historical data has one value and no history.

**LoRA.** Low-rank adapters for fine-tuning a large pretrained model by training small added matrices, merged back afterwards. `training.domino_finetune` uses it on DoMINO.

**`.mdlus`.** PhysicsNeMo's checkpoint format: a zip of the weights, the constructor arguments and the class name, so `Module.from_checkpoint` rebuilds the model without code.

**Model interface.** The setting on the point-cloud, graph and particle processes naming which upstream call convention the tensors are arranged for (`generic`, `transolver`, `figconvnet`, `meshgraphnet`, `bistride`, ...).

**Neural operator.** A model that learns a map between functions rather than fixed-size vectors - a whole field in, a whole field out, in principle at any resolution. FNO is the canonical one.

**OOD, OOD guard.** Out of distribution; a guard that flags inputs far from the training inputs so the surrogate is not silently extrapolated. `deployment.ood_guard_utils`.

**`.pmsh`.** The memory-mapped on-disk format of a PhysicsNeMo `Mesh`, the layout `MeshReader` consumes; the mesh exporter writes a series of them.

**POD, ROM.** Proper orthogonal decomposition and reduced-order models: a field represented as `u = Phi q` in a small basis. `bridges.rom_bridge` reads `RomApplication` bases.

**Process.** Kratos's unit of "something that happens at a fixed moment of the solve", attached from `ProjectParameters.json`. Every attachable thing here is one. See [Process reference](Process_Reference.html).

**Provenance map.** The record of which Kratos node or entity each simplex point and cell of a tessellation came from, so predictions on the tessellation can be written back exactly. See [Mesh Bridge](../Mesh_Bridge/Mesh_Bridge.html).

**Residual.** Three different things: the strong-form PDE residual (approximate, differentiable), the assembled solver residual (exact, a score) and the differentiable discrete residual (exact and differentiable). See [Symbolic and physics](../PhysicsNeMo_Basics/Symbolic_And_Physics.html).

**Row order.** The one array layout of the application: one row per entity, ids ascending, components minor. Everything produces and consumes it.

**ShardTensor, DTensor.** Torch's distributed tensor and PhysicsNeMo's subclass of it for domain parallelism: one tensor split across ranks, uneven shards allowed, halo exchange built into the operators.

**Simplex.** A point, segment, triangle or tetrahedron. PhysicsNeMo meshes contain nothing else, which is why the tessellation exists.

**Smallest-id diagonal.** The rule (Dompierre et al.) that splits every quadrilateral face along the diagonal through its smallest node id, so neighbouring elements triangulate a shared face identically and the tessellation is watertight.

**Sobolev training.** Training a surrogate on gradients as well as values - here, on Kratos's exact adjoint gradients. `training.sobolev_training`.

**Squeeze axis, thin axis.** The idiom for planar cases on 3-D grid models: sample the field on a grid with one axis of size 2, average that axis away before the model and duplicate it after.

**Surrogate.** A model that replaces an expensive computation with a cheap approximation of its output. Not a solver: it does not iterate or converge, and it is only as trustworthy as its training distribution.

**Tensor adaptor.** The core Kratos object giving a contiguous numpy view of entity data through a staging buffer (`CollectData`, `.data`, `StoreData`). The foundation of every bridge.

**Tessellation.** Splitting Kratos geometries (hexahedra, prisms, pyramids, quadrilaterals, quadratic elements) into simplices, with a provenance map back.

**Warm restart.** Re-initializing a trained model partially (shrink and perturb) before training on new data, to escape the loss of plasticity a fully trained network shows. `TrainModel`'s `"warm_restart"` block.

**Zero-copy.** A view that shares memory with its source. The tensor bridge is zero-copy with respect to the adaptor's staging buffer, not with respect to the Kratos entities themselves.
