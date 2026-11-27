---
title: Diffusion
keywords: diffusion corrdiff downscaling generative uncertainty
tags: [Diffusion.md]
sidebar: physicsnemo_application
summary: Conditional diffusion on grids - training with the EDM loss, ensemble deployment with uncertainty fields, DiT and volumetric U-Net denoisers, the CorrDiff two-stage recipe and FWI-style inversion.
---

# Conditional diffusion field models

CorrDiff-style downscaling with `physicsnemo.diffusion`: a denoiser learns the distribution of fine fields **conditioned** on a coarse (or otherwise partial) field, and sampling it repeatedly yields an ensemble whose mean is the prediction and whose spread is a calibrated uncertainty band.

<p align="center">
    <img src="../PhysicsNeMo_Basics/images/diffusion_split.svg" alt="Denoiser, preconditioner and sampler, and the CorrDiff two-stage recipe"/>
</p>
<p align="center">Figure 1: The three replaceable parts, and the regression-plus-residual split this page's recipe implements.</p>

## Training

`diffusion_utils.TrainDiffusionModel(model, dataset, settings)` runs the conditional EDM loss (`EDMLossSR`) over `(condition, target)` grid pairs — `CreateGridPairDataset(input_directory, target_directory)` builds them from two `GridDatasetExportProcess` outputs matched by step:

```python
from physicsnemo.diffusion.preconditioners import EDMPrecondSuperResolution
from KratosMultiphysics.PhysicsNeMoApplication.training import diffusion_utils, training_utils
from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import CreateGridPairDataset

model = EDMPrecondSuperResolution(
    img_resolution=64, img_in_channels=1, img_out_channels=1,
    model_type="SongUNet")                       # a physicsnemo Module
dataset = CreateGridPairDataset("coarse_grids", "fine_grids", squeeze_axis=2)
diffusion_utils.TrainDiffusionModel(model, dataset, Kratos.Parameters("""{ "epochs": 500 }"""))
training_utils.SaveTrainedModel(model, "downscaler.mdlus")   # regular .mdlus checkpoint
```

Settings additionally expose the EDM noise schedule (`P_mean`, `P_std`, `sigma_data`, `sigma_min`, `sigma_max`).

**Which physicsnemo API the loss runs on.** physicsnemo 2.2 replaced the diffusion stack this bridge was first written against: every class under `metrics.legacy_losses`, `samplers.legacy_deterministic_sampler` and `preconditioners.legacy` warns that it will be deprecated. `TrainDiffusionModel` now runs the conditional EDM loss on the **protocol API** - `MSEDSMLoss` over an `EDMNoiseScheduler` - and `GenerateEnsemble` on `samplers.sample`. The `"api"` setting chooses:

| `"api"` | Effect |
|---|---|
| `"auto"` (default) | protocol for `"edm_sr"`, legacy for the two CorrDiff losses |
| `"protocol"` | protocol only; a CorrDiff loss is refused with an explanation rather than silently downgraded |
| `"legacy"` | the deprecated modules, for reproducing an older run |

The two CorrDiff stages stay on the legacy losses because physicsnemo 2.2 ships **no protocol equivalent** of `RegressionLoss`/`ResidualLoss` (`diffusion.metrics.losses` has only `MSEDSMLoss` and `WeightedMSEDSMLoss`). Numbers are not bit-identical across the two paths: the legacy sampler ran its loop in float64, the protocol one runs at the latent's dtype.

**The denoiser interface is always explicit**, never inferred. `"denoiser_interface"` is one of `"edm"` (the `net(x, img_lr, sigma)` super-resolution contract - the default), `"dit"`, `"unet3d"`, or `"protocol"` for a model that already speaks `model(x, t, condition=...)`. Sniffing the signature would be wrong: both `DiT` and `DiffusionUNet3D` declare a `condition` parameter, but DiT's is a `(B, condition_dim)` label **vector** (a conditioning grid belongs in its input channels) and `DiffusionUNet3D`'s is a `TensorDict`, so a probe would feed each the wrong thing - in DiT's case without any shape error.

`"preconditioner"` (`"none"` by default, or `"edm"`) wraps a raw x0-predictor in `EDMPreconditioner`. It is **refused for `"edm"`**, whose models are preconditioned already; wrapping one again preconditions it twice and produces finite garbage. Train and deploy must agree on this choice.

## Deployment

`DiffusionInferenceProcess` samples the condition fields onto a grid, generates a reverse-diffusion ensemble (`GenerateEnsemble`: `num_samples`, `num_steps`, `solver`), and scatters the ensemble **mean** onto `output_fields` — and, when configured, the ensemble **standard deviation** onto `uncertainty_fields`:

```json
{
    "python_module" : "diffusion_inference_process",
    "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.inference",
    "Parameters"    : {
        "model_part_name"        : "CoarsePart",
        "output_model_part_name" : "FinePart",
        "model_settings"         : { "checkpoint_file" : "downscaler.mdlus", "checkpoint_type" : "physicsnemo" },
        "input_fields"           : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
        "output_fields"          : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_non_historical" } ],
        "uncertainty_fields"     : [ { "variable_name" : "NODAL_ERROR", "data_location" : "node_non_historical" } ],
        "grid_shape"             : [64, 64, 2],
        "squeeze_axis"           : 2,
        "sampler_settings"       : { "num_samples" : 16, "num_steps" : 18 }
    }
}
```

SongUNet-based denoisers are 2D: planar Kratos cases use the thin-axis idiom (`squeeze_axis`, see the Sequence Models page). Genuinely volumetric grids run through `"denoiser_interface": "unet3d"` instead — see below.

On the protocol API `"solver"` is one of `"heun"` (default), `"euler"`, `"edm_stochastic_euler"` and `"edm_stochastic_heun"`, and `"solver_options"` reaches the solver's own constructor - `S_churn` and friends live there, and the deterministic solvers reject them. Two traps the bridge turns into errors: `"num_steps"` must be at least 2 (the schedule spans `num_steps - 1` intervals, so a single step divides by zero and yields NaN timesteps), and the stochastic solvers scale their churn by their **own** `num_steps`, which defaults to 18 whatever the sampler was asked for - the bridge injects the real value.

**Model cards.** The denoiser's card (`model_settings`) carries the `"output_normalization"` of what the ensemble emits — in the two-stage recipe, regression mean included. The ensemble mean is scaled and shifted, the spread written to the `uncertainty_fields` is **scaled only**, both per channel along axis 0 (see the [Inference](../Inference/Inference.html) page).


<p align="center">
    <img src="images/diffusion_mesh.png" alt="Blurred condition, diffusion ensemble mean, sharp truth and ensemble std on the mesh"/>
</p>
<p align="center">Figure 2: Notebook 09 - what DiffusionInferenceProcess read and wrote, on the mesh: condition, ensemble mean, truth and the per-node standard deviation.</p>

## Diffusion posterior sampling (DPS) guidance

A trained denoiser encodes a prior over fields. DPS adds a **measurement** term at sampling time, so the same checkpoint can be steered toward data it never saw during training - no retraining, no second model. At every sampler step the current estimate of the clean field is scored by a differentiable term, and the gradient of that score is added to the diffusion model's own score.

The `"guidance"` block inside `sampler_settings` turns it on:

```json
"sampler_settings" : {
    "num_samples" : 16,
    "num_steps"   : 18,
    "guidance"    : {
        "type"               : "data_consistency",
        "std_y"              : 0.5,
        "observation_fields" : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_non_historical" } ],
        "mask_field"         : { "variable_name" : "NODAL_AREA", "data_location" : "node_non_historical" }
    }
}
```

| Key | Default | Meaning |
|---|---|---|
| `type` | `"none"` | `"data_consistency"` (measurements at known entries) or `"model_consistency"` (a differentiable forward operator) |
| `std_y` | `0.1` | the observation noise level, and with it the guidance strength |
| `gamma` | `0.0` | the time-dependent term; a non-zero value also needs the schedule's sigma, which the bridge supplies |
| `norm` | `2` | the exponent of the measurement misfit |
| `observation_fields` | `[]` | the measured field, read from the model part like any other |
| `mask_field` | - | a 0/1 nodal field marking which entries were measured (`data_consistency`) |
| `operator` | `"kratos_residual"` | which forward operator to use (`model_consistency`) |
| `residual_fields` | `[]` | the DOF fields the generated channels correspond to |
| `residual_model_part_name` | `""` | the part whose system is assembled, when it is not the process's own |

**`std_y` is the knob, and it bites.** The guidance step scales as `1/(2 std_y^2)`, so a value far below the observations' own scale overshoots and the sampler diverges - upstream reports that as ordinary finite numbers growing to infinity, not as an error. The bridge checks the finished sample and raises, naming `std_y`, rather than scattering `inf` onto the mesh.

Observations are read in physical units and mapped into the model's own **normalized** output space before sampling, because that is the space the sampler works in; the de-normalization happens afterwards on the finished ensemble.

### The exact discrete residual as the operator

`"type" : "model_consistency"` with `"operator" : "kratos_residual"` makes the physics itself the measurement: `physics/diffusion_residual_operator.py` maps a generated grid to the assembled residual `b(u)` of the real Kratos system, and the observation is zero - a field that solves the PDE has no residual. It is built on the same autograd Function as the [exact residual loss](../Physics_Informed/Physics_Informed.html) (forward = `BuildRHS`, backward = the consistent tangent's transpose), so the gradient is the discretely exact one.

The chain it runs is the inference process's own, in reverse order, so a guided sample is scored in the units it will be written in: de-normalize with the card, un-squeeze the thin axis, interpolate at the nodes with the **differentiable** `grid_bridge.InterpolateGridAtPointsTorch`, gather DOFs, overwrite the fixed DOFs with the model part's own Dirichlet values, assemble.

One consequence worth knowing: assembling a residual is how Kratos evaluates one, so it **writes the trial field into the model part's solution-step database**. The operator saves the incoming state and the process restores it in a `finally`, so a guided inference leaves the solve exactly as it found it.

This is the retraining-free alternative to the mask-conditioned model of the [subsurface inversion recipe](#subsurface-inversion-fwi-style) below: that one learns to honour sensors, this one enforces them at sampling time on any denoiser.

## Multi-diffusion patching

`physicsnemo.diffusion.multi_diffusion` tiles a latent larger than the resolution a model was trained at. Training draws random patches, sampling fuses a grid of overlapping ones:

```python
diffusion_utils.TrainDiffusionModel(model, dataset, Kratos.Parameters("""{
    "epochs"   : 500,
    "patching" : { "patch_shape" : [64, 64], "patch_num" : 4 }
}"""))
```

```json
"sampler_settings" : {
    "patching" : { "patch_shape" : [64, 64], "overlap_pix" : 8, "boundary_pix" : 0, "chunk_size" : 0 }
}
```

`patch_shape` must match the backbone's trained resolution, and upstream's own geometry rules apply (`patch_shape` no larger than the grid, `patch - overlap - boundary` at least 1). `chunk_size` bounds how many patches are evaluated at once.

Two limits are upstream's, not the bridge's: multi-diffusion is **2-D only**, so the volumetric `"unet3d"` path stays whole-grid and the bridge says so rather than failing deep inside; and the patch-local DPS guidance variants score each patch independently, so a globally coupled operator like the Kratos residual runs through the ordinary guidance over the fused predictor instead. A `MultiDiffusionPredictor` also **mutates** the wrapper it is given, so the bridge builds a fresh one per ensemble rather than reusing a training wrapper.

## Variations

The process is agnostic to what the condition channels mean:

- **Downscaling (CorrDiff)**: condition = the coarse solve, target = the fine field.
- **Generative design (TopoDiff-style)**: condition = constraint/mask fields, target = the design field.
- **Flow reconstruction from sparse data**: condition = masked observations (zeros where unobserved plus an indicator channel), target = the full field.

Only the training pairs differ — build them with `CreateGridPairDataset` over appropriately exported series.

The ensemble-mean + per-node uncertainty pattern this process introduced is now available for **every** deployed model through the `"uncertainty"` block (MC dropout or checkpoint ensembles) — see [Uncertainty and Governance](../Uncertainty/Uncertainty.html).

## DiT denoisers

`physicsnemo.models.dit.DiT` (diffusion transformer) plugs into the same bridge through `"denoiser_interface": "dit"` on `DiffusionInferenceProcess` (or `diffusion_utils.WrapDenoiser(dit, "dit")` in training scripts): the wrapper maps the EDM contract `net(x, img_lr, sigma)` onto `dit(x, t)` by concatenating the conditioning grid into the input channels — construct the DiT with `in_channels = C_out + C_cond` and `out_channels = C_out` — broadcasting `sigma` to the per-sample timestep tensor and casting at the float64-sampler/float32-weights boundary. RoPE, invalid-region masking and attention backends are DiT **construction** choices (`block_kwargs`/`attn_kwargs`/`attention_backend`); the wrapper only standardizes the forward call, and the raw DiT acts as the denoiser directly (no EDM pre/post-scaling), so train it through the same wrapper.

Note on upstream naming: `physicsnemo.models.diffusion` was renamed to `physicsnemo.models.diffusion_unets` (`SongUNet`, `UNet`, `DhariwalUNet`, ...). There is **no** compatibility shim: `import physicsnemo.models.diffusion` raises `ModuleNotFoundError` on 2.2.0, so the old path must be updated rather than relied on.

## Volumetric (3D) U-Net denoisers

`physicsnemo.experimental.models.diffusion_unets.DiffusionUNet3D` — a genuine volumetric diffusion U-Net, under `experimental` in 2.2 — plugs in through `"denoiser_interface": "unet3d"` (or `diffusion_utils.WrapDenoiser(unet, "unet3d")` in training scripts). Three things differ from `"dit"`:

- **Conditioning is native, not concatenated**: the wrapper passes the conditioning grid as the model's `TensorDict` `condition["volume"]` — construct the model with `x_channels = C_out` (its output width equals its latent width) and `vol_cond_channels = C_cond`. `sigma` broadcasts to the per-sample timestep tensor exactly as for DiT.
- **The grid stays 5-D**: the condition is the full `(C, D, H, W)` sample — `squeeze_axis` is rejected in this mode, and each spatial extent must be a power of 2 or a multiple of `2**(num_levels - 1)` (the model validates this itself).
- **Training is rank-agnostic on the protocol API**: the legacy `EDMLossSR` hard-codes the 4-D image rank (`randn([B, 1, 1, 1])`), so volumetric samples could not broadcast against it and the bridge used to carry a local clone of the loss whose noise draw followed the batch rank. `MSEDSMLoss` needs no such thing, so the clone is gone; `"api" : "legacy"` on a volumetric dataset now says why it cannot work. The 2D-only CorrDiff losses (`"regression"`/`"residual"`) still reject volumetric input with a clear error.

`GenerateEnsemble`, the deterministic sampler and the ensemble-mean/uncertainty scatter are rank-agnostic and work unchanged.

## CorrDiff two-stage recipe (regression + residual diffusion)

CorrDiff/StormCast-style downscaling splits the prediction: a deterministic **regression** stage learns the conditional mean, and the diffusion stage denoises only the **residual** `target − regression(condition)` — sharper ensembles, better-calibrated spread. The bridge ships the full recipe on verified 2.2.0 APIs:

- **Stage 1**: `TrainDiffusionModel(regression_model, dataset, settings)` with `"loss": "regression"` — `physicsnemo.diffusion`'s `RegressionLoss` (plain MSE against `net(zeros, img_lr)`) on a `physicsnemo.models.diffusion_unets.CorrDiffRegressionUNet`.
- **Stage 2**: `"loss": "residual"` with `regression_model=` the trained stage 1 — `ResidualLoss(regression_net=...)` drives the same `net(x, img_lr, sigma)` denoiser interface. The regression net is frozen explicitly (`.eval()` + `requires_grad_(False)` — upstream does **not** `no_grad` it), and `P_mean` defaults to upstream's `0.0` for this loss (vs `−1.2` for `edm_sr`; explicit values always win). **The residual-stage denoiser must wrap `SongUNetPosEmbd`** (`model_type="SongUNetPosEmbd"`; its `N_grid_channels`, default 4, count toward `img_in_channels` — upstream CorrDiff's own sizing): `ResidualLoss` passes positional-embedding kwargs a plain `SongUNet` rejects.
- `TrainCorrDiffPair(regression_model, diffusion_model, dataset, settings)` runs both stages (`"regression_epochs"`/`"regression_learning_rate"` override the shared values).
- **Inference**: `DiffusionInferenceProcess` gains an optional `"regression_settings"` block (`model_settings`-shaped, loaded through the model registry with card checks); when present, `RunRegressionMean`'s prediction is added to the generated ensemble **before** the mean/std are taken — the mean shifts by the regression stage, the ensemble spread stays the denoiser's, exactly CorrDiff inference. Pinned by a test asserting the with-regression run differs from the plain run (same sampler seed) by exactly the scattered regression mean.

Static/invariant conditioning — topography-like fields for `WindEngineeringApplication`/`DamApplication` cases — is just an extra `input_fields` entry (a nodal field that never changes).


<p align="center">
    <img src="images/corrdiff_mesh.png" alt="Blurred condition, CorrDiff two-stage mean, sharp truth and residual ensemble std on the mesh"/>
</p>
<p align="center">Figure 3: Notebook 15 - the two-stage recipe on the same case: the regression mean plus the residual ensemble sharpens the condition; the residual spread is the uncertainty.</p>

## TopoDiff and generative design on real compliance data

`physicsnemo.models.topodiff.TopoDiff` is built for topology optimization and takes its constraint channels NATIVELY, as a second argument - `forward(x, cons, timesteps)` - where DiT wants them concatenated into the input. It reaches the bridge as `"denoiser_interface" : "topodiff"`, or `WrapDiffusionModel(model, "topodiff", out_channels=1)` in a training script.

The data is a real solve, not a picture. `tests/kratos_solver_cases/compliance_case.py` is a SIMP cantilever: each element carries its own Young's modulus through `E = E_min + rho^p (E_0 - E_min)`, so a DENSITY FIELD is the design variable and the structure's compliance is what a design is judged by. `ConstraintChannels` supplies the conditioning a generative model needs - where the structure is held, where it is loaded, and how much material it may use - and a sampled design can be handed straight back to Kratos to be scored.

Three traps, all of them pinned:

- **TopoDiff exposes no `out_channels`**, so the bridge cannot infer the channel count and says so instead of guessing.
- **`model_channels` below 64 fails at construction** with a message about `num_heads`, not about width: the attention layer derives its head count from the width and ends up with zero.
- **A nodal `POINT_LOAD` value alone does nothing.** Without a point-load CONDITION the solve returns zero displacement and a compliance of exactly zero, which reads like a converged answer rather than a missing load.

## Subsurface inversion (FWI-style)

Inversion by conditional diffusion needs **nothing beyond the shipped bridge**: condition = sparse observations (e.g. borehole columns) plus a binary observation-mask channel, target = the subsurface property grid; `GenerateEnsemble`'s per-node standard deviation is the inversion uncertainty. The layered-earth recipe is pinned by `tests/test_corrdiff_recipe.py::TestFwiInversionRecipe`; the `GeoMechanicsApplication` pairing is availability-gated (not compiled in the reference environment).

The two ways to honour sensors are now both available, and they trade off differently. This recipe **trains** a mask-conditioned model, which learns what fields consistent with sparse readings look like and costs nothing extra at sampling time. [DPS guidance](#diffusion-posterior-sampling-dps-guidance) enforces the readings at **sampling** time on any already-trained denoiser, including one that never saw a mask - no retraining, but every sampler step pays for a gradient through the model.
