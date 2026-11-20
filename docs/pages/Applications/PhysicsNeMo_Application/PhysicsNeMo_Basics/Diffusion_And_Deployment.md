---
title: Diffusion and deployment
keywords: physicsnemo diffusion sampler preconditioner onnx metrics optim
tags: [Diffusion_And_Deployment.md]
sidebar: physicsnemo_application
summary: Generative models, the sampler/preconditioner split, metrics, optimizers and ONNX export.
---

# Diffusion and deployment

## `physicsnemo.diffusion`

A diffusion model does not predict *an* answer; it samples from a distribution of plausible answers. That is the point: run it several times and the spread between the samples is a calibrated statement about what the data does not determine.

The subpackage splits the problem into three replaceable parts, and knowing the split is most of what you need:

| Part | Answers | In physicsnemo |
|---|---|---|
| **Denoiser** | What network removes the noise? | any model implementing `diffusion.base.DiffusionModel` — the U-Nets, or `DiT` |
| **Preconditioner** | How is noise scaled and the network conditioned on the noise level? | `EDMPrecond`, `VEPrecond`, `VPPrecond`, `iDDPMPrecond`, `EDMPrecondSuperResolution` |
| **Sampler** | How do we walk from noise to a sample? | `deterministic_sampler`, `stochastic_sampler`, `EDMStochasticHeunSolver`, `EulerSolver`, `HeunSolver` |


<p align="center">
    <img src="images/diffusion_split.svg" alt="Denoiser, preconditioner and sampler as three swappable boxes, and the CorrDiff two-stage recipe below them"/>
</p>
<p align="center">Figure 1: The three replaceable parts, and CorrDiff's regression-plus-residual split.</p>

Also: `diffusion.guidance` (diffusion posterior sampling - `DataConsistencyDPSGuidance` steers samples toward masked observations, `ModelConsistencyDPSGuidance` toward a differentiable forward model, both at sampling time and without retraining), `diffusion.noise_schedulers` (EDM, VE, VP, iDDPM, one module per class since 2.2), `diffusion.multi_diffusion` (patch-based tiling of a 2-D domain larger than the training resolution), and `diffusion.metrics` with the EDM losses.

**Two APIs coexist in 2.2.** The *protocol* API - `DiffusionModel`, `Predictor` and `Denoiser` protocols, `EDMNoiseScheduler`, `EDMPreconditioner`, `MSEDSMLoss`, `samplers.sample(denoiser, latents, scheduler, solver=)` - is what guidance and multi-diffusion are written against. The *legacy* modules (`samplers.legacy_deterministic_sampler`, `metrics.legacy_losses`, `preconditioners.legacy`) each warn that they will be deprecated in a future release. This application's diffusion bridge runs on the protocol API (`"api" : "auto"`), keeping the legacy path reachable for the two CorrDiff losses, which have no protocol counterpart in 2.2, and for reproducing older runs. DPS guidance and multi-diffusion patching are bridged on top of it - see the [Diffusion](../Diffusion/Diffusion.html) page.

**CorrDiff** is the two-stage recipe that matters for physics: a *regression* model predicts the conditional mean, then a *residual* diffusion model learns what the regression could not. Predicting the mean with a deterministic model is much easier than making diffusion learn it, and the split shows.

## `physicsnemo.metrics`

- `metrics.general` — `crps` and `kcrps` (proper scoring rules for ensembles), `calibration`, `ensemble_metrics`, `entropy`, `histogram`, `mse`, `relative_error`, `power_spectrum`, `wasserstein`, `reduction`.
- `metrics.cae` — CFD-flavoured integrals and quantities.

CRPS is the one to know: it scores a whole *ensemble* against a single observation, rewarding both accuracy and honest spread. An over-confident ensemble scores badly even when its mean is right.

## `physicsnemo.optim`

`CombinedOptimizer` (different optimizers on different parameter groups) and `Muon`. Torch's own optimizers work fine; reach for these only when you need them.

## `physicsnemo.deploy.onnx`

`export_to_onnx_stream` and `run_onnx_inference`. ONNX is the portable artifact: a `.onnx` file plus ONNX Runtime needs **neither** physicsnemo **nor** torch to run, which is what makes it the right thing to hand to a production host.

Caveat: not every operator exports. The FFTs inside FNO-style models are not supported by the CPU ONNX Runtime; MLP and convolutional models export and run everywhere.

## What this application uses it for

| PhysicsNeMo API | Kratos-side module | Notes |
|---|---|---|
| `diffusion.preconditioners.EDMPrecondSuperResolution` | `training.diffusion_utils` | the conditional EDM training loss |
| `diffusion.samplers` | `processes.inference.diffusion_inference_process` | ensemble sampling, ensemble mean written to the output fields, per-node std to uncertainty fields |
| `models.diffusion_unets.CorrDiffRegressionUNet` | `training.diffusion_utils` | `TrainCorrDiffPair`, the two-stage recipe |
| `models.dit.DiT` | `training.diffusion_utils` | `denoiser_interface: "dit"` — conditioning concatenated into the input channels, sigma broadcast to timesteps |
| `experimental.models.diffusion_unets.DiffusionUNet3D` | `training.diffusion_utils` | `denoiser_interface: "unet3d"` — volumetric 5-D grids, conditioning passed natively as the `TensorDict` `"volume"` key; `MSEDSMLoss` is rank-agnostic, so no rank-generalized clone is needed |
| `models.topodiff.TopoDiff` | `training.diffusion_utils` | `denoiser_interface: "topodiff"` — constraint channels passed natively as the model's second argument rather than concatenated |
| `diffusion.guidance` (`DataConsistencyDPSGuidance`, `ModelConsistencyDPSGuidance`) | `training.diffusion_utils`, `physics.diffusion_residual_operator` | the `"guidance"` block — diffusion posterior sampling on an already-trained denoiser, the shipped operator being the exact discrete FEM residual |
| `diffusion.multi_diffusion` | `training.diffusion_utils` | the `"patching"` block — a model trained at patch resolution samples a larger grid; 2-D only upstream |
| `metrics.general.crps` | `processes.validation_metrics_process`, `deployment.uncertainty_utils` | ensemble scoring, calibration |
| `deploy.onnx.export_to_onnx_stream` | `training.training_utils.ExportOnnxModel` | writes the `.onnx` plus its model card |
| ONNX Runtime | `deployment.onnx_utils`, `processes.inference.onnx_inference_process` | cached session, GPU support |

**Two ONNX Runtime failures reported as success**, both guarded here: a silent fall back to CPU (a missing CUDA build, or a nonexistent device index, both yield a working CPU session — `require_device` turns that into an error), and a dropped device index (`"cuda:1"` running on device 0).

Serving beyond a single process is `deployment.triton_export` plus `processes.inference.triton_inference_process`: the solver host then needs neither the weights nor a GPU.

See [Diffusion](../Diffusion/Diffusion.html) and [Inference](../Inference/Inference.html).

Next: [Distributed and scale](Distributed_And_Scale.html).
