---
title: Uncertainty and guardrails
keywords: physicsnemo uncertainty quantification concrete dropout ensemble gp head calibration ood guard geometry guardrail
tags: [Uncertainty_And_Guardrails.md]
sidebar: physicsnemo_application
summary: The uncertainty methods physicsnemo provides, the one variance decomposition they all feed, the metrics that decide whether an error bar is honest, and the guardrails that catch inputs a surrogate should not be trusted on.
---

# Uncertainty and guardrails

A surrogate returns a number. Nothing in that number says whether the input was anything like the training data. Upstream treats this as two separate problems and so does this application: **uncertainty** (attach a calibrated error bar to every prediction) and **guardrails** (refuse, or warn about, inputs the model was never trained for).

<p align="center">
    <img src="images/gp_head.svg" alt="Three routes to an error bar feeding one variance decomposition, and the three questions asked of the result"/>
</p>
<p align="center">Figure 1: Three ways to an error bar, one decomposition, three questions.</p>

## The methods

| Method | Where | Training cost | Inference cost | Gives | Distance-aware |
|---|---|---|---|---|---|
| **Concrete dropout** | `physicsnemo.nn.ConcreteDropout` | low - a regularizer term | 20 to 30 stochastic passes | epistemic only | no |
| **Ensembles** | no library needed | K trainings (deep), about one (snapshot, cyclic learning rate), free (checkpoint, last K epochs), K inferences (input ensemble over remeshed variants) | K passes | epistemic only | no |
| **`FieldVariationalGPHead`** | `physicsnemo.experimental.uq` (`uq-extras`, gpytorch) | moderate - a variational GP with deep kernel learning on the backbone features | one pass | mean, total variance **and** epistemic variance per point | yes |
| **`VariationalGPHead`** | same | same | one pass | the scalar version, for a drag coefficient rather than a field | yes |
| **Diffusion ensembles** | `physicsnemo.diffusion` | a generative model | M samples | a full conditional distribution | no |

Concrete dropout replaces a hand-tuned dropout rate with a learned one, which is what makes the resulting spread calibrated rather than arbitrary. Its documented trap: **a dropout layer in eval mode is a silent no-op** - forget to switch it back on for the sampling passes and you get zero uncertainty with no error.

The GP heads are closed form: one forward pass returns mean and variance, and the *epistemic* part grows with distance from the training features, which is the property the others lack. Upstream ships the class but not the training recipe, and without the recipe the variance collapses. The recipe is: seed the inducing points from real backbone features after a warm-up, anchor the posterior mean with an auxiliary MSE, ramp the KL term, set `n_train` to the number of training *points* (not geometries) for the field head, keep the features in float32, and use `feature_norm="l2_radial"` so the feature norm still carries distance. `deployment.uncertainty_utils` implements exactly that as `uncertainty_utils.FitGpHead`.

## One decomposition

Every method estimates the same thing:

<p align="center">$$\sigma^2_{\mathrm{total}}(x) \;=\; \sigma^2_{\mathrm{epistemic}}(x) \;+\; \sigma^2_{\mathrm{noise}}(x)$$</p>

The epistemic term is what more training data would remove; it is the signal that rises where the surrogate extrapolates. The noise term is the input-dependent scatter the data itself carries (aleatoric noise, or model discrepancy). Dropout and ensembles estimate the first only; GP heads give both from one pass; a mean-variance network gives only the total.

## Is the error bar honest?

An error bar that is too narrow is worse than none. The upstream UQ guide and this application both ask three questions, with the standardized residual as the common currency:

<p align="center">$$z_k = \frac{y_k - \mu_k}{\sigma_k}$$</p>

| Question | Metric | Target |
|---|---|---|
| Is the uncertainty the right size? (calibration) | z-RMS; coverage at 95 % (fraction within 1.96 sigma); negative log predictive density; sharpness (mean sigma, smaller is better *if* calibrated) | z-RMS near 1, coverage near 0.95 |
| Does high uncertainty mark high error? (discrimination) | rank correlation between sigma and the actual error; area under the sparsification-error curve (AUSE) | correlation high, AUSE near 0 |
| Does it grow off-distribution? | growth ratio (sigma_ood / sigma_id) divided by (rmse_ood / rmse_id) | near 1; below 1 means over-confident outside the training family |

Two measured facts from this application's own Examples case, because they are typical: a four-member ensemble's nominal 95 % bars covered **50 %** of the truth (textbook over-confidence, and exactly what the calibration metrics exist to catch); and below the training range all members were wrong *together*, so the spread stayed small while the error grew tenfold. The spread's direction was right; its magnitude was not guaranteed.

Calibration and error ranking peak at different checkpoints, so pick the checkpoint for the use you have; and normalizing features away removes the distance cue a GP head needs, which is what `"l2_radial"` preserves.

## Guardrails

Uncertainty answers "how sure is the model". A guardrail answers "should the model be asked at all".

- **`OODGuard`** (`physicsnemo.experimental.guardrails.embedded`) is calibrated on the training *inputs* - a kNN density over the model's input features - and checks each inference input against it. Upstream's `check()` only logs; this application captures that log and turns it into an advisory / strict / ignore policy per process. It is calibrated by `TrainModel` and saved as a sidecar next to the checkpoint.
- **`GeometryGuardrail`** (`physicsnemo.experimental.guardrails.geometry`, new in 2.2) works on the *shape*: it extracts non-invariant descriptors of a triangular surface mesh (translation, rotation and scale are deliberately kept, because a scaled car is a different car) and fits a Gaussian-mixture or polynomial-chaos density with warn and reject percentiles. Not bridged yet - it is on the roadmap, and the bridge's outward-oriented boundary surface is its natural input.

Guardrails and uncertainty are complementary, not substitutes: a guard catches the input that is far from everything seen; a calibrated variance tells you how much to trust the answer on the inputs the guard lets through.

## What this application uses it for

| PhysicsNeMo API | Kratos-side module | Gives you |
|---|---|---|
| `nn.ConcreteDropout` and the dropout utilities | `deployment.uncertainty_utils` | `"uncertainty" : {"method": "mc_dropout"}` on any inference process |
| checkpoint ensembles (no upstream API) | `deployment.uncertainty_utils` | `"method": "ensemble"`, and `"retain_ensemble"` to keep the members for CRPS |
| `experimental.uq.FieldVariationalGPHead` | `deployment.uncertainty_utils` (`uncertainty_utils.FitGpHead`, `uncertainty_utils.PredictWithGpHead`, saved as a `.gp_head.pt` sidecar since gpytorch modules cannot be TorchScripted) | `"method": "gp"` with calibrated epistemic variance |
| `metrics.general.calibration`, `crps`, `kcrps` | `deployment.uncertainty_utils`, `processes.validation_metrics_process` | the `"uncertainty_comparisons"` (coverage, NLL, sharpness, calibration error) and `"ensemble_comparisons"` (CRPS over named members) blocks |
| `experimental.guardrails.embedded.OODGuard` | `deployment.ood_guard_utils` | the `"ood_guard"` block, calibrated by `TrainModel` |
| diffusion ensembles | `processes.inference.diffusion_inference_process` | ensemble mean to the output fields, per-node std to `uncertainty_fields` |

Every uncertainty method writes its mean to the ordinary output fields and its standard deviation to dedicated uncertainty fields, so downstream exporters and outputs carry both without knowing which method produced them. A spread is de-normalized by scale only, never shifted - the one place the model-card machinery treats a field differently.

See [Uncertainty](../Uncertainty/Uncertainty.html) for the settings and the measured numbers.

Next: [Active learning](Active_Learning_Concepts.html).
