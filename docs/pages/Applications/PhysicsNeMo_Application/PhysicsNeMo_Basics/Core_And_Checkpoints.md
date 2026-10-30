---
title: Core and checkpoints
keywords: physicsnemo module checkpoint mdlus model card
tags: [Core_And_Checkpoints.md]
sidebar: physicsnemo_application
summary: physicsnemo.Module, the .mdlus format, and the model card this application layers on top.
---

# Core and checkpoints

`physicsnemo.core` is small and you will use two things from it constantly: `Module` and the checkpoint format it defines.

## `physicsnemo.Module`

A `physicsnemo.Module` **is** a `torch.nn.Module` — it subclasses it, so everything you know about `forward`, `parameters()`, `.to(device)` and `state_dict()` applies unchanged. What it adds is that it **records the arguments it was constructed with**.

That single addition is the reason the checkpoint format works. A plain torch `state_dict` is a bag of tensors: to load it you must first reconstruct the architecture yourself, in code, with exactly the hyperparameters it was trained with. A `physicsnemo.Module` writes those hyperparameters into the checkpoint, so loading needs no architecture code at all:

```python
import physicsnemo
model = physicsnemo.Module.from_checkpoint("surrogate.mdlus")   # that is all
```

`ModelMetaData` is the descriptor a model class carries (name, whether it supports CUDA graphs, AMP, ONNX export, and so on). You only touch it when writing a new architecture.

## What is inside a `.mdlus` file

A zip archive containing:

- the `state_dict` — the weights;
- the constructor arguments — the architecture;
- the class's registry name, so `from_checkpoint` can find the class again.

<p align="center">
    <img src="images/mdlus_and_card.svg" alt="A .mdlus archive, a TorchScript file and the model-card sidecar side by side, and the gather, normalize, forward, de-normalize, scatter sequence a deployment process runs with them"/>
</p>
<p align="center">Figure 1: The two checkpoint formats, the card that gives their tensors meaning, and what a deployment process does with all three every step.</p>

The consequences matter in practice:

- **A `.mdlus` needs physicsnemo installed to load.** The class it names has to exist. If you want an artifact that runs without physicsnemo, export to ONNX (see [Diffusion and deployment](Diffusion_And_Deployment.html)).
- **A `.mdlus` does not need your training script.** This is the point.
- **It does not record what the numbers mean.** Nothing in the format says which field is in channel 0, what units it is in, or whether the targets were normalized. That gap is what the model card below fills.

`physicsnemo.core.registry` is the class registry `from_checkpoint` looks names up in; `physicsnemo.compat` handles loading checkpoints written by older releases.

## TorchScript, the other format

A TorchScript file (`.pt` written by `torch.jit.script`/`trace`) is the alternative: it stores the *traced computation*, not the architecture, so it loads with only torch installed and no physicsnemo at all. In exchange it is frozen — no `torch.compile`, no easy surgery, and anything not traceable (a gpytorch head, for instance) simply cannot go in one.

Both formats are supported everywhere in this application. Rule of thumb: **`.mdlus` while you are still training, TorchScript or ONNX when you ship.**

## The model card — this application's addition

A checkpoint that says nothing about its fields is a checkpoint you can silently misuse. Writing a model's raw output onto `TEMPERATURE` when it was trained on targets scaled to zero mean and unit variance produces finite, plausible-looking, completely wrong numbers.

So `deployment.model_registry` writes a **model card**: a JSON sidecar next to the checkpoint naming the input and output fields, their order, and — the part that prevents the failure above — `"output_normalization"`, the scaling the targets were trained under, and its mirror `"input_normalization"`, the scaling the features were trained under. Every deployment process validates its configuration against the card, standardizes what it feeds the model, and inverts the output normalization before touching a Kratos variable. An absent key is exactly the identity, so cards are optional; a wrong channel count raises rather than producing silent NaNs.

This is the single most useful thing to know about the application: **if you train a model here, save it with a card.**

## What this application uses it for

| PhysicsNeMo API | Kratos-side module | Used by |
|---|---|---|
| `physicsnemo.Module.from_checkpoint` | `deployment.model_registry` | every deployment process |
| `ModelMetaData` | — | only when writing a new architecture |
| `physicsnemo.core.registry` | `deployment.model_registry` | resolving a checkpoint's class |

Kratos-side entry points: `model_registry.LoadModel`, and `SaveTrainedModel` in `training.training_utils`, which writes `.mdlus` or TorchScript and gathers an FSDP2 model's shards to full tensors first.

Next: [Models](Models.html) — what to put in the checkpoint.
