---
title: From scratch
keywords: tutorial getting started surrogate workflow
tags: [From_Scratch.md]
sidebar: physicsnemo_application
summary: One surrogate, end to end - export from a real solve, train, save with a card, deploy in the loop, validate.
---

# From scratch

This walks the whole path once, naming the exact module at every step. It is deliberately the *simplest* version of each: a plain MLP, nodal fields, no uncertainty, no distribution. Everything else in this documentation is a variation on these five steps.

If PhysicsNeMo itself is new to you, read [PhysicsNeMo Basics](../PhysicsNeMo_Basics/Overview.html) first. If you want to know where a module lives, [Where things live](Module_Map.html).

**Prerequisites:** a Kratos build with this application and a solver application (the examples below use `ConvectionDiffusionApplication`), plus `pip install torch nvidia-physicsnemo`.

## The shape of it

<p align="center">
    <img src="images/lifecycle.svg" alt="Export, train, save, deploy, validate, with the Kratos, application and PhysicsNeMo pieces under each step"/>
</p>
<p align="center">Figure 1: The five steps and who owns what at each.</p>

```
 1. export     a real solve  ---> sample_0.npz, sample_1.npz, ...
 2. train      those files   ---> a torch model
 3. save       the model     ---> surrogate.pt + surrogate.pt.card.json
 4. deploy     the model     ---> predictions on Kratos nodes, in the loop
 5. validate   the loop      ---> numbers saying whether to believe it
```

## 1. Export training data from a real solve

Attach `dataset_export_process` to a solve you already have. It writes one `.npz` per step with the fields you name, and imports no torch — you can run this step on a machine with no ML stack at all.

```json
{
    "python_module" : "dataset_export_process",
    "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.export",
    "Parameters"    : {
        "model_part_name" : "ThermalModelPart",
        "list_of_fields"  : [
            { "variable_name" : "HEAT_FLUX",   "data_location" : "node_historical" },
            { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" }
        ],
        "output_path"     : "dataset",
        "output_interval" : 1
    }
}
```

Run the solve — ideally several, sweeping whatever you want the surrogate to generalize over (a boundary condition, a material parameter, a load). One solve gives you one sample; a surrogate needs a family.

The field keys in the resulting files are `"<VARIABLE>__<location>"`, e.g. `"TEMPERATURE__node_historical"`. You will need them in the next step.

*Variations:* `grid_dataset_export_process` resamples onto a voxel grid; `mesh_export_process` writes a `.pmsh` series; `streaming_dataset_export_process` skips the files entirely and trains out of the running solve.

## 2. Train

```python
import torch
import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.training import torch_dataset, training_utils

dataset = torch_dataset.CreateNpzDataset(
    "dataset",
    input_keys=["HEAT_FLUX__node_historical"],
    output_keys=["TEMPERATURE__node_historical"])

model = torch.nn.Sequential(
    torch.nn.Linear(1, 64), torch.nn.Tanh(),
    torch.nn.Linear(64, 64), torch.nn.Tanh(),
    torch.nn.Linear(64, 1))

history = training_utils.TrainModel(model, dataset, Kratos.Parameters("""{
    "epochs"        : 500,
    "batch_size"    : 8,
    "learning_rate" : 1e-3,
    "echo_interval" : 50
}"""))
```

`CreateNpzDataset` yields `(inputs, targets)` as `(n_entities, width)` float32 tensors — **the same layout `inference_process` will feed the model at deployment**. That correspondence is the whole reason the pieces fit together; if you build a dataset by hand, match it.

A `torch.nn.Module` is fine here. Use a `physicsnemo.Module` (see [Models](../PhysicsNeMo_Basics/Models.html)) when you want the architecture to travel in the checkpoint.

*Variations:* `extra_loss_terms=` adds a physics residual to the objective ([Symbolic and physics](../PhysicsNeMo_Basics/Symbolic_And_Physics.html)); `epoch_callbacks=` lets you score the surrogate against the real PDE residual while it trains.

## 3. Save it — with a card

```python
card = {
    "input_fields"  : [{"variable_name": "HEAT_FLUX",   "data_location": "node_historical"}],
    "output_fields" : [{"variable_name": "TEMPERATURE", "data_location": "node_historical"}],
    "notes"         : "swept HEAT_FLUX over [0, 1e4], 40 solves",
}
checkpoint_type = training_utils.SaveTrainedModel(model, "surrogate.pt", card=card)
```

**Do not skip the card.** It writes `surrogate.pt.card.json` next to the checkpoint, recording what the channels mean, and every deployment process validates its configuration against it. If you trained on normalized targets, add `"output_normalization"` to the card — the deployment path then inverts the scaling for you, instead of writing normalized numbers onto a physical variable where they look plausible and are wrong by the training scaling. If you standardized the inputs too, add `"input_normalization"` the same way (`model_registry.MakeMeanStdNormalization(mean, std)` builds either entry) and the deployment path standardizes what it feeds the model.

## 4. Deploy it inside a solve

Attach `inference_process` to the analysis you want to accelerate. It gathers the input fields, runs a no-grad forward pass, and writes the output back onto real Kratos variables — after which nothing downstream can tell a predicted value from a solved one.

```json
{
    "python_module" : "inference_process",
    "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes.inference",
    "Parameters"    : {
        "model_part_name" : "ThermalModelPart",
        "model_settings"  : { "checkpoint_file" : "surrogate.pt", "device" : "cpu" },
        "input_fields"    : [ { "variable_name" : "HEAT_FLUX",   "data_location" : "node_historical" } ],
        "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
        "execution_point" : "finalize_solution_step"
    }
}
```

*Variations:* `hybrid_initialization_process` warm-starts the solver from the prediction instead of replacing it — often the better trade, since you keep the solver's convergence guarantee and only spend fewer iterations getting there. For meshes, grids, point clouds or time series, swap in the matching process from [Where things live](Module_Map.html).

## 5. Decide whether to believe it

A surrogate that runs is not a surrogate that works.

```json
{
    "python_module" : "validation_metrics_process",
    "kratos_module" : "KratosMultiphysics.PhysicsNeMoApplication.processes",
    "Parameters"    : {
        "model_part_name"   : "ThermalModelPart",
        "list_of_comparisons" : [ {
            "predicted_variable" : "TEMPERATURE",
            "predicted_location" : "node_historical",
            "reference_variable" : "TEMPERATURE_REFERENCE",
            "reference_location" : "node_historical",
            "metrics"            : [ "mse", "rmse", "max_abs_error" ]
        } ],
        "output_file"       : "validation.json"
    }
}
```

Then, when the answer matters:

- **error bars** — the `"uncertainty"` block on any deployment process (MC dropout, a checkpoint ensemble, or a GP head);
- **a guard** — `"ood_guard"`, calibrated on the training inputs by `TrainModel` and checked per inference, so an input far outside the training distribution is flagged rather than silently extrapolated;
- **calibration metrics** — whether those error bars are honest at all.

See [Uncertainty](../Uncertainty/Uncertainty.html).

## Where to go next

| If you want | Go to |
|---|---|
| The same path, executable | `examples/notebooks/03_training_a_surrogate.ipynb` |
| The same path on a real solver | notebook 07, and the [thermal surrogate lifecycle](https://github.com/KratosMultiphysics/Examples/tree/master/physics_nemo_application/use_cases/thermal_surrogate_lifecycle) use case |
| Unstructured meshes rather than nodal vectors | [Mesh Bridge](../Mesh_Bridge/Mesh_Bridge.html), [GNNs](../Graph_Neural_Networks/Graph_Neural_Networks.html) |
| Physics in the loss | [Physics-Informed](../Physics_Informed/Physics_Informed.html) |
| The model to pick its own training data | [Active Learning](../Active_Learning/Active_Learning.html) |
| To run across ranks | [Distributed](../Distributed/Distributed.html) |
