---
title: Active learning
keywords: physicsnemo active learning driver query strategy label strategy metrology protocols
tags: [Active_Learning_Concepts.md]
sidebar: physicsnemo_application
summary: The physicsnemo.active_learning loop - its four phases, the protocols you implement, the Driver that runs them, and the seat Kratos takes in it.
---

# Active learning

Active learning inverts the usual order. Instead of generating a dataset and then training, the model trains, decides where it is least certain, asks for labels *there*, trains again. When a label costs a finite-element solve, choosing which solves to run is the whole game - and the solver is Kratos.

`physicsnemo.active_learning` is the bundled framework for that loop. This page explains it on its own terms; [Active Learning](../Active_Learning/Active_Learning.html) documents the Kratos-side implementation.

## The loop

<div class="mermaid">
flowchart LR
    pool[Unlabeled pool of candidate cases]
    subgraph driver [Driver, one active_learning_step]
        direction LR
        train[1. Train or fine-tune the surrogate]
        metro[2. Metrology, measure progress]
        query[3. Query strategy picks the next cases]
        label[4. Label strategy produces ground truth]
        train --> metro --> query --> label
    end
    pool --> query
    query -->|query queue| label
    label -->|serialize queue, labeled samples| train
    train -.-> ckpt[.mdlus checkpoint and driver_log.json]
    subgraph kratos [Kratos, the labeler]
        direction TB
        backend[execution backend, in-process or subprocess]
        solve[one AnalysisStage per case, with dataset_export_process]
        backend --> solve
    end
    label --- kratos
</div>

One `active_learning_step()` runs the four phases in order; `Driver.run()` repeats it until `max_steps`. Any phase can be switched off (`skip_training`, `skip_metrology`, `skip_labeling`), which is how you validate the plumbing with pre-computed data before a single solve is launched.

Phases talk through **queues**. The query strategy puts samples on a queue; the label strategy takes them off, labels them, and puts the results on a second queue; the driver drains that into the training pool. The queue type is a protocol too (`AbstractQueue`), so a multiprocessing or distributed queue drops in.

## Protocols, not base classes

Every pluggable piece is a `typing.Protocol`: as long as an object has the right methods and attributes at run time, it fits - no inheritance required. That is what lets this application define its label strategy lazily against physicsnemo without importing it at module scope.

| Protocol | You provide | Called |
|---|---|---|
| `QueryStrategy` | `max_samples`, `sample(query_queue)` - reads `driver.unlabeled_pool`, scores candidates, enqueues the chosen ones | once per step, phase 3 |
| `LabelStrategy` | `label(queue_to_label, serialize_queue)`, plus `__is_external_process__` (a solver is involved) and `__provides_fields__` (the field names it adds) | once per step, phase 4 |
| `MetrologyStrategy` | any measurement of progress beyond the validation loss | once per step, phase 2, optional |
| `TrainingProtocol`, `ValidationProtocol`, `InferenceProtocol` | one training step, one validation step, one inference step | inside the `TrainingLoop` |
| `TrainingLoop` | the epoch loop; `DefaultTrainingLoop` is the ready-made one | phase 1 |

All of them share `ActiveLearningProtocol`: a `__protocol_name__`, an `attach(driver)` that gives the strategy the driver's scope (`driver.learner`, the pools, the configs, the checkpoint directories), a configured `logger`, and `checkpoint_dir`/`strategy_dir` for anything it wants to persist.

## The Driver and its configs

`Driver` is the orchestrator. It is configured by dataclasses:

| Config | Holds |
|---|---|
| `DriverConfig` | infrastructure - batch size, logging, distributed settings, `max_steps`, `checkpoint_interval`, the skip flags |
| `StrategiesConfig` | the strategy *instances*: query, label, metrology, training loop |
| `TrainingConfig` | the training components - datapools, optimizer, scheduler, epochs per step |
| `OptimizerConfig` | the optimizer's own parameters |

It wraps the learner in `DistributedDataParallel` when ranks exist, checkpoints configurations, weights, optimizer state and queue contents at `checkpoint_interval`, injects the step index and phase into every log line, and resumes from `load_checkpoint(checkpoint_dir)`. Artifacts of a run: `.mdlus` checkpoints per step, `driver_log.json` with the timeline, and whatever the metrology strategy writes.

One contract that is easy to miss: **`TrainingConfig` needs a `train_datapool` even when training is skipped**, because labeled samples are drained into it after phase 4.

## Query strategies worth knowing

| Strategy | Signal | Cost | When |
|---|---|---|---|
| random | none | none | the baseline everything must beat; competitive early on |
| ensemble disagreement | variance across K trained models | K trainings | when you already train an ensemble for error bars |
| predictive entropy | the model's own output distribution | one pass | classifiers, or a diffusion ensemble's spread |
| physics residual | how badly the prediction violates the governing equations | one residual assembly per candidate | when a solver exists - this is the one only a solver-backed setup can offer |
| hybrid | e.g. 60 % uncertainty, 40 % random | | guards against the training set collapsing onto one region |

The honest finding from this application's own Examples case: on smooth one- and two-parameter thermal families, two or three solves saturate the surrogate and active learning is indistinguishable from random. The machinery pays off when the family is rough or high-dimensional, not before.

## Where Kratos sits

Kratos is the **label strategy**. `active_learning.kratos_label_strategy` provides `CreateKratosLabelStrategy`, which turns each queried sample into a real solve through an execution backend:

- the **in-process** backend runs the `AnalysisStage` in the current interpreter (small problems, notebooks);
- the **subprocess** backend (recommended) launches one Kratos process per case, keeps Kratos's MPI ranks and torch's distributed ranks in separate OS processes, fans out over `max_parallel_jobs`, and becomes an HPC job submission with an `srun`/`sbatch --wait` prefix. Results come back as the `.npz` files the case's `dataset_export_process` writes.

The query side is `active_learning.query_strategies` (ensemble disagreement, entropy, and the solver-residual score assembled through `physics.solver_residuals`), and the metrology side reuses the validation-metrics machinery. A sample is a `KratosALSample`, serializable to and from the queue.

<p align="center">
    <img src="../General/images/lifecycle.svg" alt="The surrogate lifecycle, of which active learning is the loop closing step 5 back to step 1"/>
</p>
<p align="center">Figure 2: Active learning is the five-step lifecycle with step 5 wired back to step 1, and the choice of the next solves made by the model.</p>

Upstream's own worked example is a two-moons classifier (`ClassifierUQQuery`, a pass-through `DummyLabelStrategy`, an F1 metrology); the pattern for an external CFD solver in the user guide is exactly the subprocess backend here, with JSON in and out.

Next: [Training utilities and performance](Training_Utilities_And_Performance.html).
