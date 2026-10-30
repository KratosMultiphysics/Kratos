---
title: Active Learning
keywords: active learning label strategy driver
tags: [Active_Learning.md]
sidebar: physicsnemo_application
summary: Kratos as the ground-truth solver inside a physicsnemo.active_learning loop - the sample type, the in-process and subprocess execution backends, wiring into a Driver, query strategies including the solver residual, and metrology.
---

# Kratos as the ground-truth solver in an active-learning loop

`physicsnemo.active_learning` orchestrates train → metrology → query → label cycles through its `Driver`, and defines a `LabelStrategy` protocol explicitly designed to wrap an expensive external solver. This application implements that protocol so **Kratos is the labeling oracle**: the query strategy proposes design points, and each one is labeled by an actual Kratos solve.

<div class="mermaid">
flowchart LR
    query[query strategy picks cases] -->|query queue| label[CreateKratosLabelStrategy]
    subgraph backends [execution backends]
        direction TB
        inproc[InProcessBackend, this interpreter]
        sub[SubprocessBackend, one Kratos process per case]
        hpc[SubprocessBackend with an srun or sbatch prefix]
    end
    label --> backends
    backends -->|dataset_export_process writes .npz| harvest[labeled KratosALSample]
    harvest -->|serialize queue| driver[Driver drains into the training pool]
    driver --> train[fine-tune] --> metro[metrology] --> query
</div>
<p align="center">Figure 1: The seat Kratos takes in the loop. The framework itself is explained in <a href="../PhysicsNeMo_Basics/Active_Learning_Concepts.html">Active learning</a> under PhysicsNeMo Basics.</p>

## The sample type

`KratosALSample` is the queue item flowing through the driver's queues:

- `parameters`: the design point, as JSON-path → value overrides applied to the case's `ProjectParameters.json` (e.g. `"processes/loads_process_list/0/Parameters/modulus": 42.0`),
- `fields`: the labeled result, keys `<VARIABLE>__<location>` (the `DatasetExportProcess` npz convention),
- `sample_id`, `metadata`.

## Execution backends

**SubprocessBackend (recommended).** Copies a case template directory per sample, patches its `ProjectParameters.json`, runs `run_command` in the isolated case directory with retries and a timeout, and harvests results from the `.npz` files written by the case's own `DatasetExportProcess` (the template must list it among its processes). Running externally keeps Kratos MPI ranks and `torch.distributed` ranks in **separate OS processes** — and a `["srun", ...]` or `["sbatch", "--wait", ...]` prefix turns labeling into HPC job submission with no extra code.

**InProcessBackend.** Instantiates and runs an `AnalysisStage` in the current interpreter (resolving the stage class from `analysis_stage_module`, like CoSimulation's Kratos wrapper) and extracts output fields directly from the in-memory `Model`. Convenient for small problems; shares the GIL with training.

### Parallel labeling

`SubprocessBackend` fans a labeling batch out over concurrent subprocesses with the `"max_parallel_jobs"` setting (default 1 = serial). Every sample already runs in its own isolated case directory, so parallel runs cannot collide; results come back in submission order, and failures are reported per sample without aborting the batch (the `RunCases` contract on the backend base class). The label strategy drains the whole queue into one batch before calling the backend, so a `max_parallel_jobs` of N labels N cases at a time — with an `srun` prefix in `run_command`, that is N simultaneous HPC jobs.


<p align="center">
    <img src="images/labels.png" alt="Three rows of nodes, one per labelled case, coloured by the PRESSURE field each Kratos run produced"/>
</p>
<p align="center">Figure 2: The labels Kratos produced in notebook 04 - one row of nodes per queried case.</p>

## Wiring it into a Driver

```python
import queue
import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.kratos_label_strategy import CreateKratosLabelStrategy
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.subprocess_backend import SubprocessBackend
from physicsnemo.active_learning.driver import Driver
from physicsnemo.active_learning.config import DriverConfig, StrategiesConfig, TrainingConfig

backend = SubprocessBackend(Kratos.Parameters("""{
    "template_directory" : "my_case_template",
    "run_command"        : ["python3", "MainKratos.py"]
}"""))
label_strategy = CreateKratosLabelStrategy(backend, provides_fields={"VELOCITY__node_historical"})

driver = Driver(
    config=DriverConfig(batch_size=8),
    learner=my_model,
    strategies_config=StrategiesConfig(
        query_strategies=[my_query_strategy],
        queue_cls=queue.Queue,
        label_strategy=label_strategy),
    training_config=TrainingConfig(train_datapool=my_pool, max_training_epochs=50))
driver.run()
```

Notes:

- The driver appends labeled samples to `training_config.train_datapool` after the labeling phase (a plain `list` satisfies the `DataPool` protocol).
- A `train_datapool` is required whenever labeling is enabled, even with `skip_training=True`.
- Failed solves are logged and skipped (`label_strategy.failed_samples` counts them); one diverged case does not abort the labeling phase.

## Query strategies

`active_learning.query_strategies` provides ready-made implementations of the `QueryStrategy` protocol. Each draws a candidate pool from a user callable (`candidate_sampler(n)` → parameter dicts in `KratosALSample.parameters` format), scores it, and enqueues the top `max_samples` candidates:

| Factory | Score | Extra inputs |
|---|---|---|
| `CreateEnsembleDisagreementStrategy` | prediction variance across an ensemble of checkpoints (`"ensemble_checkpoints"`: list of `LoadModel` settings, ≥ 2) | `encode_candidate(dict)` → model input |
| `CreateEntropyStrategy` | Gaussian entropy proxy over `"num_stochastic_passes"` forward passes of one model with dropout kept active | `encode_candidate` |
| `CreateSolverResidualStrategy` | user-supplied `residual_evaluator(dict)` → float — typically the PDE residual of the surrogate's prediction (below) | `residual_evaluator` |

Common settings: `"max_samples"` and `"candidate_pool_size"`. `SelectTopCandidates(scores, candidates, max_samples)` is exposed for custom strategies (pure numpy).

## Solver residuals: the physics scores the surrogate

`solver_residuals.BuildResidualEvaluator(model_part)` assembles the **actual PDE residual** of whatever state the model part currently holds, through the same scheme/builder-and-solver machinery a solve would use (`BuildRHS` with fixed-DOF rows zeroed). Write the surrogate's prediction into the unknown variable, then:

```python
from KratosMultiphysics.PhysicsNeMoApplication.physics import solver_residuals
evaluator = solver_residuals.BuildResidualEvaluator(model_part)  # after solver Initialize / a solve
norm = evaluator.ComputeResidualNorm()          # scalar score
per_dof = evaluator.ComputeNodalResiduals()     # {(node_id, variable): |r|}
```

The residual is assembled in C++ **outside any autodiff graph: it is not differentiable** and can never be a gradient-carrying loss term — use it for query scoring, monitoring (see the Training page's epoch callbacks) and validation.

## Metrology

`active_learning.metrology.CreateValidationMetricsMetrology(settings, evaluation_callable)` implements the `MetrologyStrategy` protocol on the same metric backend as `ValidationMetricsProcess`: per iteration, the user callable returns named `(predicted, reference)` array pairs (e.g. surrogate predictions against `KratosALSample.fields` references), the configured `"metrics"` are recorded, and records serialize to `"output_file"` as JSON. Pass it via `StrategiesConfig(metrology_strategies=[...])` and leave `skip_metrology=False`.
