---
title: Overview
keywords: physicsnemo machine learning surrogate
tags: [Overview.md]
sidebar: physicsnemo_application
summary: The bridge between Kratos Multiphysics and NVIDIA PhysicsNeMo - what it covers, its dependency policy, and where everything is documented.
---

# PhysicsNeMo Application

The *PhysicsNeMo Application* bridges *Kratos Multiphysics* with [NVIDIA PhysicsNeMo](https://docs.nvidia.com/physicsnemo/latest/index.html), a PyTorch-based framework for physics-ML. It covers the full loop between a finite-element solver and a machine-learning workflow:

<p align="center">
    <img src="images/architecture.svg" alt="Kratos on the left, PhysicsNeMo on the right, the application's packages in between, artifacts along the bottom"/>
</p>
<p align="center">Figure 1: The application sits between the two frameworks; <a href="Architecture.html">Architecture</a> explains the columns.</p>

- **Data out**: export simulation fields as training datasets ([Tensor Bridge](../Tensor_Bridge/Tensor_Bridge.html)).
- **Geometry out**: convert arbitrary Kratos meshes into PhysicsNeMo's simplicial mesh representation, with provenance to map predictions back ([Mesh Bridge](../Mesh_Bridge/Mesh_Bridge.html)).
- **Predictions in**: run trained models inside a Kratos solution loop, including warm-starting solves ([Inference](../Inference/Inference.html)).
- **Kratos as oracle**: use Kratos as the ground-truth solver inside a `physicsnemo.active_learning` loop ([Active Learning](../Active_Learning/Active_Learning.html)).

Built on those four, the rest of the application covers:

- **Training** ([Training](../Training/Training.html)) — datasets, schedules, streaming and warm restarts; **Reduced-order models** ([ROM](../Reduced_Order_Models/Reduced_Order_Models.html)); **Sequence models** ([Sequence Models](../Sequence_Models/Sequence_Models.html)); **Super-resolution** ([Super-Resolution](../Super_Resolution/Super_Resolution.html)).
- **Model families**: graph networks ([GNNs](../Graph_Neural_Networks/Graph_Neural_Networks.html)), point clouds and transformers ([Point Clouds](../Point_Clouds/Point_Clouds.html)), diffusion and downscaling ([Diffusion](../Diffusion/Diffusion.html)), Lagrangian particles ([Particle Methods](../Particle_Methods/Particle_Methods.html)), and external-aero datapipes ([CAE Datapipes](../CAE_Datapipes/CAE_Datapipes.html)).
- **Physics in the loss** ([Physics-Informed](../Physics_Informed/Physics_Informed.html)) — strong-form residual terms, PINN solves, the exact discrete residual through the real FEM assembly, and adjoint sensitivities cross-validated against Kratos's own adjoint stack.
- **Trust and scale**: prediction uncertainty and governance ([Uncertainty](../Uncertainty/Uncertainty.html)), MPI-distributed training and inference ([Distributed](../Distributed/Distributed.html)), and deploying a surrogate as a co-simulation solver ([CoSimulation](../CoSimulation/CoSimulation.html)).

Serving (ONNX export and Triton) is covered in [Inference](../Inference/Inference.html); adaptive remeshing in [Mesh Bridge](../Mesh_Bridge/Mesh_Bridge.html); streaming datasets and warm restarts in [Training](../Training/Training.html).

## Documentation map

| Section | Read it for |
|---|---|
| **General** | [Architecture](Architecture.html), [Kratos concepts for ML readers](Kratos_Concepts_For_ML.html), [Where things live](Module_Map.html), [From scratch](From_Scratch.html), the [Process reference](Process_Reference.html), [Installation and environment](Installation_And_Environment.html), [Performance](Performance.html), [Testing and contributing](Testing_And_Contributing.html), [Troubleshooting and traps](Troubleshooting_And_Traps.html), a [Glossary](Glossary.html), and the [Retrospectives](Retrospectives.html) of how the roadmap moved |
| **PhysicsNeMo Basics** | NVIDIA PhysicsNeMo itself, in fourteen pages from [what it is](../PhysicsNeMo_Basics/Overview.html) to [versions and compatibility](../PhysicsNeMo_Basics/Versions_And_Compatibility.html) |
| **Tensor Bridge**, **Mesh Bridge** | data and geometry out of Kratos, and predictions back |
| **Training**, **Inference** | the loop, the checkpoints and cards, every way to run a model inside a solve |
| **Active Learning** | Kratos as the labeling oracle |
| **Super-Resolution**, **Graph Neural Networks**, **Particle Methods**, **Sequence Models**, **Point Clouds**, **CAE Datapipes**, **Reduced Order Models**, **Diffusion** | one page per data shape and model family |
| **Adjoint**, **Physics-Informed** | gradients and residuals through the real solver |
| **Uncertainty**, **CoSimulation**, **Distributed** | trust, coupling, ranks |
| **Examples** | the gallery of twenty-six use cases and the twenty notebooks |

The README's roadmap lists what is pending, with its gate; nothing there is duplicated here.

## Dependency policy

`torch` and `nvidia-physicsnemo` are **optional, pure-Python runtime dependencies**:

- `import KratosMultiphysics.PhysicsNeMoApplication` always succeeds without them; only the specific submodules that need them fail — lazily, at call time, with an actionable error message.
- There is intentionally **no CMake gating** (no `find_package(Torch)`, no `USE_PHYSICSNEMO` option): nothing in the C++ core links against torch or CUDA. Do not add such a gate.
- Contributions must keep all `import torch` / `import physicsnemo` statements inside lazy `_TryImport*()` helpers.
- **Packaging**: the application ships like any per-application Kratos wheel (`scripts/wheels/build_wheel.py` derives app wheels automatically from what was compiled). `torch`/`nvidia-physicsnemo` are never wheel dependencies — they stay optional runtime `pip install`s, and the import-contract test guarantees a wheel installed without them remains fully importable.

## Working with other Kratos applications

Some tests and examples run **real Kratos solves** (e.g. a stationary heat-conduction case through `ConvectionDiffusionApplication`). They are conditioned on the compiled applications following the standard pattern (as used by CoSimulationApplication):

```python
import KratosMultiphysics.kratos_utilities as kratos_utils
have_convection_diffusion = kratos_utils.CheckIfApplicationsAvailable(
    "ConvectionDiffusionApplication", "LinearSolversApplication")

@KratosUnittest.skipUnless(have_convection_diffusion, "Missing required applications: ...")
class TestRealSolver(...): ...
```

Tests are always registered and skip cleanly (with an explanatory message) on builds that lack the required applications. The reusable in-memory thermal case powering these tests lives in `tests/kratos_solver_cases/thermal_case.py` (meshed with the core `StructuredMeshGeneratorProcess`, no `.mdpa` fixtures, solver fed through `"input_type": "use_input_model_part"`).

## Examples

Twenty runnable notebooks live in `applications/PhysicsNeMoApplication/examples/notebooks/`, and twenty-six documented use cases against real solves in the Examples repository; both are indexed on the [Examples](../Examples/Examples.html) page.

## Installation

Compile like any other Kratos application:

```bash
add_app ${KRATOS_APP_DIR}/PhysicsNeMoApplication
```

To use the ML features additionally install the runtime dependencies:

```bash
pip install torch             # CPU or CUDA build, your choice
pip install nvidia-physicsnemo
pip install nvidia-physicsnemo-cfd   # optional: hybrid-initialization recipes + CFD metric registry
```
