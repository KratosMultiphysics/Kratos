# OptimizationApplication

The Kratos OptimizationApplication is a framework for solving optimization problems in continuum mechanics. It is supposed to handle both gradient-based (adjoint-based) and gradient-free methods.

## Main Features

- State-of-the-art techniques and algorithms for shape, thickness and material/topology optimization.
- Efficient and consistent filtering techniques for parametrization-free shape, thickness and material/topology optimization.
- Abstract problem formulation which enables concurrent and nested multilevel-multi-scale optimization problems.
- Adaptive gradient-projection technique, developed specially for problems with an arbitrary large number of design variables of different scales.
- Modular implementation which enables analysis and optimization of multi-physics problems. sdfsdfsd
- Realization and implementation of additive manufacturing constraints, e.g. hangover conditions (support structures), stackability and geometric limitations.

## Understand the structure and flow of the application

- Take a look at the tests
- Work from the json to the main python file
- Examine each level of classes step by step

Diagram with the general structure and some comments

```Mermaid
%%{init: {'theme':'neutral',"sequence": { "wrap": true}}}%%

sequenceDiagram
	Optimizer(AnalysisStage)->>OptimizationSolver(PythonSolver): for every iteration:<br> solve time step / solve the next iteration step
	OptimizationSolver(PythonSolver)->>Algorithm: for all algorithms: <br> optimize the current parameters
```

```Mermaid
%%{init: {'theme':'neutral',"sequence": { "wrap": true}}}%%
sequenceDiagram
	Algorithm->>Objectives: for all Objectives:<br> what is your current value and your sensitivities
        Objectives->>ExecutionPolicyWrapper: passes request on
        Algorithm->>Constrains: for all Constrains:<br> are you active and if: what is your current value and your sensitivities
        Constrains->>ExecutionPolicyWrapper: passes request on

        Algorithm->>ControlWrapper: Modify the sensitivities
        ControlWrapper->>Modifier: for all modifier:<br>Modify the sensitivities

        Algorithm->>Algorithm: calculate next parameter update

        Algorithm->>ControlWrapper: Modify the parameter update
        ControlWrapper->>Modifier: for all modifier:<br>Modify the parameter update

        Algorithm->>ControlWrapper: apply update to model
        ControlWrapper->>Control: passes request on

        Algorithm->>ControlWrapper: Modify the final updated parameters
        ControlWrapper->>Modifier: for all modifier:<br>Modify the final updated parameters
```

## Understand the JSON config of the application

```Mermaid
%%{init: {'theme':'neutral'}}%%
flowchart LR
        JSON --> ProblemData
        JSON --> OptimizationSettings

        OptimizationSettings --> Meshes
        subgraph meshes
        Meshes -- Array_of --> D1[Dict]
        D1 --> M1[Module]
        D1 --> T1[Type]
        D1 --> S1[Settings]
        end

        OptimizationSettings --> Analyses
        subgraph analyses
        Analyses -- Array_of --> D2[Dict]
        D2 --> ExecutionPolicySettings
        end

        OptimizationSettings --> Responses
        subgraph responses
        Responses -- Array_of --> D3[Dict]
        D3 --> N3[Name]
        D3 --> M3[Module]
        D3 --> T3[Type]
        D3 --> O3[Objective]
        D3 --> C3[Constrains]
        D3 --> S3[Settings]
        end

        OptimizationSettings --> Controls
        subgraph controls
        Controls -- Array_of --> D4[Dict]
        D4 --> N4[Name]
        D4 --> T4[Type]
        D4 --> S4[Settings]
        D4 --> M4[ModifiersList]
        end

        OptimizationSettings --> Algorithm
        subgraph algorithm
        Algorithm -- Array_of --> D5[Dict]
        D5 --> T5[Type]
        D5 --> S5[Settings]
        end
```
