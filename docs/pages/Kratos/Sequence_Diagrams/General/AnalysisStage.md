---
keywords: AnalysisStage
tags: []
sidebar: kratos_sequence_diagrams
title: Analysis Stage
summary: 
---

<div class="mermaid">
sequenceDiagram
    Orchestrator->>+AnalysisStage: __init__
    AnalysisStage->>+Solver: CreateSolver
    AnalysisStage->>Solver: AddVariables
    AnalysisStage-->>Orchestrator: A.Stage Created

    Orchestrator->>AnalysisStage: Run

    AnalysisStage->>+AnalysisStage: Initialize
    AnalysisStage->>AnalysisStage: CreateModelers
    AnalysisStage->>AnalysisStage: _ModelersSetupGeometryModel
    AnalysisStage->>AnalysisStage: _ModelersPrepareGeometryModel
    AnalysisStage->>AnalysisStage: _ModelersSetupModelPart
    AnalysisStage->>Solver: ImportModelPart
    AnalysisStage->>Solver: PrepareModel
    AnalysisStage->>Solver: AddDofs
    AnalysisStage->>AnalysisStage: ModifyInitialProperties
    AnalysisStage->>AnalysisStage: ModifyInitialGeometry
    AnalysisStage->>AnalysisStage: _CreateListOfProcesses

    loop EveryProcess
      AnalysisStage->>Process: ExecuteInitialize
    end

    AnalysisStage->>Solver: Initialize

    loop EveryProcess
      AnalysisStage->>Process: ExecuteBeforeSolutionLoop
    end

    AnalysisStage-->>-AnalysisStage: A.Stage Initialized

    loop EvaluateStopCriteria
      AnalysisStage->>AnalysisStage: InitializeSolutionStep
      AnalysisStage->>Solver: Predict
      AnalysisStage->>Solver: SolveSolutionStep
      AnalysisStage->>AnalysisStage: FinalizeSolutionStep
    end

    Solver-->>-AnalysisStage: Destroy
    AnalysisStage-->>-Orchestrator: Run Finalized
</div>