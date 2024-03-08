---
keywords: BuilderAndSolver
tags: []
sidebar: kratos_sequence_diagrams
title: Builder And Solver
summary: 
---

## Run Sequence
<div class="mermaid">
sequenceDiagram
    autonumber
    Orchestrator->>+AnalysisStage: __init__
    create participant Solver
    AnalysisStage->>Solver: CreateSolver
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
    create participant Process
    AnalysisStage->>Process: _CreateProcesses(Process)
    create participant OutputProcess
    AnalysisStage->>OutputProcess: _CreateProcesses(OutputProcess)

    loop EveryProcess
      AnalysisStage->>Process: ExecuteInitialize
    end

    AnalysisStage->>Solver: InitializeSolver

    loop EveryProcess
      AnalysisStage->>Process: ExecuteBeforeSolutionLoop
    end

    AnalysisStage-->>-AnalysisStage: A.Stage Initialized

    AnalysisStage->>+AnalysisStage: RunSolutionLoop

    loop EvaluateStopCriteria

      AnalysisStage->>+AnalysisStage: InitializeSolutionSteep

      loop EveryProcess
        AnalysisStage->>Process: ExecuteInitializeSolutionStep
      end

      AnalysisStage-->>-AnalysisStage: Solution Step Initialized

      AnalysisStage->>Solver: Predict
      AnalysisStage->>Solver: SolveSolutionStep
      AnalysisStage->>AnalysisStage: CheckConvergence
      
      AnalysisStage->>+AnalysisStage: FinalizeSolutionStep

      loop EveryProcess
        AnalysisStage->>Process: ExecuteFinalizeSolutionStep
      end

      AnalysisStage-->>-AnalysisStage: Solution Step Finalize

      AnalysisStage->>+AnalysisStage: OutputSolutionStep

      loop EveryOutputProcess
        alt IsOutputStep
          loop EveryProcess
            AnalysisStage->>Process: ExecuteBeforeOutputStep
          end

          AnalysisStage->>OutputProcess: PrintOutput

          loop EveryProcess
            AnalysisStage->>Process: ExecuteAfterOutputStep
          end
        end
      end

      AnalysisStage-->>-AnalysisStage: Output Finalized
    
    end

    AnalysisStage-->>-AnalysisStage: Solution Loop Finalized
    AnalysisStage-->>-Orchestrator: Run Finalized
</div>