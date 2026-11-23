---
type: Python Module
title: analysis_stage
resource: python_scripts/analysis_stage.py#L1-L465
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-process-factory
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-model-parameters-factory
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-python-null-solver
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-modeler-factory
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [AnalysisStage](../../classes/python_scripts/analysis_stage/AnalysisStage.md)
- [__init__](../../functions/python_scripts/analysis_stage/AnalysisStage/__init__.md)
- [Run](../../functions/python_scripts/analysis_stage/AnalysisStage/Run.md)
- [KeepAdvancingSolutionLoop](../../functions/python_scripts/analysis_stage/AnalysisStage/KeepAdvancingSolutionLoop.md)
- [RunSolutionLoop](../../functions/python_scripts/analysis_stage/AnalysisStage/RunSolutionLoop.md)
- [Initialize](../../functions/python_scripts/analysis_stage/AnalysisStage/Initialize.md)
- [Finalize](../../functions/python_scripts/analysis_stage/AnalysisStage/Finalize.md)
- [GetComputingModelPart](../../functions/python_scripts/analysis_stage/AnalysisStage/GetComputingModelPart.md)
- [GetFinalData](../../functions/python_scripts/analysis_stage/AnalysisStage/GetFinalData.md)
- [InitializeSolutionStep](../../functions/python_scripts/analysis_stage/AnalysisStage/InitializeSolutionStep.md)
- [PrintAnalysisStageProgressInformation](../../functions/python_scripts/analysis_stage/AnalysisStage/PrintAnalysisStageProgressInformation.md)
- [SolveSolutionStep](../../functions/python_scripts/analysis_stage/AnalysisStage/SolveSolutionStep.md)
- [FinalizeSolutionStep](../../functions/python_scripts/analysis_stage/AnalysisStage/FinalizeSolutionStep.md)
- [OutputSolutionStep](../../functions/python_scripts/analysis_stage/AnalysisStage/OutputSolutionStep.md)
- [Check](../../functions/python_scripts/analysis_stage/AnalysisStage/Check.md)
- [Clear](../../functions/python_scripts/analysis_stage/AnalysisStage/Clear.md)
- [ModifyInitialProperties](../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyInitialProperties.md)
- [ModifyInitialGeometry](../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyInitialGeometry.md)
- [ModifyBeforeSolutionLoop](../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyBeforeSolutionLoop.md)
- [ModifyAfterSolverInitialize](../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyAfterSolverInitialize.md)
- [ApplyBoundaryConditions](../../functions/python_scripts/analysis_stage/AnalysisStage/ApplyBoundaryConditions.md)
- [ChangeMaterialProperties](../../functions/python_scripts/analysis_stage/AnalysisStage/ChangeMaterialProperties.md)
- [Save](../../functions/python_scripts/analysis_stage/AnalysisStage/Save.md)
- [Load](../../functions/python_scripts/analysis_stage/AnalysisStage/Load.md)
- [_GetSolver](../../functions/python_scripts/analysis_stage/AnalysisStage/_GetSolver.md)
- [_CreateSolver](../../functions/python_scripts/analysis_stage/AnalysisStage/_CreateSolver.md)
- [_HasSolver](../../functions/python_scripts/analysis_stage/AnalysisStage/_HasSolver.md)
- [_AdvanceTime](../../functions/python_scripts/analysis_stage/AnalysisStage/_AdvanceTime.md)
- [_AddVariables](../../functions/python_scripts/analysis_stage/AnalysisStage/_AddVariables.md)
- [_PrepareModelPart](../../functions/python_scripts/analysis_stage/AnalysisStage/_PrepareModelPart.md)
- [_AddDofs](../../functions/python_scripts/analysis_stage/AnalysisStage/_AddDofs.md)
- [_InitializeInternals](../../functions/python_scripts/analysis_stage/AnalysisStage/_InitializeInternals.md)
- [_ModelersSetupGeometryModel](../../functions/python_scripts/analysis_stage/AnalysisStage/_ModelersSetupGeometryModel.md)
- [_ModelersPrepareGeometryModel](../../functions/python_scripts/analysis_stage/AnalysisStage/_ModelersPrepareGeometryModel.md)
- [_ModelersSetupModelPart](../../functions/python_scripts/analysis_stage/AnalysisStage/_ModelersSetupModelPart.md)
- [_GetListOfModelers](../../functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfModelers.md)
- [_CreateModelers](../../functions/python_scripts/analysis_stage/AnalysisStage/_CreateModelers.md)
- [__BackwardCompatibleModelersCreation](../../functions/python_scripts/analysis_stage/AnalysisStage/__BackwardCompatibleModelersCreation.md)
- [_GetListOfProcesses](../../functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfProcesses.md)
- [_GetListOfOutputProcesses](../../functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfOutputProcesses.md)
- [_CreateProcesses](../../functions/python_scripts/analysis_stage/AnalysisStage/_CreateProcesses.md)
- [_GetOrderOfProcessesInitialization](../../functions/python_scripts/analysis_stage/AnalysisStage/_GetOrderOfProcessesInitialization.md)
- [_GetOrderOfOutputProcessesInitialization](../../functions/python_scripts/analysis_stage/AnalysisStage/_GetOrderOfOutputProcessesInitialization.md)
- [_CheckDeprecatedOutputProcesses](../../functions/python_scripts/analysis_stage/AnalysisStage/_CheckDeprecatedOutputProcesses.md)
- [_GetSimulationName](../../functions/python_scripts/analysis_stage/AnalysisStage/_GetSimulationName.md)
- [__CreateListOfProcesses](../../functions/python_scripts/analysis_stage/AnalysisStage/__CreateListOfProcesses.md)

# Imports

- `KratosMultiphysics`
- `KratosMultiphysics.process_factory`
- `KratosMultiphysics.kratos_utilities`
- `KratosMultiphysics.model_parameters_factory`
- `KratosMultiphysics.python_null_solver`
- `KratosMultiphysics.modeler_factory`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)