---
type: Python Method
title: Initialize
resource: python_scripts/analysis_stage.py#L70-L115
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_CreateModelers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_ModelersSetupGeometryModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_ModelersPrepareGeometryModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_ModelersSetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_PrepareModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_AddDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/ModifyInitialProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/ModifyInitialGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/__CreateListOfProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_InitializeInternals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/ModifyAfterSolverInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/ModifyBeforeSolutionLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/_GetSimulationName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def Initialize(self):`

# Calls

- [_CreateModelers](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_CreateModelers.md)
- [_ModelersSetupGeometryModel](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_ModelersSetupGeometryModel.md)
- [_ModelersPrepareGeometryModel](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_ModelersPrepareGeometryModel.md)
- [_ModelersSetupModelPart](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_ModelersSetupModelPart.md)
- [_PrepareModelPart](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_PrepareModelPart.md)
- [_AddDofs](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_AddDofs.md)
- [ModifyInitialProperties](../../../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyInitialProperties.md)
- [ModifyInitialGeometry](../../../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyInitialGeometry.md)
- [__CreateListOfProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/__CreateListOfProcesses.md)
- [_GetListOfProcesses](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetListOfProcesses.md)
- [_InitializeInternals](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_InitializeInternals.md)
- [ModifyAfterSolverInitialize](../../../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyAfterSolverInitialize.md)
- [ModifyBeforeSolutionLoop](../../../../functions/python_scripts/analysis_stage/AnalysisStage/ModifyBeforeSolutionLoop.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)
- [_GetSimulationName](../../../../functions/python_scripts/analysis_stage/AnalysisStage/_GetSimulationName.md)